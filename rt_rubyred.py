#!/usr/bin/env python3
"""Q16Seq / RT-RubyRed: real-time RubyRed analysis of a MinKNOW output directory.

Watches a directory for new *.fastq.gz files during sequencing, demultiplexes
each new batch of files with torchlex (inner BRK barcodes), classifies it with
RubyRed, merges the per-batch feature tables and keeps an HTML report up to date.
Samples spiked with the lambda internal standard (lambda1-3) are converted to
absolute 16S rRNA gene copies (Q16Seq).

Author: Lukasz Krych, University of Copenhagen
"""
import argparse
import datetime as dt
import glob
import gzip
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

import pandas as pd
import psutil

TOOL_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, TOOL_DIR)
import report  # noqa: E402

VERSION = "1.1.0"
DEFAULT_CLASSIFIER = "{rr_home}/resources/classifiers/MIMt_16s_311025_species"
HOME = os.path.expanduser("~")

# ---------------------------------------------------------------- terminal
C = {
    "red": "\033[38;5;160m", "ruby": "\033[1;38;5;161m", "grey": "\033[38;5;245m",
    "green": "\033[38;5;71m", "yellow": "\033[38;5;178m", "blue": "\033[38;5;74m",
    "gold": "\033[1;38;5;214m",
    "bold": "\033[1m", "off": "\033[0m",
}
if not sys.stdout.isatty() and not os.environ.get("RT_COLOR"):
    C = {k: "" for k in C}

LOG_FH = None


def now():
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def say(msg, colour=None, plain=None):
    """Print a timestamped line to the terminal and the run log."""
    ts = dt.datetime.now().strftime("%H:%M:%S")
    col = C.get(colour, "") if colour else ""
    print(f"{C['grey']}[{ts}]{C['off']} {col}{msg}{C['off'] if col else ''}", flush=True)
    if LOG_FH:
        LOG_FH.write(f"[{now()}] {plain if plain is not None else msg}\n")
        LOG_FH.flush()


def banner():
    with open(os.path.join(TOOL_DIR, "assets", "banner.txt")) as fh:
        art = fh.read().rstrip("\n")
    with open(os.path.join(TOOL_DIR, "assets", "banner_q16seq.txt")) as fh:
        q16 = fh.read().rstrip("\n").split("\n")
    width = max(len(line) for line in art.split("\n"))
    pad = " " * max(0, (width - max(len(line) for line in q16)) // 2)
    print()
    print(C["ruby"] + art + C["off"])
    print(C["gold"] + "\n".join(pad + line for line in q16) + C["off"])
    print(f"{C['grey']}  real-time rRNA amplicon classification + Q16Seq λ quantification  ·  v{VERSION}  ·  "
          f"Lukasz Krych, University of Copenhagen{C['off']}\n")


def human_bytes(n):
    for unit in ("B", "K", "M", "G", "T"):
        if abs(n) < 1024 or unit == "T":
            return f"{n:.1f}{unit}" if unit != "B" else f"{n}B"
        n /= 1024.0


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---------------------------------------------------------------- resources
class Resources:
    """Keeps the tool within a CPU / RAM budget and reports usage."""

    def __init__(self, cpu_frac, mem_frac, threads=None):
        n_cpu = os.cpu_count() or 1
        self.threads = threads or max(1, int(n_cpu * cpu_frac))
        # Pin this process (and therefore every child) to the upper half of the cores:
        # a hard CPU cap that does not depend on cgroup cpu delegation.
        self.cpus = list(range(n_cpu - self.threads, n_cpu))
        try:
            os.sched_setaffinity(0, self.cpus)
        except (AttributeError, OSError):
            pass
        try:
            os.nice(10)
        except OSError:
            pass
        self.mem_total = psutil.virtual_memory().total
        self.mem_limit = int(self.mem_total * mem_frac)
        self.scope_ok = self._check_scope()
        self.proc = psutil.Process()
        psutil.cpu_percent(None, percpu=True)

    def _check_scope(self):
        if not shutil.which("systemd-run"):
            return False
        r = subprocess.run(["systemd-run", "--user", "--scope", "--quiet",
                            "-p", "MemoryMax=1G", "true"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return r.returncode == 0

    def wrap(self, cmd):
        """Prefix a command so it runs niced, low-IO-priority and memory-capped."""
        pre = []
        if self.scope_ok:
            pre = ["systemd-run", "--user", "--scope", "--quiet",
                   "-p", f"MemoryMax={self.mem_limit}", "-p", "MemorySwapMax=0"]
        if shutil.which("ionice"):
            pre += ["ionice", "-c2", "-n7"]
        return pre + cmd

    def usage(self):
        # CPU: load on the cores the tool is pinned to (short-lived qiime workers would
        # be missed by per-process sampling); expressed as % of the whole machine.
        n_cpu = os.cpu_count() or 1
        per = psutil.cpu_percent(None, percpu=True)
        pinned = [per[i] for i in self.cpus if i < len(per)]
        rss = 0
        for p in [self.proc] + self._children():
            try:
                rss += p.memory_info().rss
            except psutil.Error:
                pass
        vm = psutil.virtual_memory()
        return {
            "tool_cpu": sum(pinned) / n_cpu,
            "tool_rss": rss,
            "sys_cpu": sum(per) / len(per) if per else 0.0,
            "sys_mem": vm.percent,
        }

    def _children(self):
        try:
            return self.proc.children(recursive=True)
        except psutil.Error:
            return []

    def line(self):
        u = self.usage()
        return (f"CPU tool {u['tool_cpu']:4.1f}% (cap 50%) · sys {u['sys_cpu']:4.1f}%  |  "
                f"RAM tool {human_bytes(u['tool_rss'])} / cap {human_bytes(self.mem_limit)} · "
                f"sys {u['sys_mem']:.0f}%")


# ---------------------------------------------------------------- run state
class Run:
    def __init__(self, args):
        self.args = args
        self.watch = os.path.realpath(args.watch)
        self.date = dt.datetime.now().strftime("%Y%m%d")
        self.run_name = args.name or derive_run_name(self.watch)
        self.dir_name = os.path.basename(self.watch.rstrip("/"))
        if args.resume:
            self.out = os.path.realpath(args.resume)
            meta = self.load_meta()
            self.date = meta.get("date", self.date)
            self.run_name = meta.get("run_name", self.run_name)
        else:
            stamp = dt.datetime.now().strftime("%Y%m%d_%H%M")
            self.out = os.path.join(os.path.expanduser(args.outdir),
                                    f"Q16Seq_{self.run_name}_{self.dir_name}_{stamp}")
        self.demux = os.path.join(self.out, "fastq_demultiplexed")
        self.tables = os.path.join(self.out, "feature_tables")
        self.rr_dir = os.path.join(self.out, "rubyred_batches")
        self.logs = os.path.join(self.out, "logs")
        self.rep_seqs = os.path.join(self.out, "rep_seqs")
        self.index = os.path.join(self.out, "input_file_index.tsv")
        self.stats_file = os.path.join(self.out, "batch_stats.json")
        self.merged = os.path.join(
            self.out, f"{self.run_name}_{self.dir_name}_{self.date}_feature_table.tsv")
        self.report = os.path.join(
            self.out, f"{self.run_name}_{self.dir_name}_{self.date}_report.html")
        for d in (self.out, self.demux, self.tables, self.rr_dir, self.logs, self.rep_seqs):
            os.makedirs(d, exist_ok=True)
        if not os.path.exists(self.index):
            with open(self.index, "w") as fh:
                fh.write("file\tsize_bytes\tmtime\tbatch\tadded\n")
        self.stats = self.load_stats()
        self.save_meta()

    # ---- persistence
    def meta_path(self):
        return os.path.join(self.out, "run_meta.json")

    def load_meta(self):
        try:
            with open(self.meta_path()) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}

    def save_meta(self):
        meta = {"tool": "RT-RubyRed", "version": VERSION, "watch_dir": self.watch,
                "run_name": self.run_name, "dir_name": self.dir_name, "date": self.date,
                "output_dir": self.out, "started": self.load_meta().get("started", now())}
        with open(self.meta_path(), "w") as fh:
            json.dump(meta, fh, indent=2)
        return meta

    def load_stats(self):
        try:
            with open(self.stats_file) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {"batches": []}

    def save_stats(self):
        tmp = self.stats_file + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(self.stats, fh, indent=1)
        os.replace(tmp, self.stats_file)

    def indexed_files(self):
        with open(self.index) as fh:
            next(fh)
            return {line.split("\t", 1)[0] for line in fh if line.strip()}

    def append_index(self, files, batch):
        with open(self.index, "a") as fh:
            for f in files:
                st = os.stat(f)
                fh.write(f"{f}\t{st.st_size}\t{int(st.st_mtime)}\t{batch}\t{now()}\n")

    def next_batch_no(self):
        done = [b["batch"] for b in self.stats["batches"]]
        return (max(done) + 1) if done else 1


def derive_run_name(path):
    """MinKNOW layout: data/<experiment>/<sample>/<flowcell_run>/fastq_pass[/barcodeNN]."""
    parts = path.rstrip("/").split("/")
    if "fastq_pass" in parts:
        i = parts.index("fastq_pass")
        if i >= 3:
            return parts[i - 3]
    return parts[-2] if len(parts) > 1 else parts[-1]


# ---------------------------------------------------------------- child processes
CURRENT = {"proc": None}
STOP = {"flag": False, "count": 0}


def on_sigint(signum, frame):
    STOP["count"] += 1
    STOP["flag"] = True
    if STOP["count"] == 1:
        print(f"\n{C['yellow']}  Interrupt received → stopping current step and finalising "
              f"(merging BRK files across batches).{C['off']}", flush=True)
    else:
        print(f"{C['yellow']}  Already finalising, please wait...{C['off']}", flush=True)
    kill_current()


def kill_current():
    p = CURRENT["proc"]
    if p and p.poll() is None:
        try:
            os.killpg(p.pid, signal.SIGTERM)
        except OSError:
            pass


def run_child(cmd, log_path, res, label, status_every=30):
    """Run a command in its own process group, logging output; show periodic status."""
    with open(log_path, "a") as log:
        log.write(f"\n[{now()}] $ {' '.join(cmd)}\n")
        log.flush()
        p = subprocess.Popen(res.wrap(cmd), stdout=log, stderr=subprocess.STDOUT,
                             start_new_session=True)
        CURRENT["proc"] = p
        t0 = time.time()
        last = t0
        while p.poll() is None:
            time.sleep(1)
            if STOP["flag"]:
                kill_current()
                try:
                    p.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
                break
            if time.time() - last >= status_every:
                last = time.time()
                mins = (last - t0) / 60
                step = last_log_line(log_path)
                say(f"   … {label} {mins:4.1f} min  |  {res.line()}", "grey",
                    plain=f"   ... {label} {mins:.1f} min")
                if step:
                    print(f"{C['grey']}             ↳ {step[:110]}{C['off']}", flush=True)
        CURRENT["proc"] = None
        return p.wait()


def last_log_line(path):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, 2)
            fh.seek(max(0, fh.tell() - 4000))
            lines = fh.read().decode(errors="replace").replace("\r", "\n").splitlines()
        for ln in reversed(lines):
            ln = ln.strip()
            if ln and not ln.startswith(("perl:", "LC_", "LANG", "are supported")):
                return ln
    except OSError:
        pass
    return ""


# ---------------------------------------------------------------- watching
def scan_new(run):
    """Return settled FASTQ files in the watch dir that are not yet indexed."""
    pat = "**/*.fastq*" if run.args.recursive else "*.fastq*"
    seen = run.indexed_files()
    cutoff = time.time() - run.args.settle
    new = []
    for f in glob.glob(os.path.join(run.watch, pat), recursive=run.args.recursive):
        if not (f.endswith(".fastq.gz") or f.endswith(".fastq")):
            continue
        if f in seen:
            continue
        try:
            if os.path.getmtime(f) > cutoff:
                continue  # possibly still being written
        except OSError:
            continue
        new.append(f)
    return sorted(new, key=natural_key)


def count_fastq_reads(path):
    n = 0
    with open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n // 4


# ---------------------------------------------------------------- batch processing
def process_batch(run, res, files):
    b = run.next_batch_no()
    tag = f"batch{b:02d}"
    bdir = os.path.join(run.demux, tag)
    rr_out = os.path.join(run.rr_dir, tag)
    stage = os.path.join(run.rr_dir, f".{tag}_input")
    blog = os.path.join(run.logs, f"{tag}.log")
    for d in (bdir, rr_out, stage):
        shutil.rmtree(d, ignore_errors=True)
    t0 = time.time()
    size = sum(os.path.getsize(f) for f in files)
    say(f"{C['bold']}▶ {tag}{C['off']}{C['blue']}: {len(files)} new file(s) "
        f"({human_bytes(size)}) → torchlex demultiplexing", "blue",
        plain=f"> {tag}: {len(files)} new files ({human_bytes(size)}) -> torchlex")

    # 1. torchlex (inner BRK barcodes); same settings RubyRed uses internally
    lst = os.path.join(run.logs, f"{tag}_files.txt")
    with open(lst, "w") as fh:
        fh.write("\n".join(files) + "\n")
    cmd = [run.args.torchlex, "-i", *files, "-o", bdir, "-t", str(res.threads),
           "--barcodes", run.args.barcodes, "--min-qual", "0", "--min-len", "1",
           "--min-reads", "1", "--max-len", "10000"]
    rc = run_child(cmd, blog, res, "torchlex")
    if STOP["flag"] or rc != 0:
        if not STOP["flag"]:
            say(f"✗ {tag}: torchlex failed (exit {rc}); see {blog}. "
                f"Files stay unindexed and will be retried.", "red")
            time.sleep(run.args.poll)
        shutil.rmtree(bdir, ignore_errors=True)
        return None

    brk = sorted(glob.glob(os.path.join(bdir, "BRK*.fastq")), key=natural_key)
    raw = {os.path.basename(f)[:-6]: count_fastq_reads(f) for f in brk}
    n_raw = sum(raw.values())
    say(f"   torchlex: {n_raw:,} reads assigned to {len(raw)} BRK barcode(s)", "grey")

    # 2. RubyRed on the batch (hard-linked staging copy: RubyRed moves its inputs)
    table_out = os.path.join(run.tables, f"feature_table_{tag}.tsv")
    added = 0
    status = "ok"
    if brk:
        os.makedirs(stage, exist_ok=True)
        for f in brk:
            os.link(f, os.path.join(stage, os.path.basename(f)))
        say(f"   RubyRed classification ({res.threads} threads)…", "grey")
        cmd = ["bash", run.args.rubyred, "-i", stage, "-d", "-o", rr_out,
               "-p", str(res.threads), "-m", "1", "-t", str(run.args.subsample),
               "-f", "1", "-s", os.path.join(run.args.rr_home, "scripts"),
               "-w", os.path.join(run.args.rr_home, "resources"),
               *classifier_args(run.args), *run.args.rubyred_args.split()]
        rc = run_child(cmd, blog, res, "RubyRed")
        shutil.rmtree(stage, ignore_errors=True)
        shutil.rmtree(os.path.join(rr_out, "demux_data"), ignore_errors=True)
        if STOP["flag"]:
            shutil.rmtree(bdir, ignore_errors=True)
            shutil.rmtree(rr_out, ignore_errors=True)
            return None
        ft = os.path.join(rr_out, "final_feature_table", "feature_table.tsv")
        if rc == 0 and os.path.exists(ft):
            df = read_table(ft)
            df.to_csv(table_out, sep="\t", index_label="#OTU ID")
            added = int(df.values.sum())
            n_seq = merge_rep_seqs(run, rr_out, b)
            say(f"   rep_seqs: {n_seq:,} reads appended to rep_seqs/ (per-taxon FASTA)", "grey")
        else:
            status = "rubyred_failed"
            say(f"✗ {tag}: RubyRed failed (exit {rc}); reads are kept in {bdir} "
                f"but not added to the table. Log: {blog}", "red")
    else:
        status = "no_reads"

    # 3. bookkeeping: files are indexed only once the batch is complete
    run.append_index(files, b)
    run.stats["batches"].append({
        "batch": b, "tag": tag, "files": len(files), "bytes": size,
        "raw_reads": raw, "raw_total": n_raw, "classified_added": added,
        "status": status, "finished": now(), "minutes": round((time.time() - t0) / 60, 2),
        "table": os.path.basename(table_out) if status == "ok" else None,
    })
    run.save_stats()

    total = sum(x["classified_added"] for x in run.stats["batches"])
    say(f"{C['bold']}✔ {tag} done{C['off']}{C['green']} in {(time.time()-t0)/60:.1f} min: "
        f"{C['bold']}{added:,} reads added{C['off']}{C['green']} to the feature table "
        f"(run total {total:,})", "green",
        plain=f"OK {tag}: {added:,} reads added (total {total:,})")
    print(f"             input : {run.watch}")
    print(f"             output: {run.out}")
    if status == "ok":
        print(f"             table : feature_tables/{os.path.basename(table_out)}")
    merge_and_report(run, min_batches=2)
    return b


def merge_rep_seqs(run, rr_out, batch):
    """Append a batch's per-taxon reads to rep_seqs/<taxon>.fasta.gz (strains pooled).

    RubyRed numbers reads from 1 in every run, so IDs are made unique across batches:
    >BRK185_89 in batch02 becomes >BRK185_b02_89.
    """
    marker = os.path.join(run.rep_seqs, f".merged_batch{batch:02d}")
    if os.path.exists(marker):  # already appended (resumed run)
        return 0
    tag = f"_b{batch:02d}_"
    n = 0
    for src in sorted(glob.glob(os.path.join(rr_out, "rep_seqs", "*.fasta.gz"))):
        # species-level classifier adds a T__<strain> rank: pool strains per species
        dest = os.path.join(run.rep_seqs,
                            re.sub(r"_T__.*(?=\.fasta\.gz$)", "", os.path.basename(src)))
        with gzip.open(src, "rt") as fin, open(dest, "ab") as raw, \
                gzip.GzipFile(fileobj=raw, mode="wb") as fout:
            for line in fin:
                if line.startswith(">"):
                    sample, _, idx = line[1:].rstrip("\n").rpartition("_")
                    line = f">{sample}{tag}{idx}\n" if sample else line
                    n += 1
                fout.write(line.encode())
    open(marker, "w").close()
    return n


def classifier_args(args):
    """RubyRed -c/-r/-y for the chosen classifier folder (<dir>/<name>_classifier.qza ...)."""
    if args.classifier in ("", "default", "original"):
        return []
    d = os.path.realpath(os.path.expanduser(args.classifier))
    name = os.path.basename(d)
    files = {"-c": f"{name}_classifier.qza", "-r": f"{name}_seqs.qza", "-y": f"{name}_tax.qza"}
    out = []
    for flag, f in files.items():
        path = os.path.join(d, f)
        if os.path.exists(path):
            out += [flag, path]
    return out


def read_table(path):
    df = pd.read_csv(path, sep="\t", comment=None, skiprows=_skip_rows(path), index_col=0)
    df.index.name = "#OTU ID"
    return df.astype(float).round().astype(int)


def _skip_rows(path):
    with open(path) as fh:
        first = fh.readline()
    return 1 if first.startswith("# Constructed") else 0


def merge_and_report(run, min_batches=2):
    tables = sorted(glob.glob(os.path.join(run.tables, "feature_table_batch*.tsv")),
                    key=natural_key)
    if len(tables) < min_batches or not tables:
        return
    merged = None
    for t in tables:
        df = read_table(t)
        merged = df if merged is None else merged.add(df, fill_value=0)
    merged = merged.fillna(0).astype(int)
    merged = merged[sorted(merged.columns, key=natural_key)]
    merged = merged.loc[merged.sum(axis=1).sort_values(ascending=False).index]
    tmp = run.merged + ".tmp"
    merged.to_csv(tmp, sep="\t", index_label="#OTU ID")
    os.replace(tmp, run.merged)
    say(f"   merged {len(tables)} batch tables → {os.path.basename(run.merged)} "
        f"({merged.shape[0]} taxa × {merged.shape[1]} samples, "
        f"{int(merged.values.sum()):,} reads)", "grey")
    try:
        report.build(run.report, merged, run.stats, run.save_meta(), TOOL_DIR, lambda_cfg(run.args))
        say(f"   report updated → {os.path.basename(run.report)}", "grey")
    except Exception as exc:  # the report must never stop the analysis
        say(f"   report generation failed: {exc}", "red")
    write_quant(run.args, merged, run.merged)


def lambda_cfg(args):
    return {"l1": args.lambda1, "map": args.lambda_map_data}


def write_quant(args, table, merged_path):
    """<prefix>_quantitative_table.tsv + _lambda_qc.tsv when --lambda1 and --lambda-map are given."""
    lam = report.lambda_rows(table)
    if not lam:
        say("   ⚠ no λ internal-standard reads in the table (S__lambda_standard)", "yellow")
        return
    if not (args.lambda1 and args.lambda_map_data):
        return
    try:
        copies, qc = report.quantify(table, args.lambda1, args.lambda_map_data)
        base = merged_path[:-len("_feature_table.tsv")] if merged_path.endswith("_feature_table.tsv") \
            else os.path.splitext(merged_path)[0]
        copies.to_csv(base + "_quantitative_table.tsv", sep="\t", index_label="#OTU ID",
                      float_format="%.4g")
        qc.to_csv(base + "_lambda_qc.tsv", sep="\t")
        n_ok = int((qc["status"] == "ok").sum())
        say(f"   Q16Seq: copies for {copies.shape[1]} sample(s) → "
            f"{os.path.basename(base)}_quantitative_table.tsv ({n_ok} within ROQ)", "grey")
    except Exception as exc:
        say(f"   quantification failed: {exc}", "red")


# ---------------------------------------------------------------- finalisation
def finalise(run):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    say("■ Finalising: concatenating BRK files across batches", "yellow")
    batches = sorted(glob.glob(os.path.join(run.demux, "batch*")), key=natural_key)
    done_tags = {b["tag"] for b in run.stats["batches"]}
    for bd in batches:  # drop any half-finished batch (its files were never indexed)
        if os.path.basename(bd) not in done_tags:
            shutil.rmtree(bd, ignore_errors=True)
    batches = [b for b in batches if os.path.basename(b) in done_tags]
    names = sorted({os.path.basename(f) for bd in batches
                    for f in glob.glob(os.path.join(bd, "BRK*.fastq"))}, key=natural_key)
    for name in names:
        dest = os.path.join(run.demux, name)
        tmp = dest + ".partial"
        with open(tmp, "wb") as out:
            if os.path.exists(dest):  # resumed run that was finalised before
                with open(dest, "rb") as src:
                    shutil.copyfileobj(src, out, 16 << 20)
            for bd in batches:
                src_path = os.path.join(bd, name)
                if os.path.exists(src_path):
                    with open(src_path, "rb") as src:
                        shutil.copyfileobj(src, out, 16 << 20)
        os.replace(tmp, dest)
    for bd in batches:
        shutil.rmtree(bd, ignore_errors=True)
    say(f"   {len(names)} BRK file(s) written to {run.demux}; {len(batches)} batch "
        f"folder(s) cleared", "grey")
    merge_and_report(run, min_batches=1)
    s = run.stats["batches"]
    say(f"{C['bold']}RT-RubyRed finished{C['off']}{C['green']}: {len(s)} batch(es), "
        f"{sum(x['files'] for x in s)} input files, "
        f"{sum(x['raw_total'] for x in s):,} demultiplexed reads, "
        f"{sum(x['classified_added'] for x in s):,} classified reads", "green")
    print(f"             output: {run.out}")
    if os.path.exists(run.merged):
        print(f"             table : {os.path.basename(run.merged)}")
        print(f"             report: {os.path.basename(run.report)}")


# ---------------------------------------------------------------- main
def default_rr_home():
    """RUBYRED_HOME, else the pinned Q16Seq copy of RubyRed, else ~/my_scripts/RubyRed."""
    if os.environ.get("RUBYRED_HOME"):
        return os.environ["RUBYRED_HOME"]
    pinned = os.path.join(HOME, "my_scripts", "RubyRed_Q16Seq")
    return pinned if os.path.isdir(pinned) else os.path.join(HOME, "my_scripts", "RubyRed")


def parse_args(argv=None):
    rr_home = default_rr_home()
    p = argparse.ArgumentParser(
        prog="RT-RubyRed",
        description="Real-time RubyRed: watch a MinKNOW fastq directory, demultiplex new "
                    "files with torchlex in batches, classify with RubyRed and keep a merged "
                    "feature table + HTML report up to date. Ctrl-C stops and finalises.")
    p.add_argument("-w", "--watch",
                   help="directory MinKNOW writes *.fastq.gz into (e.g. .../fastq_pass/barcode94)")
    p.add_argument("--from-table", metavar="TSV",
                   help="no watching: build the report (and Q16Seq tables) from an existing "
                        "RubyRed/RT-RubyRed feature table, written next to it")
    q = p.add_argument_group("Q16Seq λ internal standard (can also be set in the HTML report)")
    q.add_argument("--lambda1", type=float, metavar="COPIES",
                   help="λ1 copies per reaction (λ2 = λ1/10, λ3 = λ1/100)")
    q.add_argument("--lambda-map", metavar="TSV",
                   help="sample<TAB>λ level (1/2/3, L1, λ1)[<TAB>group]; with --lambda1 writes "
                        "<prefix>_quantitative_table.tsv and _lambda_qc.tsv")
    p.add_argument("-n", "--name", help="run name (default: MinKNOW experiment name from the path)")
    p.add_argument("-o", "--outdir", default=os.path.join(HOME, "Desktop"),
                   help="where the results folder is created (default: ~/Desktop)")
    p.add_argument("--resume", metavar="DIR", help="continue an earlier RT-RubyRed results folder")
    p.add_argument("--poll", type=int, default=60, help="seconds between directory checks (60)")
    p.add_argument("--settle", type=int, default=60,
                   help="ignore files modified in the last N seconds (60)")
    p.add_argument("--min-files", type=int, default=1,
                   help="minimum number of new files to start a batch (1)")
    p.add_argument("--idle-stop", type=float, default=0,
                   help="finalise automatically after N hours without new files (0 = never)")
    p.add_argument("-R", "--recursive", action="store_true",
                   help="also pick up files in sub-directories of the watch directory")
    p.add_argument("--cpu-frac", type=float, default=0.5, help="share of CPU cores to use (0.5)")
    p.add_argument("--mem-frac", type=float, default=0.5, help="share of RAM to allow (0.5)")
    p.add_argument("--subsample", type=int, default=1000000000,
                   help="RubyRed -t per batch (default: no subsampling)")
    p.add_argument("--classifier", default=DEFAULT_CLASSIFIER.format(rr_home=rr_home),
                   help="classifier folder (default: MIMt_16s_311025_species, strains pooled "
                        "per species); 'original' = RubyRed's own default")
    p.add_argument("--rubyred", default=os.path.join(rr_home, "RubyRed"),
                   help="RubyRed script (default: $RUBYRED_HOME/RubyRed, "
                        "RUBYRED_HOME defaults to ~/my_scripts/RubyRed)")
    p.add_argument("--rubyred-args", default="",
                   help='extra RubyRed options, e.g. "-q 12 -l 900 -x 1500"')
    p.add_argument("--barcodes", default=os.path.join(rr_home, "resources", "torchlex_barcodes.fasta"),
                   help="torchlex barcode FASTA (default: RubyRed's torchlex_barcodes.fasta)")
    p.add_argument("--torchlex", default=shutil.which("torchlex") or "torchlex")
    p.add_argument("-v", "--version", action="version", version=f"RT-RubyRed {VERSION} (Q16Seq)")
    args = p.parse_args(argv)
    args.rr_home = os.path.dirname(os.path.realpath(args.rubyred))
    if not args.watch and not args.from_table:
        p.error("one of -w/--watch or --from-table is required")
    args.lambda_map_data = report.read_lambda_map(args.lambda_map) if args.lambda_map else {}
    if args.lambda_map and not args.lambda_map_data:
        p.error(f"no sample/λ-level lines found in {args.lambda_map}")
    return args


def from_table(args):
    """Report + Q16Seq tables for an existing feature table (no RubyRed run)."""
    src = os.path.realpath(args.from_table)
    table = read_table(src)
    table = table[sorted(table.columns, key=natural_key)]
    base = os.path.splitext(src)[0]
    if base.endswith("_feature_table"):
        base = base[:-len("_feature_table")]
    name = args.name or os.path.basename(base)
    meta = {"run_name": name, "dir_name": "feature table", "watch_dir": src,
            "output_dir": os.path.dirname(src), "started": now(), "version": VERSION}
    stats = {"batches": [{"tag": "table", "files": 1, "raw_total": int(table.values.sum()),
                          "raw_reads": {s: int(table[s].sum()) for s in table.columns}}]}
    out = base + "_report.html"
    report.build(out, table, stats, meta, TOOL_DIR, lambda_cfg(args))
    lam = report.lambda_rows(table)
    say(f"report → {out}", "green")
    say(f"{table.shape[0]} taxa × {table.shape[1]} samples; λ reads: "
        f"{int(table.loc[lam].values.sum()) if lam else 0:,}", "grey")
    write_quant(args, table, base + "_feature_table.tsv")


def check_q16seq_setup(args):
    """Warn when λ reads would be lost: primer trimming in RubyRed or no λ in the classifier."""
    try:
        with open(args.rubyred, encoding="utf-8", errors="replace") as fh:
            if "--discard-untrimmed" in fh.read():
                print(f"{C['yellow']}⚠ {args.rubyred} trims reads to the 341F–1390R primers and discards "
                      f"the rest: λ internal-standard reads will be lost. Use the pinned "
                      f"~/my_scripts/RubyRed_Q16Seq (see README).{C['off']}")
    except OSError:
        pass
    clf = classifier_args(args)
    tax = clf[clf.index("-y") + 1] if "-y" in clf else None
    if tax:
        import zipfile
        try:
            with zipfile.ZipFile(tax) as z:
                name = next(n for n in z.namelist() if n.endswith("data/taxonomy.tsv"))
                has = b"lambda_standard" in z.read(name)
            print(f"{C['bold']}  λ standard{C['off']}   "
                  f"{'in classifier (S__lambda_standard)' if has else C['yellow'] + 'NOT in classifier taxonomy' + C['off']}")
        except (OSError, StopIteration, zipfile.BadZipFile):
            pass


def main(argv=None):
    global LOG_FH
    args = parse_args(argv)
    banner()
    if args.from_table:
        from_table(args)
        return
    if not os.path.isdir(args.watch):
        sys.exit(f"{C['red']}Watch directory not found: {args.watch}{C['off']}")
    for tool in (args.torchlex, "qiime"):
        if not shutil.which(tool):
            sys.exit(f"{C['red']}'{tool}' not on PATH - activate the qiime2 env "
                     f"(the RT-RubyRed launcher does this).{C['off']}")
    if args.classifier not in ("", "default", "original") and not classifier_args(args):
        if args.classifier.endswith("/MIMt_16s_311025_species"):  # default not installed
            print(f"{C['yellow']}Species classifier not found ({args.classifier}); using RubyRed's "
                  f"default. See database/README.md to install it.{C['off']}")
            args.classifier = "original"
        else:
            sys.exit(f"{C['red']}Classifier not found: {args.classifier}{C['off']}")
    res = Resources(args.cpu_frac, args.mem_frac)
    run = Run(args)
    LOG_FH = open(os.path.join(run.out, "RT-RubyRed.log"), "a")

    k = C["bold"]
    o = C["off"]
    print(f"{k}  Run{o}          {run.run_name}")
    print(f"{k}  Watching{o}     {run.watch}  "
          f"{C['grey']}(owned by {owner(run.watch)}; read-only use){o}")
    print(f"{k}  Results{o}      {run.out}")
    print(f"{k}  Merged table{o} {os.path.basename(run.merged)}")
    print(f"{k}  Report{o}       {os.path.basename(run.report)}")
    print(f"{k}  Resources{o}    {res.threads}/{os.cpu_count()} cores "
          f"(pinned {res.cpus[0]}-{res.cpus[-1]}), RAM cap {human_bytes(res.mem_limit)} of "
          f"{human_bytes(res.mem_total)} "
          f"{'(cgroup enforced)' if res.scope_ok else '(NOT enforced: systemd-run unavailable)'}"
          f", nice 10 / ionice idle-ish")
    clf = classifier_args(args)
    print(f"{k}  RubyRed{o}      {args.rubyred}")
    print(f"{k}  Classifier{o}   {os.path.basename(clf[1]) if clf else 'RubyRed default'}")
    check_q16seq_setup(args)
    if args.lambda1:
        print(f"{k}  Q16Seq{o}       λ1 = {args.lambda1:.3g} copies (λ2 {args.lambda1/10:.3g}, "
              f"λ3 {args.lambda1/100:.3g}); {len(args.lambda_map_data)} sample(s) mapped")
    print(f"{k}  Tasks{o}        every {args.poll}s: find new *.fastq.gz → torchlex demux → "
          f"fastq_demultiplexed/batchNN → RubyRed → feature_table_batchNN.tsv")
    print(f"               ≥2 batches: merge tables → {os.path.basename(run.merged)} + HTML report")
    print(f"               Ctrl-C: stop, concatenate BRK*.fastq across batches, clear batch folders")
    print()
    say(f"RT-RubyRed started (watch={run.watch}, out={run.out})", "grey")
    if run.stats["batches"]:
        say(f"resuming after batch {run.stats['batches'][-1]['batch']:02d}", "yellow")

    signal.signal(signal.SIGINT, on_sigint)
    signal.signal(signal.SIGTERM, on_sigint)
    last_new = time.time()
    last_status = 0
    while not STOP["flag"]:
        new = scan_new(run)
        if len(new) >= args.min_files:
            last_new = time.time()
            process_batch(run, res, new)
            continue  # check again straight away for files that arrived meanwhile
        if args.idle_stop and time.time() - last_new > args.idle_stop * 3600:
            say(f"no new files for {args.idle_stop} h → finishing", "yellow")
            break
        if time.time() - last_status >= args.poll:
            last_status = time.time()
            n_idx = len(run.indexed_files())
            say(f"⏳ waiting for files · {n_idx} processed, {len(new)} pending · "
                f"{len(run.stats['batches'])} batch(es)  |  {res.line()}", "grey",
                plain=f"waiting: {n_idx} processed, {len(new)} pending")
        for _ in range(min(args.poll, 10)):
            if STOP["flag"]:
                break
            time.sleep(1)
    finalise(run)


def owner(path):
    try:
        import pwd
        import grp
        st = os.stat(path)
        return f"{pwd.getpwuid(st.st_uid).pw_name}:{grp.getgrgid(st.st_gid).gr_name}"
    except (KeyError, OSError, ImportError):
        return "?"


if __name__ == "__main__":
    main()
