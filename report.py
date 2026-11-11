"""Self-contained HTML report for Q16Seq / RT-RubyRed (plotly.js inlined, works offline)."""
import datetime as dt
import html
import json
import os
import re

RANKS = {"K": "kingdom", "P": "phylum", "C": "class", "O": "order", "F": "family",
         "G": "genus", "S": "species"}


def taxon_label(lineage):
    """Species name if resolved, otherwise the lowest named rank."""
    named = []
    for part in lineage.split(";"):
        m = re.match(r"\s*([A-Za-z])__(.*)", part)
        if m and m.group(2).strip():
            named.append((m.group(1).upper(), m.group(2).strip().replace("_", " ")))
    if not named:
        return "Unassigned"
    rank, name = named[-1]
    if rank == "S":
        return name
    if rank == "G":
        return f"{name} sp."
    return f"{name} ({RANKS.get(rank, rank)})"


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


# ---------------------------------------------------------------- lambda internal standard
LAMBDA_RE = re.compile(r"lambda", re.I)
LEVELS = {"1": 1, "L1": 1, "L01": 1, "λ1": 1, "LAMBDA1": 1,
          "2": 2, "L2": 2, "L02": 2, "λ2": 2, "LAMBDA2": 2,
          "3": 3, "L3": 3, "L03": 3, "λ3": 3, "LAMBDA3": 3}
ROQ = (0.10, 0.90)  # lambda reads must be 10-90% of all classified reads (Q16Seq ROQ)


def lambda_rows(table):
    """Lineages of the lambda internal standard (S__lambda_standard) in a feature table."""
    return [lin for lin in table.index if LAMBDA_RE.search(str(lin))]


def lambda_copies(level, l1):
    """Copies per reaction of standard lambda1/2/3: lambda2 = lambda1/10, lambda3 = lambda1/100."""
    return float(l1) / (10 ** (level - 1))


def read_lambda_map(path):
    """TSV/CSV: sample, lambda level (1/2/3, L1, λ1, ...), optional group. '#' lines skipped."""
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = [c.strip() for c in re.split(r"[\t,;]", line.rstrip("\n"))]
            if len(cols) < 2:
                continue
            lvl = LEVELS.get(cols[1].upper().replace(" ", ""))
            if lvl is None:  # header line or unknown level
                continue
            out[cols[0]] = {"level": lvl, "group": cols[2] if len(cols) > 2 else ""}
    return out


def quantify(table, l1, mapping):
    """Absolute 16S copies per reaction for every taxon and sample with a lambda level.

    copies(taxon, sample) = lambda copies added x reads(taxon) / reads(lambda)
    Returns (copies DataFrame without the lambda rows, per-sample QC DataFrame).
    """
    import pandas as pd
    lam = lambda_rows(table)
    lam_reads = table.loc[lam].sum(axis=0) if lam else pd.Series(0, index=table.columns)
    bact = table.drop(index=lam)
    qc, cols = [], {}
    for s in table.columns:
        m = mapping.get(s)
        total = int(table[s].sum())
        lr = int(lam_reads[s])
        row = {"sample": s, "group": (m or {}).get("group", ""),
               "lambda_level": f"L{m['level']}" if m else "", "lambda_copies": None,
               "lambda_reads": lr, "bacterial_reads": total - lr,
               "lambda_fraction": round(lr / total, 4) if total else None, "status": ""}
        if not m:
            row["status"] = "no lambda level"
        elif lr == 0:
            row["status"] = "no lambda reads"
        else:
            row["lambda_copies"] = lambda_copies(m["level"], l1)
            cols[s] = bact[s] * row["lambda_copies"] / lr
            f = lr / total
            row["status"] = "ok" if ROQ[0] <= f <= ROQ[1] else "outside ROQ (lambda 10-90%)"
        qc.append(row)
    copies = pd.DataFrame(cols, index=bact.index)
    copies.loc["Total bacterial 16S copies"] = copies.sum(axis=0)
    return copies, pd.DataFrame(qc).set_index("sample")


def build(out_path, table, stats, meta, tool_dir, lam=None):
    labels, seen = [], {}
    for lin in table.index:
        lab = taxon_label(lin)
        if lab in seen:
            seen[lab] += 1
            lab = f"{lab} [{seen[lab]}]"
        else:
            seen[lab] = 1
        labels.append(lab)

    batches = [b for b in stats.get("batches", [])]
    raw = {}
    for b in batches:
        for s, n in b.get("raw_reads", {}).items():
            raw.setdefault(s, {})[b["tag"]] = n
    samples = sorted(set(table.columns) | set(raw), key=natural_key)
    t = table.reindex(columns=samples, fill_value=0)

    data = {
        "samples": samples,
        "taxa": labels,
        "lineages": list(table.index),
        "counts": t.values.astype(int).tolist(),
        "raw": raw,
        "batches": [b["tag"] for b in batches],
        "lambda": {
            "rows": [i for i, lin in enumerate(table.index) if LAMBDA_RE.search(str(lin))],
            "l1": (lam or {}).get("l1"),
            "map": (lam or {}).get("map") or {},
            "roq": list(ROQ),
        },
        "meta": {
            "run": meta.get("run_name", ""),
            "dir": meta.get("dir_name", ""),
            "watch": meta.get("watch_dir", ""),
            "out": meta.get("output_dir", ""),
            "started": meta.get("started", ""),
            "updated": dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "files": sum(b.get("files", 0) for b in batches),
            "raw_total": sum(b.get("raw_total", 0) for b in batches),
            "version": meta.get("version", ""),
        },
    }

    assets = os.path.join(tool_dir, "assets")

    def asset(name):
        with open(os.path.join(assets, name), encoding="utf-8") as fh:
            return fh.read()

    page = (asset("report_template.html")
            .replace("<!--__RUBYRED_LOGO__-->", strip_xml(asset("rubyred_logo.svg")))
            .replace("<!--__KU_LOGO__-->", strip_xml(asset("ku_seal.svg")))
            .replace("__TITLE__", html.escape(f"Q16Seq · {data['meta']['run']} · {data['meta']['dir']}"))
            .replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
            .replace("/*__PLOTLY__*/", asset("plotly.min.js").replace("</script", "<\\/script")))
    tmp = out_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(page)
    os.replace(tmp, out_path)


def strip_xml(svg):
    return re.sub(r"<\?xml[^>]*\?>|<!DOCTYPE[^>]*>|<!--.*?-->", "", svg, flags=re.S).strip()
