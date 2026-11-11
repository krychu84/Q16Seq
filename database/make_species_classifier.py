#!/usr/bin/env python
"""Make the MIMt 16S sklearn classifier consistent and species-level (no retraining).

1. Lineage clean-up (all ranks): curated synonyms/spelling variants are unified
   (SYNONYMS below), species-name typos of their genus are fixed, and any taxon
   filed under more than one parent lineage is moved to one lineage (most complete,
   then most common). Otherwise classify-sklearn splits the probability between
   the two spellings and stops early (e.g. Akkermansia reads stopping at phylum
   because of C__Verrucomicrobiae vs C__Verrucomicrobiia).
2. Species instead of strains (plus SPECIES_MERGE for species 16S cannot separate,
   e.g. Akkermansia massiliensis/biwaensis -> A. muciniphila, Shigella spp. -> E. coli;
   genera Escherichia + Shigella -> Escherichia-Shigella via SYNONYMS):
   Every reference keeps its own trained k-mer profile. Species names
are cleaned (strain, culture-collection, subsp. parts removed) and, where several
references share a species, the strain part moves to an 8th rank:
    ...;S__Escherichia_coli_strain_JCM_1649  ->  ...;S__Escherichia_coli;T__JCM_1649
classify-sklearn sums the probabilities of all labels under a taxonomy node, so the
species node now pools all its strains (a read is called E. coli when all E. coli
strains together reach the confidence threshold). Genus and higher ranks are
unchanged. RubyRed collapses to level 7, so strains never reach the feature table.

Usage: make_species_classifier.py <original classifier dir> <output prefix>
Output: <prefix>_classifier.qza, <prefix>_tax.qza (species only), label_map.tsv,
        lineage_changes.tsv
"""
import os
import re
import sys
import time
from collections import Counter, defaultdict
from difflib import SequenceMatcher

import numpy as np
import pandas as pd
import qiime2
from sklearn.pipeline import Pipeline


# Same taxon, two names in MIMt -> one name (current name, used by most references)
SYNONYMS = {
    "P__Euryarchaeota": "P__Methanobacteriota",
    "P__Candidatus_Thermoplasmatota": "P__Thermoplasmatota",
    "P__Nanoarchaeota": "P__Nanobdellota",
    "P__Candidatus_Nanohaloarchaeota": "P__Candidatus_Nanohalarchaeota",
    "C__Verrucomicrobiae": "C__Verrucomicrobiia",
    "C__Opitutae": "C__Opitutia",
    "C__Gemmatimonadetes": "C__Gemmatimonadia",
    "O__Candidatus_Prometheoarchaeales": "O__Promethearchaeales",
    "F__Candidatus_Prometheoarchaeaceae": "F__Promethearchaeaceae",
    "F__Nocardiopsaceae": "F__Nocardiopsidaceae",
    "F__Candidatus_Criblamydiaceae": "F__Criblamydiaceae",
    "G__Candidatus_Prometheoarchaeum": "G__Promethearchaeum",
    "G__Phytoplasma": "G__Candidatus_Phytoplasma",
    "G__Candidatus_Nitrosocosmicus": "G__Nitrosocosmicus",
    "G__Nanobsidianus": "G__Candidatus_Nanobsidianus",
    "S__Candidatus_Prometheoarchaeum_syntrophicum": "S__Promethearchaeum_syntrophicum",
    # 16S cannot separate Escherichia and Shigella (Shigella is genomically E. coli);
    # one genus, named as in SILVA
    "G__Escherichia": "G__Escherichia-Shigella",
    "G__Shigella": "G__Escherichia-Shigella",
}

# Species that 16S cannot separate, reported as one species (the others become
# T__ "strains" of the target, so their probabilities are pooled). Add pairs here.
SPECIES_MERGE = {
    "S__Akkermansia_massiliensis": "S__Akkermansia_muciniphila",
    "S__Akkermansia_biwaensis": "S__Akkermansia_muciniphila",
    "S__Shigella_flexneri": "S__Escherichia_coli",
    "S__Shigella_sonnei": "S__Escherichia_coli",
    "S__Shigella_boydii": "S__Escherichia_coli",
    "S__Shigella_dysenteriae": "S__Escherichia_coli",
}


def merge_species(label):
    ranks = label.split(";")
    if ranks[-1] in SPECIES_MERGE:
        ranks[-1] = SPECIES_MERGE[ranks[-1]]
    return ";".join(ranks)


def harmonise(labels):
    """Return cleaned lineages (same order) and a Counter of what was changed."""
    stats = Counter()
    rows = []
    for lab in labels:
        r = [x.strip() for x in lab.split(";")]
        for i, x in enumerate(r):
            if x in SYNONYMS:
                r[i] = SYNONYMS[x]
                stats[f"synonym {x} -> {SYNONYMS[x]}"] += 1
        rows.append(r)
    # species named after a misspelled or older genus: use the genus of G__
    by_epithet = defaultdict(set)
    for r in rows:
        if len(r) >= 7 and len(r[5]) > 3 and len(r[6]) > 3:
            w = r[6][3:].split("_", 1)
            if len(w) == 2 and w[0] == r[5][3:]:
                by_epithet[(r[5], w[1].split("_strain_")[0])].add(w[0])
    for r in rows:
        if len(r) >= 7 and len(r[5]) > 3 and len(r[6]) > 3:
            g = r[5][3:]
            w = r[6][3:].split("_", 1)
            if len(w) == 2 and w[0] != g and not w[0].startswith(("Candidatus", "[")):
                epi = w[1].split("_strain_")[0]
                if SequenceMatcher(None, w[0], g).ratio() >= 0.9 or (r[5], epi) in by_epithet:
                    stats[f"species genus word {w[0]} -> {g}"] += 1
                    r[6] = f"S__{g}_{w[1]}"
    # every named taxon gets exactly one parent lineage (top-down, repeat until stable)
    for _ in range(5):
        changed = 0
        for i in range(1, 7):
            parents = defaultdict(Counter)
            for r in rows:
                if len(r) > i and len(r[i]) > 3:
                    parents[r[i]][tuple(r[:i])] += 1
            best = {}
            for name, c in parents.items():
                if len(c) > 1:
                    best[name] = min(c, key=lambda p: (sum(len(x) <= 3 for x in p), -c[p], p))
            for r in rows:
                if len(r) > i and r[i] in best and tuple(r[:i]) != best[r[i]]:
                    stats[f"parent of {r[i]}"] += 1
                    r[:i] = list(best[r[i]])
                    changed += 1
        if not changed:
            break
    return [";".join(r) for r in rows], stats


def species(label):
    """S__Genus_species_strain_X / _DSM_123 / _subsp._y -> S__Genus_species."""
    ranks = label.split(";")
    if not ranks[-1].startswith("S__") or ranks[-1] == "S__":
        return label
    s = ranks[-1][3:]
    s = re.split(r"_(?:strain|str\.)_", s)[0]
    t = s.split("_")
    if t[0] == "Candidatus" and len(t) >= 3:
        s = "_".join(t[:3])
    elif len(t) >= 2 and t[1] in ("sp.", "spp.", "sp"):
        s = t[0] + "_sp."
    elif len(t) >= 3 and t[1] in ("genomosp.", "genosp."):
        s = "_".join(t[:3])
    elif len(t) >= 2 and re.fullmatch(r"[\['\"]?[A-Z][A-Za-z\-]*[\]'\"]?", t[0]) \
            and re.fullmatch(r"[a-z][a-z\-]+'?", t[1]):
        s = "_".join(t[:2])
    ranks[-1] = "S__" + s
    return ";".join(ranks)


def with_strain_rank(old, new):
    """Species label, plus a T__ strain rank when a species has several references."""
    n_per_sp = Counter(new)
    seen = defaultdict(int)
    out = []
    for o, sp in zip(old, new):
        if n_per_sp[sp] == 1:
            out.append(sp)
            continue
        strain = o.split(";")[-1][3:]
        base = sp.split(";")[-1][3:]
        strain = strain[len(base):] if strain.startswith(base) else strain
        strain = re.sub(r"^_(?:strain_|str\._)?", "", strain) or "ref"
        lab = f"{sp};T__{strain}"
        seen[lab] += 1
        out.append(lab if seen[lab] == 1 else f"{lab}_{seen[lab]}")
    return np.array(out)


def main():
    src, pfx = sys.argv[1].rstrip("/"), sys.argv[2]
    name = os.path.basename(src)
    t0 = time.time()
    pipe = qiime2.Artifact.load(f"{src}/{name}_classifier.qza").view(Pipeline)
    nb = pipe.steps[-1][1]
    old = np.asarray(nb.classes_)
    clean, stats = harmonise(old)
    new = np.array([merge_species(species(c)) for c in clean])
    final = with_strain_rank(clean, new)
    assert len(set(final)) == len(final)
    order = np.argsort(final, kind="stable")  # q2 requires sorted labels
    nb.classes_ = final[order]
    for attr in ("feature_count_", "feature_log_prob_", "class_count_", "class_log_prior_"):
        setattr(nb, attr, getattr(nb, attr)[order])
    qiime2.Artifact.import_data("TaxonomicClassifier", pipe).save(f"{pfx}_classifier.qza")
    print(f"{len(old)} references -> {len(set(new))} species "
          f"({sum(';T__' in f for f in final)} with a T__ strain rank)  ({time.time()-t0:.0f}s)")

    tax = qiime2.Artifact.load(f"{src}/{name}_tax.qza").view(pd.DataFrame)
    tax_clean, _ = harmonise(tax["Taxon"])
    tax["Taxon"] = [merge_species(species(t)) for t in tax_clean]
    qiime2.Artifact.import_data("FeatureData[Taxonomy]", tax).save(f"{pfx}_tax.qza")
    pd.DataFrame({"original": old, "classifier_label": final}) \
        .query("original != classifier_label") \
        .to_csv(os.path.join(os.path.dirname(pfx), "label_map.tsv"), sep="\t", index=False)
    lin = lambda t: ";".join(t.split(";")[:6])
    pd.DataFrame({"original": old, "cleaned": clean})[[lin(a) != lin(b) or a.split(";")[-1] != b.split(";")[-1]
                                                       for a, b in zip(old, clean)]] \
        .to_csv(os.path.join(os.path.dirname(pfx), "lineage_changes.tsv"), sep="\t", index=False)
    print(f"lineage clean-up: {sum(a != b for a, b in zip(old, clean))} references changed")
    for k, v in sorted(stats.items(), key=lambda kv: -kv[1]):
        print(f"   {v:5d}  {k}")
    print("done")


if __name__ == "__main__":
    main()
