# MIMt 16S species-level classifier

`make_species_classifier.py` turns RubyRed's MIMt 16S sklearn classifier (`MIMt_16s_311025`) into a cleaned, species-level classifier (`MIMt_16s_311025_species`) **without retraining**. It takes about 1 minute and needs about 9 GB of RAM.

```bash
conda activate qiime2
C=~/my_scripts/RubyRed/resources/classifiers
mkdir -p $C/MIMt_16s_311025_species
python make_species_classifier.py $C/MIMt_16s_311025 $C/MIMt_16s_311025_species/MIMt_16s_311025_species
ln -s ../MIMt_16s_311025/MIMt_16s_311025_seqs.qza $C/MIMt_16s_311025_species/MIMt_16s_311025_species_seqs.qza
```

Outputs, in the target folder:

| File | Content |
|---|---|
| `MIMt_16s_311025_species_classifier.qza` | the cleaned classifier (for `classify-sklearn`, RubyRed `-c`) |
| `MIMt_16s_311025_species_tax.qza` | cleaned species-level taxonomy (RubyRed `-y`, used by vsearch/BLAST) |
| `label_map.tsv` | every classifier label that changed: original → new |
| `lineage_changes.tsv` | references whose lineage was corrected |

Ready-made `.qza` files are attached to the [GitHub releases](https://github.com/krychu84/RT-RubyRed/releases). They were built from `MIMt_16s_311025` with QIIME 2 2024.5 and scikit-learn 1.4.2. Use them only with the same scikit-learn version; with any other version, build the classifier yourself.

## Why

QIIME's `classify-sklearn` works down the ranks (kingdom → … → species) and reports the deepest rank where the summed probability of one name reaches the confidence threshold (RubyRed: 0.9). A read stops early whenever the probability is **split between two labels that mean the same thing**:

- **Strain-level species names.** 23,593 of 31,436 MIMt species labels carry a strain name (`S__Escherichia_coli_strain_JCM_1649`) or a culture-collection number. Each strain is a separate species label, so a read that fits a species with several strains is split between them and falls back to genus.
- **The same taxon under two names or lineages.** For example, *Akkermansia* was filed under both `C__Verrucomicrobiae` and `C__Verrucomicrobiia`. The classifier was 100% sure the reads were *Akkermansia*, but split 16/84 between the two class spellings, so every such read stopped at **phylum**.
- **Species 16S cannot resolve.** *A. muciniphila* reads shared ~16% of their probability with *A. massiliensis*, and *E. coli* reads were split between the genera *Escherichia* and *Shigella*.

## How (no retraining)

Naive Bayes stores one trained k-mer profile per label. The script keeps every profile unchanged and edits only the labels:

1. **Lineage clean-up (`harmonise`)**
   - **`SYNONYMS`**: a curated list of spelling variants and renamed taxa, mapped to the current name used by most references:

     | Rank | Original → new |
     |---|---|
     | phylum | `Euryarchaeota` → `Methanobacteriota` |
     | phylum | `Candidatus_Thermoplasmatota` → `Thermoplasmatota` |
     | phylum | `Nanoarchaeota` → `Nanobdellota` |
     | phylum | `Candidatus_Nanohaloarchaeota` → `Candidatus_Nanohalarchaeota` |
     | class | `Verrucomicrobiae` → `Verrucomicrobiia` |
     | class | `Opitutae` → `Opitutia` |
     | class | `Gemmatimonadetes` → `Gemmatimonadia` |
     | family | `Nocardiopsaceae` → `Nocardiopsidaceae` |
     | family | `Candidatus_Criblamydiaceae` → `Criblamydiaceae` |
     | order → species | `Candidatus_Prometheoarchae*` → `Promethearchae*` |
     | genus | `Phytoplasma` → `Candidatus_Phytoplasma` |
     | genus | `Candidatus_Nitrosocosmicus` → `Nitrosocosmicus` |
     | genus | `Nanobsidianus` → `Candidatus_Nanobsidianus` |
     | genus | `Escherichia` + `Shigella` → `Escherichia-Shigella` |
   - **Species named after a misspelled or older genus** are renamed to the genus they are filed under, but only when that genus already holds the same species under its current name (e.g. `Flaviaesturariibacter` → `Flaviaestuariibacter`, `Kinneretia asaccharophila` → `Roseateles asaccharophila`).
   - **Every taxon gets one parent lineage:** the most complete one (fewest empty ranks), then the most common (e.g. *Fictibacillus*, *Hoeflea*, *Hungatella*, ~25 archaeal genera).
   - Result: every name at every rank has exactly one parent. 198 of 27,677 labels were corrected.

2. **Species instead of strains (`species`)**
   - Strain names, culture-collection numbers (DSM, ATCC, JCM, …), `subsp.` and isolate codes are removed, giving *Genus species*.
   - `Genus_sp._X` → `Genus sp.`. `Candidatus` names keep three words, genomospecies keep their designation, and bracketed genera (`[Clostridium]_nexile`) are handled.

3. **Species 16S cannot separate (`SPECIES_MERGE`)**

   | Original | Reported as |
   |---|---|
   | *Akkermansia massiliensis*, *A. biwaensis* | *Akkermansia muciniphila* |
   | *Shigella flexneri*, *sonnei*, *boydii*, *dysenteriae* | *Escherichia coli* |

4. **Strains become an 8th rank (`with_strain_rank`)**
   - Where several references share a species, the strain part becomes an extra rank: `…;S__Escherichia_coli;T__JCM_1649`.
   - `classify-sklearn` adds up the probabilities of all labels below a taxon, so the species now pools all its strains. A read is called *E. coli* when all *E. coli* references together reach the threshold.
   - Genus and higher ranks keep exactly the probabilities they had before.
   - RubyRed collapses tables at level 7, so strains never reach the feature table. RT-RubyRed also pools the per-strain `rep_seqs` files by species.
   - Labels stay unique and sorted, as `classify-sklearn` requires. **Do not** simply give strains identical names: duplicate labels break its confidence calculation.
   - Adding up the per-strain counts into one profile per species (equivalent to retraining on species labels) was tried and rejected: the averaged profile made thousands of *Akkermansia* reads drop to phylum.

## Test

The same 49,775 reads were classified with both classifiers (confidence 0.9):

| Classified to | original | cleaned |
|---|---|---|
| species | 8.3 % | 41.5 % |
| genus | 58.8 % | 44.7 % |
| above genus | 32.9 % | 13.8 % |
| stopping at phylum | 9,544 | 1 |
| *Enterobacteriaceae* (family only) | 6,759 | 876 |
| strain names in the output | 4,129 reads | 0 |

- No read became less specific, and every read called to species by both classifiers got the same species.
- Most remaining reads at class or order level were PCR chimeras: each half of the read clearly belonged to a different genus (e.g. *Staphylococcus* + *Limosilactobacillus*).

## Adding your own merges

To report two species (or genera) as one, add them to `SPECIES_MERGE` (species) or `SYNONYMS` (any rank) at the top of the script, then re-run it:

```python
SPECIES_MERGE = {
    "S__Akkermansia_massiliensis": "S__Akkermansia_muciniphila",
    "S__Some_species_b": "S__Some_species_a",   # new
}
```

Check a merge in the data first. For example, look at `predict_proba` for the reads that stop at genus to see which species share their probability.
