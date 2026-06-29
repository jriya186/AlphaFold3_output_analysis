# af3toolkit

Post-AlphaFold3 structural contact analysis for LDL-receptor-family / Reelin
complexes (ApoER2, LDLR, VLDLR).

AF3 inputs for this kind of analysis are usually **spliced fragments**, not
full-length proteins (e.g. only exons 2-8 of a receptor, or only repeats 5-6
of Reelin). That means raw PDB residue numbers do **not** correspond to real
positions in the full-length protein. This toolkit re-maps every contact
back to:

- the correct full-sequence position
- the exon (receptor) or repeat (Reelin) it falls in
- the structural domain it belongs to (receptor side)
- whether that region has experimental structural support or is AF3
  prediction only (Reelin side)

and explicitly flags which domains weren't even part of the modeled
construct, so "no contacts in domain X" is never confused with "no binding
in domain X."

## Install

```bash
pip install -e .
# optional, for nicer CLI tables:
pip install -e ".[cli-pretty]"
```

## CLI usage

```bash
af3toolkit analyze \
  --pdb apoer2_reelin56_complex.pdb \
  --receptor ApoER2 \
  --receptor-exons 2,3,4,5,6,7,8 \
  --reelin-repeats 5,6 \
  --save-csv results/apoer2_repeat56_contacts.csv
```

`--save-csv` is optional — omit it and the full table just prints to
terminal without writing a file.

Discover what's already configured:

```bash
af3toolkit list-receptors
af3toolkit list-reelin-repeats
```

## Output

Every contact row includes **both** the receptor and Reelin side, fully
re-mapped:

| Receptor Residue | Receptor PDB Pos | Receptor Full-Seq Pos | Receptor Exon | Receptor Domain | Reelin Residue | Reelin PDB Pos | Reelin Full-Seq Pos | Reelin Repeat | Reelin Validated |
|---|---|---|---|---|---|---|---|---|---|

Above the table, the CLI also prints:
- a **domain-level contact summary** (counts per receptor domain)
- a list of **domains excluded** from the construct (not modeled, not evidence of no binding)
- **calcium-coordinating residues**, since LA repeats are calcium-dependent

## Why receptor and Reelin sides are treated differently

Only Reelin repeats 5-6 have a published crystal structure of the actual
receptor-bound interface. Every other repeat's contacts in an AF3 model are
**predictions with no experimental structure to check them against**. The
`Reelin Validated` column makes that distinction explicit in the output,
rather than implying uniform confidence across the whole protein.

## Before using this on your own data

- `af3toolkit/receptors.py`: exon→domain labels are approximate (based on
  general LDL-receptor-family architecture) and marked `TODO: verify` —
  check against each receptor's UniProt feature table before citing exact
  boundaries.
- `af3toolkit/reelin.py`: `full_seq_range` for each repeat is `None` by
  default — fill these in from your own records before running an
  analysis; the pipeline will raise a clear error if you try to use a
  repeat that isn't filled in yet, rather than silently mis-mapping it.

## Repo structure

```
af3toolkit/          # the package
  receptors.py        # exon coordinates + domain labels (ApoER2, LDLR, VLDLR)
  reelin.py            # repeat coordinates + validated/prediction-only tags
  sequence_mapping.py  # generic local-PDB-position -> full-sequence mapper
  contacts.py          # Biopython-based contact + ion-coordination detection
  confidence.py        # best-AF3-model selection by iPTM/pTM
  pipeline.py           # ties it all together: analyze_complex()
  cli.py                # command-line interface
notebooks/archive/    # original exploratory notebooks this package was
                       # built from (kept for reference, not maintained)
tests/                 # pytest tests for the position-mapping logic
```

## Development

```bash
pip install -e ".[dev]"
pytest tests/
```
