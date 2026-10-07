# Provenance and attribution — `inst/extdata`

These files ship inside the installed package, which makes them redistribution. Two of
the upstream sources ask for attribution in return; this file provides it.

| File | Derived from | Terms |
|---|---|---|
| `uniprot_annotated_protein_info.txt` | UniProtKB, plus GO / KEGG / eggNOG / CAZy / Pfam / InterPro cross-references and NCBI Taxonomy lineages | CC BY 4.0 — attribution required |
| `taxonomy.txt` | NCBI Taxonomy | US Government public domain — no restrictions |
| `report.pg_matrix.tsv` | DIA-NN output, lab-generated | Ours; see note below |
| `report.pr_matrix.txt` | DIA-NN output, lab-generated | Ours; see note below |
| `diann.parquet` | DIA-NN output, lab-generated | Ours; see note below |

## Attribution

**UniProtKB.** `uniprot_annotated_protein_info.txt` holds 29,157 protein records —
accessions, sequences, protein names, organism names and functional cross-references —
retrieved from UniProtKB and reshaped into a flat table. UniProt is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), which permits redistribution
of derived work with attribution.

> The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase.*
> <https://www.uniprot.org/> — © UniProt Consortium, CC BY 4.0.

The GO terms in the `go` column come from the [Gene Ontology](https://geneontology.org/)
(CC BY 4.0, © the Gene Ontology Consortium); InterPro and Pfam accessions are
[CC0](https://creativecommons.org/publicdomain/zero/1.0/). The `xref_kegg` and
`xref_cazy` columns are **identifiers only** — cross-reference keys carried through from
UniProt, not KEGG or CAZy content — so neither resource's academic-use restriction
attaches to this file.

**NCBI Taxonomy.** `taxonomy.txt` maps 545 organism identifiers to their ranked lineage.
NCBI Taxonomy is a work of the US Government and is in the public domain; no permission
is needed and none is claimed here. Citation is nonetheless appreciated:

> Schoch CL *et al.* *NCBI Taxonomy: a comprehensive update on curation, resources and
> tools.* Database (2020). <https://www.ncbi.nlm.nih.gov/taxonomy>

## The DIA-NN outputs

`report.pg_matrix.tsv`, `report.pr_matrix.txt` and `diann.parquet` are search results
from an in-house metaproteomics experiment, produced by running DIA-NN over lab-acquired
spectra. The measurements are the lab's own.

Two things about them are worth stating plainly:

- They are **output of** DIA-NN, not any part of DIA-NN itself. No DIA-NN code, binary
  or model is present here, and none of DIA-NN's proprietary terms reach these files.
- Their `Protein.Group`, `Protein.Names`, `Genes` and `First.Protein.Description` columns
  are identifiers and descriptions carried over from the UniProt FASTA the search ran
  against, so the UniProt attribution above covers them too.

Raw file paths visible in the column headers are the acquisition paths from the
instrument; they are provenance, not a reference to anything shipped.

## Scope

This file covers `inst/extdata` only. Test fixtures under `tests/testthat/fixtures/`
draw on the same three sources and no others — UniProt protein records, NCBI taxonomy
names and lineages, and lab-generated search output — so the attribution above applies to
them equally. Fixtures are not installed with the package; they ship only in the source
tarball and the repository.

If a fixture is ever added from a source not listed above, add it to the table.
