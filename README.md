# Scripts for genbank submissions

The following code has custom scripts that convert sequence names to a Genbank-ready 
format in three different file types: GFF (exported from Geneious), FASTA (exported from Geneious), TSV/Excel (user's annotations).

The submitters are Tillie Dunham and Mollie Burton, whose sequences and annotation formats are similar but require individual attention.

```mermaid
graph TD;
    GFF-->Py_in;
    TSV/Excel-->Py_in;
    Fasta-->Py_in;
    Py_out-->fsa;
    Py_out-->tbl;
    Py_out-->sqn;
    Py_in-->User-data-module;
    User-data-module-->Py_out;
```
