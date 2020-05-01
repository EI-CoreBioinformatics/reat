# Transcriptome module workflows

The transcriptome module is divided in three general steps: Read alignment, transcript assembly and Mikado runs. These steps are represented in the files `align.wdl` and `mikado.wdl`.

* `align.wdl` - Defines the alignment of reads and assembling of transcripts.
* `mikado.wdl` - Defines the evaluation and scoring of transcripts, along with the inclusion of external evidence such as species related proteins.

## Align workflow (`align.wdl`)

Define in words what this is more specifically

Mention all the possible choices for aligners (short, and long) and assemblers (short, and long) along with the filtering options (merge, merge_collapse and filter_gff) 

A diagram would be helpful here.

## Mikado workflow (`mikado.wdl`)

Define in words what this is. Link to the Mikado paper as a reference.

Add a diagram of mikado + ancilliary processes (ORFs + Homology).
