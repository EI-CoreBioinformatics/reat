# ORF Calling

As part of generating accurate gene models from the initial transcriptomic data, REAT incorporates the process of ORF calling as one of it's steps. This can be achieved specifically by using either `Transdecoder` (TODO add citation) or `Prodigal` (TODO add citation).

## Transdecoder

The original Transdecoder software has been reorganised in this workflow by using it's parts to build a workflow that allows for better paralelisation of the different steps. More specifically: Longest ORF identification, ORF Scoring, BLAST alignments, ORF Selection and start codon identification can now be distributed to multiple workers. This can be done using a user provided chunking parameter that allows users to easily adjust the resource usage of this workflow to their available infrastructure.

## Prodigal

This workflow uses the prodigal software without modifications to annotate potential ORFs on the transcripts.