# Incorporating homology data

Previously characterised genomic information from closely related species can help the process of annotation by identifying transcripts which are very similar to previously identified models.

REAT can make use of a collection of proteins from closely related species to help the process of annotation. This is done by using a curated database to align the prepared transcripts resulting from the `prepare` step of Mikado. The transcripts are then linked to different proteins in the database based on their similarity. This information is then finally used in Mikado's `pick` stage to select the most informative gene models.

Given the potential large number of trancripts that need to be aligned to the database, this step allows the user to provide a chunking parameter, which will be used to specify the number of transcripts to be aligned by a single worker in this scatter step. Selecting the chunking parameter will allow users to tune this step to the available hardware ensuring best usage of their infrastructure.