# Subworkflows

All the trancriptome module tasks and subworkflows are included in this directory. Using subworkflows helps maximise the reusability of the different tasks under different conditions minimising the need to repeat code and providing more opportunity for the tasks calls to be cached and hence reused.

## Structure

* align_long - Long read alignment workflow
* assembly_long - Long read assembly workflow
* align_short - Short read alignment workflow
* assembly_short - Short read assembly workflow
* common - Contains tasks and structs (helpful property groupings) which are reused across different workflows.
* exonerate - *TODO: Move to a different workflow*
* mikado - Defines the Mikado workflow, also includes ORF calling and Homology mapping workflows which are used as part of Mikado inputs.
* portcullis - Portcullis workflow, includes a regrouping of samples. More details found in the portcullis subdirectory.
* repeat_masker - *TODO: Move to a different workflow*
* sanitise - Several tasks in this workflow have specific requirements for the information on the reference genome and if present a base annotation, this workflow ensures files are valid inputs to all tasks.

More details for each subworkflow can be found within each subdirectory.