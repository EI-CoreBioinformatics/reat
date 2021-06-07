
# REAT - Robust and Extendable eukaryotic Annotation Toolkit

REAT is a robust easy-to-use genome annotation toolkit for turning high-quality genome assemblies into usable and informative resources. REAT makes use of state-of-the-art annotation tools and is robust to varying quality and sources of molecular evidence.

REAT provides an integrated environment that comprises both 
 * a set of workflows geared towards integrating multiple sources of evidence into a genome annotation, and 
 * an execution environment (Docker, Singularity) for these workflows.

## Getting Started

To configure a development environment for REAT you need:

* the REAT source code
* Software dependencies
* A WDL compatible workflow engine (REAT development is mainly geared towards using Cromwell)

To obtain the REAT source code, please:

```
git clone https://github.com/ei-corebioinformatics/reat
```

While an exhaustive list of software dependencies is available in the appendix, a simpler way of obtaining all the required binaries for REAT is to build a Singularity container based on the definition file (`reat_singularity.def`) included in the repository.

```
singularity build reat.img reat_singularity.def
```

Cromwell can be obtained from the original repository at https://github.com/broadinstitute/cromwell/releases.

### Prerequisites

For ease of development, Singularity is recommended. In case this is not available (you are working on MacOS and haven't got https://sylabs.io/guides/3.5/admin-guide/installation.html#mac available), please follow the instructions in the `REAT-REAT-Singularity.def` for building and installing the binary dependencies.

#### Java 8

#### Singularity (optional)


#### Cromwell >= v48

You will need to download this from the Cromwell [repository](https://github.com/broadinstitute/cromwell/releases)

### Installing

First obtain the source code using

```
git clone https://github.com/ei-corebioinformatics/reat
cd reat
```

You should now have the latest development version, in which you will find the `reat_singularity.def`. This allows you to build a Singularity container. 

We also provide a conda environment which you can use to install dependencies if you prefer.

```shell
conda env create -f reat.yaml --name reat
```

Alternatively, install the following software dependencies on your system, so that the executables are available  in your `$PATH` environment variable.

* DIAMOND - 0.9.31
* BioPerl
* FullLengtherNext - 1.0.1
* Install BLAST - 2.7.1
* Install MagicBlast - 1.5.0
* LibDeflate - master
* HTSLib - 1.9
* Samtools - 1.9
* BCFTools - 1.9
* BAMtools - 2.5.1
* gclib - 54917d0
* gffread - ba7535f 
* GMAP - 2019-02-15
* MiniMap2 - 2.17
* GenomeTools - 1.5.10
* HISAT2 - 2.1.0
* STAR - 2.7.3a
* seqtk - master
* Stringtie2 - v2.1.1
* Scallop - 0.10.4
* Scallop-lr - 0.9.2
* Prodigal - 2.6.3
* Transdecoder - 5.5.0
* Portcullis - 1.2.2
* Mikado - 2.0

Instructions for building each of the tools can be found in the `reat_singularity.def`.

It is recommended to make use of `virtualenv` when installing python packages.
An alternative to manual compilation of all tools listed is using conda.

Please install `reat` to ensure the scripts are findable during execution. To install, simply use from your current pip/conda environment:

```shell
pip install ./reat
```

## Running REAT

Cromwell supports a range of execution backends such as cloud services and batch schedulers. 
These can be configured in files such as those present in the `cromwell_configuration` directory.
For more examples you can also visit https://github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends

### Local

Reat finds and executes a binary named `cromwell` which can be configured by placing a file similar to:

```sh
#!/bin/sh
java -Dconfig.file=reat/cromwell_configuration/local.conf -jar cromwell.jar
```

Where `local.conf` is the file in this repository under `cromwell_configuration` and `cromwell.jar` is a cromwell release jar file (https://github.com/broadinstitute/cromwell/releases)


Finally, for the homology module you can run:
```
reat --workflow_options_file workflow_options/options.json --computational_resources compute_inputs.json homology <homology_args>
```

And, for the transcriptome module:
```
reat --workflow_options_file workflow_options/options.json --computational_resources compute_inputs.json transcriptome <transcriptome_args>
```

### SLURM

To run REAT in a SLURM backend the cromwell instance needs to be configured using the configuration file equivalent 
for your environment to `cromwell_configuration/slurm.conf`. This file contains the configuration required by Cromwell 
to execute jobs using a SLURM backend. It includes details such as how to submit, check and kill a job running on a 
SLURM backend.

Running reat under slurm requires changing the `-Dconfig.file` argument from the local example of our `cromwell` command to:

```sh
#!/bin/sh
java -Dconfig.file=reat/cromwell_configuration/slurm.conf -jar cromwell.jar
```

### Inputs

#### Transcriptome module
```
usage: reat transcriptome [-h] [--samples SAMPLES [SAMPLES ...]] [--csv_paired_samples CSV_PAIRED_SAMPLES] [--csv_long_samples CSV_LONG_SAMPLES] --reference REFERENCE [--annotation ANNOTATION]
                          [--parameters_file PARAMETERS_FILE] [--all_extra_config ALL_EXTRA_CONFIG] [--long_extra_config LONG_EXTRA_CONFIG] [--lq_extra_config LQ_EXTRA_CONFIG] --all_scoring_file
                          ALL_SCORING_FILE [--long_scoring_file LONG_SCORING_FILE] [--long_lq_scoring_file LONG_LQ_SCORING_FILE] [--homology_proteins HOMOLOGY_PROTEINS]
                          [--separate_mikado_LQ SEPARATE_MIKADO_LQ] [--short_reads_aligner {hisat,star}] [--HQ_aligner {minimap2,gmap}] [--LQ_aligner {minimap2,gmap}] [--min_identity MIN_IDENTITY]
                          [--min_intron_len MIN_INTRON_LEN] [--max_intron_len MAX_INTRON_LEN] [--max_intron_len_ends MAX_INTRON_LEN_ENDS] [--PR_hisat_extra_parameters PR_HISAT_EXTRA_PARAMETERS]
                          [--PR_star_extra_parameters PR_STAR_EXTRA_PARAMETERS] [--HQ_aligner_extra_parameters HQ_ALIGNER_EXTRA_PARAMETERS] [--LQ_aligner_extra_parameters LQ_ALIGNER_EXTRA_PARAMETERS]
                          [--skip_scallop] [--HQ_assembler {filter,merge,stringtie,stringtie_collapse}] [--LQ_assembler {filter,merge,stringtie,stringtie_collapse}]
                          [--HQ_assembler_extra_parameters HQ_ASSEMBLER_EXTRA_PARAMETERS] [--LQ_assembler_extra_parameters LQ_ASSEMBLER_EXTRA_PARAMETERS]
                          [--PR_stringtie_extra_parameters PR_STRINGTIE_EXTRA_PARAMETERS] [--PR_scallop_extra_parameters PR_SCALLOP_EXTRA_PARAMETERS] [--extra_parameters EXTRA_PARAMETERS]
                          [--orf_caller {prodigal,transdecoder,none}] [--orf_calling_proteins ORF_CALLING_PROTEINS]

optional arguments:
  -h, --help            show this help message and exit
  --samples SAMPLES [SAMPLES ...]
                        Reads organised in the input specification for REAT, for more information please look at https://github.com/ei-corebioinformatics/reat for an example (default: None)
  --csv_paired_samples CSV_PAIRED_SAMPLES
                        CSV formatted input paired read samples. Without headers. The CSV fields are as follows name, strand, files (because this is an array that can contain one or more pairs, this
                        fields' values are separated by semi-colon and space. Files in a pair are separated by semi-colonpairs are separated by a single space), merge, score, is_ref, exclude_redundant
                        sample_strand takes values 'fr-firststrand', 'fr-unstranded', 'fr-secondstrand' merge, is_ref and exclude_redundant are boolean and take values 'true', 'false' Example: PR1,fr-
                        secondstrand,A_R1.fq;A_R2.fq /samples/paired/B1.fq;/samples/paired/B2.fq,false,2 (default: None)
  --csv_long_samples CSV_LONG_SAMPLES
                        CSV formatted input long read samples. Without headers. The CSV fields are as follows name, strand, files (space separated if there is more than one), quality, score, is_ref,
                        exclude_redundant sample_strand takes values 'fr-firststrand', 'fr-unstranded', 'fr-secondstrand' quality takes values 'low', 'high' is_ref and exclude_redundant are booleans
                        and take values 'true', 'false' Example: Sample1,fr-firststrand,A.fq /samples/long/B.fq ./inputs/C.fq,low,2 (default: None)
  --reference REFERENCE
                        Reference FASTA to annotate (default: None)
  --annotation ANNOTATION
                        Annotation of the reference, this file will be used as the base for the new annotation which will incorporate from the available evidence new gene models or update existing ones
                        (default: None)
  --parameters_file PARAMETERS_FILE
                        Base parameters file, this file can be the output of a previous REAT run which will be used as the base for a new parameters file written to the output_parameters_file argument
                        (default: None)

Mikado:
  Parameters for Mikado runs

  --all_extra_config ALL_EXTRA_CONFIG
                        External configuration file for Paired and Long reads mikado (default: None)
  --long_extra_config LONG_EXTRA_CONFIG
                        External configuration file for Long reads mikado run (default: None)
  --lq_extra_config LQ_EXTRA_CONFIG
                        External configuration file for Low-quality long reads only mikado run (this is only applied when 'separate_mikado_LQ' is enabled) (default: None)
  --all_scoring_file ALL_SCORING_FILE
                        Mikado long and short scoring file (default: None)
  --long_scoring_file LONG_SCORING_FILE
                        Mikado long scoring file (default: None)
  --long_lq_scoring_file LONG_LQ_SCORING_FILE
                        Mikado low-quality long scoring file (default: None)
  --homology_proteins HOMOLOGY_PROTEINS
                        Homology proteins database, used to score transcripts by Mikado (default: None)
  --separate_mikado_LQ SEPARATE_MIKADO_LQ
                        Specify whether or not to analyse low-quality long reads separately from high-quality, this option generates an extra set of mikado analyses including low-quality data (default:
                        None)

Alignment:
  Parameters for alignment of short and long reads

  --short_reads_aligner {hisat,star}
                        Choice of short read aligner (default: hisat)
  --HQ_aligner {minimap2,gmap}
                        Choice of aligner for high-quality long reads (default: gmap)
  --LQ_aligner {minimap2,gmap}
                        Choice of aligner for low-quality long reads (default: minimap2)
  --min_identity MIN_IDENTITY
                        Minimum alignment identity to retain transcript (default: 0.9)
  --min_intron_len MIN_INTRON_LEN
                        Where available, the minimum intron length allowed will be specified for the aligners (default: 20)
  --max_intron_len MAX_INTRON_LEN
                        Where available, the maximum intron length allowed will be specified for the aligners (default: 200000)
  --max_intron_len_ends MAX_INTRON_LEN_ENDS
                        Where available, the maximum *boundary* intron length allowed will be specified for the aligner, when specified this implies max_intron_len only applies to the *internal*
                        introns and this parameter to the *boundary* introns (default: 100000)
  --PR_hisat_extra_parameters PR_HISAT_EXTRA_PARAMETERS
                        Extra command-line parameters for the selected short read aligner, please note that extra parameters are not validated and will have to match the parameters available for the
                        selected read aligner (default: None)
  --PR_star_extra_parameters PR_STAR_EXTRA_PARAMETERS
                        Extra command-line parameters for the selected short read aligner, please note that extra parameters are not validated and will have to match the parameters available for the
                        selected read aligner (default: None)
  --HQ_aligner_extra_parameters HQ_ALIGNER_EXTRA_PARAMETERS
                        Extra command-line parameters for the selected long read aligner, please note that extra parameters are not validated and will have to match the parameters available for the
                        selected read aligner (default: None)
  --LQ_aligner_extra_parameters LQ_ALIGNER_EXTRA_PARAMETERS
                        Extra command-line parameters for the selected long read aligner, please note that extra parameters are not validated and will have to match the parameters available for the
                        selected read aligner (default: None)

Assembly:
  Parameters for assembly of short and long reads

  --skip_scallop
  --HQ_assembler {filter,merge,stringtie,stringtie_collapse}
                        Choice of long read assembler. - filter: Simply filters the reads based on identity and coverage- merge: cluster the input transcripts into loci, discarding "duplicated"
                        transcripts (those with the same exact introns and fully contained or equal boundaries). This option also discards contained transcripts- stringtie: Assembles the long reads
                        alignments into transcripts- stringtie_collapse: Cleans and collapses long reads but does not assemble them (default: filter)
  --LQ_assembler {filter,merge,stringtie,stringtie_collapse}
                        Choice of long read assembler. - filter: Simply filters the reads based on identity and coverage- merge: cluster the input transcripts into loci, discarding "duplicated"
                        transcripts (those with the same exact introns and fully contained or equal boundaries). This option also discards contained transcripts- stringtie: Assembles the long reads
                        alignments into transcripts- stringtie_collapse: Cleans and collapses long reads but does not assembles them (default: stringtie_collapse)
  --HQ_assembler_extra_parameters HQ_ASSEMBLER_EXTRA_PARAMETERS
                        Extra parameters for the long reads assembler, please note that extra parameters are not validated and will have to match the parameters available for the selected assembler
                        (default: None)
  --LQ_assembler_extra_parameters LQ_ASSEMBLER_EXTRA_PARAMETERS
                        Extra parameters for the long reads assembler, please note that extra parameters are not validated and will have to match the parameters available for the selected assembler
                        (default: None)
  --PR_stringtie_extra_parameters PR_STRINGTIE_EXTRA_PARAMETERS
                        Extra parameters for stringtie, please note that extra parameters are not validated and will have to match the parameters available for stringtie (default: None)
  --PR_scallop_extra_parameters PR_SCALLOP_EXTRA_PARAMETERS
                        Extra parameters for scallop, please note that extra parameters are not validated and will have to match the parameters available for scallop (default: None)

Portcullis:
  Parameters specific to portcullis

  --extra_parameters EXTRA_PARAMETERS
                        Extra parameters for portcullis execution (default: None)

ORF Caller:
  Parameters for ORF calling programs

  --orf_caller {prodigal,transdecoder,none}
                        Choice of available orf calling softwares (default: none)
  --orf_calling_proteins ORF_CALLING_PROTEINS
                        Set of proteins to be aligned to the genome for orf prediction by Transdecoder (default: None)
```


#### Homology module

```
usage: reat homology [-h] --genome GENOME --alignment_species ALIGNMENT_SPECIES --annotations_csv ANNOTATIONS_CSV
                     [--annotation_filters {all,none,exon_len,intron_len,internal_stop,aa_len,splicing} [{all,none,exon_len,intron_len,internal_stop,aa_len,splicing} ...]] --mikado_config MIKADO_CONFIG
                     --mikado_scoring MIKADO_SCORING [--junctions JUNCTIONS] [--utrs UTRS] [--pick_extra_config PICK_EXTRA_CONFIG] [--min_cdna_length MIN_CDNA_LENGTH]
                     [--max_intron_length MAX_INTRON_LENGTH] [--filter_min_cds FILTER_MIN_CDS] [--filter_max_intron FILTER_MAX_INTRON] [--filter_min_exon FILTER_MIN_EXON]
                     [--alignment_min_exon_len ALIGNMENT_MIN_EXON_LEN]
                     [--alignment_filters {all,none,exon_len,intron_len,internal_stop,aa_len,splicing} [{all,none,exon_len,intron_len,internal_stop,aa_len,splicing} ...]]
                     [--alignment_min_identity ALIGNMENT_MIN_IDENTITY] [--alignment_min_coverage ALIGNMENT_MIN_COVERAGE] [--alignment_max_per_query ALIGNMENT_MAX_PER_QUERY]
                     [--alignment_show_intron_length] [--exon_f1_filter EXON_F1_FILTER] [--junction_f1_filter JUNCTION_F1_FILTER]

optional arguments:
  -h, --help            show this help message and exit
  --genome GENOME       Fasta file of the genome to annotate (default: None)
  --alignment_species ALIGNMENT_SPECIES
                        Available aligner species, for more information, please look at URL (default: None)
  --annotations_csv ANNOTATIONS_CSV
                        CSV file with reference annotations to extract proteins/cdnas for spliced alignments in csv format. The CSV fields are as follows genome_fasta,annotation_gff e.g
                        Athaliana.fa,Athaliana.gff (default: None)
  --annotation_filters {all,none,exon_len,intron_len,internal_stop,aa_len,splicing} [{all,none,exon_len,intron_len,internal_stop,aa_len,splicing} ...]
                        Filter annotation coding genes by the filter types specified (default: ['none'])
  --mikado_config MIKADO_CONFIG
                        Base configuration for Mikado consolidation stage. (default: None)
  --mikado_scoring MIKADO_SCORING
                        Scoring file for Mikado pick at consolidation stage. (default: None)
  --junctions JUNCTIONS
                        Validated junctions BED file for use in Mikado consolidation stage. (default: None)
  --utrs UTRS           Gene models that may provide UTR extensions to the homology based models at the mikado stage (default: None)
  --pick_extra_config PICK_EXTRA_CONFIG
                        Extra configuration for Mikado pick stage (default: None)
  --min_cdna_length MIN_CDNA_LENGTH
                        Minimum cdna length for models to consider in Mikado consolidation stage (default: 100)
  --max_intron_length MAX_INTRON_LENGTH
                        Maximum intron length for models to consider in Mikado consolidation stage (default: 1000000)
  --filter_min_cds FILTER_MIN_CDS
                        If 'aa_len' filter is enabled for annotation coding features, any CDS smaller thanthis parameter will be filtered out (default: 20)
  --filter_max_intron FILTER_MAX_INTRON
                        If 'intron_len' filter is enabled, any features with introns longer than this parameter will be filtered out (default: 200000)
  --filter_min_exon FILTER_MIN_EXON
                        If 'exon_len' filter is enabled, any features with exons shorter than this parameter will be filtered out (default: 20)
  --alignment_min_exon_len ALIGNMENT_MIN_EXON_LEN
                        Minimum exon length, alignment parameter (default: 20)
  --alignment_filters {all,none,exon_len,intron_len,internal_stop,aa_len,splicing} [{all,none,exon_len,intron_len,internal_stop,aa_len,splicing} ...]
                        Filter alignment results by the filter types specified (default: ['none'])
  --alignment_min_identity ALIGNMENT_MIN_IDENTITY
                        Minimum identity filter for alignments (default: 50)
  --alignment_min_coverage ALIGNMENT_MIN_COVERAGE
                        Minimum coverage filter for alignments (default: 80)
  --alignment_max_per_query ALIGNMENT_MAX_PER_QUERY
                        Maximum number of alignments per input query protein (default: 4)
  --alignment_show_intron_length
                        Add an attribute to the alignment gff with the maximum intron len for each mRNA (default: False)
  --exon_f1_filter EXON_F1_FILTER
                        Filter alignments scored against its original structure with a CDS exon f1 lower than this value (default: None)
  --junction_f1_filter JUNCTION_F1_FILTER
                        Filter alignments scored against its original structure with a CDS junction f1 lower than this value (default: None)
```

