.. REAT documentation master file, created by
   sphinx-quickstart on Fri May  7 16:42:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

REAT - Robust and Extendable eukaryotic Annotation Toolkit
===========================================================

REAT is a robust easy-to-use genome annotation toolkit for turning high-quality genome assemblies into usable and informative resources. REAT makes use of state-of-the-art annotation tools and is robust to varying quality and sources of molecular evidence.

REAT provides an integrated environment that comprises both
a set of workflows geared towards integrating multiple sources of evidence into a genome annotation, and
an execution environment for these workflows.


Installation
============

To install REAT you can:

.. code-block:: bash

  git clone https://github.com/ei-corebioinformatics/reat
  wget https://github.com/broadinstitute/cromwell/releases/download/62/cromwell-62.jar
  conda env create -f reat/reat.yml

These commands will download the cromwell binary required to execute the workflows and make REAT available in the 'reat' conda environment which can be activated using:

.. code-block:: bash

  conda activate reat


Each task in the workflow is configured with default resource requirements appropriate for most tasks, but these can be overriden by user provided ones.
For an example of this file see::

   {
       "ei_annotation.wf_align.long_read_alignment_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.long_read_assembly_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.long_read_indexing_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_alignment_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_alignment_sort_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_merge_resources": {
           "cpu_cores": 4,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_scallop_assembly_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_stringtie_assembly_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_align.short_read_stats_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 8
       },
       "ei_annotation.wf_main_mikado.homology_alignment_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_main_mikado.homology_index_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 8
       },
       "ei_annotation.wf_main_mikado.orf_calling_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 8
       },
       "ei_annotation.wf_main_mikado.protein_alignment_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       },
       "ei_annotation.wf_main_mikado.protein_index_resources":
       {
           "cpu_cores": 6,
           "max_retries": 1,
           "mem_gb": 16
       }
   }


To configure the cromwell engine, there are two relevant files, the cromwell runtime options and the workflow options files.

The cromwell engine can be configured to run in your environment using a file such as:

.. highlight:: none
.. include:: ../cromwell_configuration/slurm.conf
  :literal:


The workflow options can be used to activate the caching behaviour in cromwell, i.e:

.. include:: ../workflow_options/options.json
  :literal:

Running REAT
=============

There are several workflows that make REAT, here we will describe 'transcriptome' and 'homology'.


Transcriptome Workflow
-----------------------

The intention of the transcriptome workflow is to use a variety of data types, from short reads to long reads of varied quality and length.

The data input for the workflow can be defined through the use of comma separated files one for short read samples and another for long read samples. These samples are then processed in several steps, first they are aligned to the genome, then assembled into transcripts, junctions are determined from the data and finally they are combined into a consolidated set of gene models.

The aligner and assembly programs used for short and long read samples can be selected through command line arguments. There are also command line arguments to select extra options to be applied at each step.

In case an annotation is available, this can be provided for junctions and reference models to be extracted and these can then be augmented using the evidence present in the data.

.. highlight:: none
.. include:: ./transcriptome_help.txt
  :literal:

Sample files
^^^^^^^^^^^^^

The way samples are organised in the input files reflects how the files that correspond to the sample will be processed.
Data can be combined or kept separate at different stages of the workflow in accordance with the configuration provided
and the characteristics of the data.


Short read data
""""""""""""""""

Each line corresponds to a sample.
There are four required fields: Sample name, strandness, RNA-seq paired data, merge.
Followed by three optional fields: score, is_reference, exclude_redundant.
Previous fields to an optional field must be present in the line.
Files within a pair are separated by semi-colon and where there are multiple pairs in a sample, these are separated by spaces.

.. code-block:: bash

  Ara0,fr-firststrand,data/Ara1.1.fastq.gz;data/Ara1.2.fastq.gz,true,20
  Ara1,fr-firststrand,data/Ara1.1.fastq.gz;data/Ara1.2.fastq.gz data/Ara2.1.fastq.gz;data/Ara2.2.fastq.gz,true,20
  Ara2,fr-firststrand,data/Ara3.1.fastq.gz;data/Ara3.2.fastq.gz data/Ara5.1.fastq.gz;data/Ara5.2.fastq.gz data/Ara6.1.fastq.gz;data/Ara6.2.fastq.gz,false


Sample RNA-seq data can be merged in different places, the options for controlling when the merging happens are as follows:
All transcripts assembled from paired reads within a sample are combined after assembling, paired read alignments can be merged before assembly using the 'merge' parameter in the CSV file.

Junctions
++++++++++

Junctions from RNA-seq data can be determined in several ways.
By default junctions are collected for all the RNA-seq fastq pair as defined in the 'RNA-seq paired data' section of the CSV file for each sample.
Alternatively, samples can be combined where appropriate using the 'ei_annotation.wf_align.group_to_samples' parameter in the input.json file.
This parameter will define arbitrary groupings of the samples present in the short read CSV, with the following format::

  "ei_annotation.wf_align.group_to_samples": {
    "group1": ["Sample1", "Sample2"],
    "group2": ["Sample3", "Sample4"]
  }

These groups will be validated against the samples in the CSV files, group names should be unique, samples can only belong to a single group and all samples should be part of a group.

Long read data
"""""""""""""""

Each line corresponds to a sample.
There are four required fields: Sample name, strandness, RNA-seq long read data, merge.
Followed by three optional fields: score, is_reference and exclude_redundant.
Previous fields to an optional field must be present in the line.
Where multiple read files correspond to a single sample (this implies they result in a single set of transcripts), the third column will contain all the files separated by spaces.

.. code-block:: bash

   A01_1,fr-firststrand,data/A1_1.fastq.gz,low
   A01_2,fr-firststrand,data/A1_2.fastq.gz,low
   B01,fr-firststrand,data/B1.fastq.gz,low,10,true,true
   C01,fr-firststrand,data/C1.fastq.gz,low
   ALL,fr-firststrand,data/D1_1.fastq.gz data/D1_2.fastq.gz data/D1_3.fastq.gz data/D1_4.fastq.gz,low
   CCS,fr-firststrand,data/CCS.fastq.gz,high
   polished,fr-firststrand,data/polished.fastq.gz,high


.. warning::

   The 'reference' sample name is reserved for internal use.
   If this name is being used in any of the sample input CSV files, you will be notified with an error message.


Homology workflow
-----------------

When there is protein evidence available from related species, the homology workflow can be used to generate gene models based on this evidence.
This is achieved by aligning the proteins provided through a set of related species annotations and evaluating these alignments to generate a score.
Protein alignments are evaluated in two ways: Coherence of the alignment structure with respect to the original model's structure and consensus structure from the multiple species.
These scores are then used by Mikado to group and filter models, generating a set of predicted models.

.. highlight:: none
.. include:: ./homology_help.txt
  :literal:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
