Transcriptome Workflow
========================

The intention of the transcriptome workflow is to use a variety of data types, from short reads to long reads of varied quality and length.

The data input for the workflow can be defined through the use of comma separated files one for short read samples and another for long read samples. These samples are then processed in several steps, first they are aligned to the genome, then assembled into transcripts, junctions are determined from the data and finally they are combined into a consolidated set of gene models.

The aligner and assembly programs used for short and long read samples can be selected through command line arguments. There are also command line arguments to select extra options to be applied at each step.

In case an annotation is available, this can be provided for junctions and reference models to be extracted and these can then be augmented using the evidence present in the data.

.. highlight:: none
.. include:: ./transcriptome_help.txt
  :literal:

Sample files
------------------

The way samples are organised in the input files reflects how the files that correspond to the sample will be processed.
Data can be combined or kept separate at different stages of the workflow in accordance with the configuration provided
and the characteristics of the data.


Short read data
^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^

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

.. image:: /_static/transcriptome_workflow.png
  :alt: Transcriptome workflow diagram

Configurable computational resources available::

  "ei_annotation.wf_align.long_read_alignment_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.long_read_assembly_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.long_read_indexing_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.long_read_twopass_merge_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.long_read_twopass_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.portcullis_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.sanitise_reference_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_alignment_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_alignment_sort_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_indexing_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_merge_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_scallop_assembly_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_stats_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_align.short_read_stringtie_assembly_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.homology_alignment_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.homology_index_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_all_pick_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_all_prepare_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_all_serialise_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_lq_pick_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_lq_prepare_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_lq_serialise_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_pick_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_prepare_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.mikado_long_serialise_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.orf_calling_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.protein_alignment_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)",
  "ei_annotation.wf_main_mikado.protein_index_resources": " {
                 cpu_cores -> Int?
                max_retries -> Int?
                boot_disk_gb -> Int?
                queue -> String?
                disk_gb -> Int?
                constraints -> String?
                mem_gb -> Float?
                preemptible_tries -> Int?
                }? (optional)"

