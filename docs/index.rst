.. REAT documentation master file, created by
   sphinx-quickstart on Fri May  7 16:42:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to REAT's documentation!
================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

REAT - Robust and Extendable eukaryotic Annotation Toolkit
===========================================================

REAT is a robust easy-to-use genome annotation toolkit for turning high-quality genome assemblies into usable and informative resources. REAT makes use of state-of-the-art annotation tools and is robust to varying quality and sources of molecular evidence.

REAT provides an integrated environment that comprises both

a set of workflows geared towards integrating multiple sources of evidence into a genome annotation, and
an execution environment for these workflows.



Running REAT
=============

There are several workflows that make REAT, here we will describe 'transcriptome' and 'homology'.


Transcriptome Workflow
-----------------------

The intention of the transcriptome workflow is to use a variety of data types, from short reads to long reads of varied quality and length.

The data input for the workflow can be defined through the use of comma separated files one for short read samples and another for long read samples. These samples are then processed in several steps, first they are aligned to the genome, then assembled into transcripts, junctions are determined from the data and finally they are combined into a consolidated set of gene models.

The aligner and assembly programs used for short and long read samples can be selected through command line arguments. There are also command line arguments to select extra options to be applied at each step.

In case an annotation is available, this can be provided for junctions and reference models to be extracted and these can then be augmented using the evidence present in the data.

Sample files
^^^^^^^^^^^^^

The way samples are organised in the input files reflects how the files that correspond to the sample will be processed. Data can be combined or kept separate at different stages of the workflow in accordance with the configuration provided and the characteristics of the data.


Short read data
""""""""""""""""

::

   Ara0,fr-firststrand,data/Ara1.1.fastq.gz;data/Ara1.2.fastq.gz,true,20
   Ara1,fr-firststrand,data/Ara1.1.fastq.gz;data/Ara1.2.fastq.gz data/Ara2.1.fastq.gz;data/Ara2.2.fastq.gz,true,20
   Ara2,fr-firststrand,data/Ara3.1.fastq.gz;data/Ara3.2.fastq.gz data/Ara5.1.fastq.gz;data/Ara5.2.fastq.gz data/Ara6.1.fastq.gz;data/Ara6.2.fastq.gz,false



Long read data
"""""""""""""""

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
