.. REAT documentation master file, created by
   sphinx-quickstart on Fri May  7 16:42:55 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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


Each task in the workflow is configured with default resource requirements appropriate for most tasks, but these can be
overriden by user provided ones. For examples of resource configuration files, refer to each module's description.

To configure the cromwell engine, there are two relevant files, the cromwell runtime options and the workflow options files.

The cromwell engine can be configured to run in your environment using a file such as:

.. highlight:: none
.. include:: ../cromwell_configuration/slurm.conf
  :literal:


The workflow options can be used to activate the caching behaviour in cromwell, i.e:

.. include:: ../workflow_options/options.json
  :literal:


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules/transcriptome/index
   modules/homology/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
