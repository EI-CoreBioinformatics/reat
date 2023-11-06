
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

For ease of development, Singularity is recommended. In case this is not available (you are working on MacOS and haven't got https://sylabs.io/guides/3.5/admin-guide/installation.html#mac available), please follow the instructions in the `REAT/reat_singularity.def` for building and installing the binary dependencies.

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

Please install `reat` to ensure the scripts are findable during execution. To install, simply use from your current pip/conda environment:

```shell
pip install ./reat
```

## Running REAT

Cromwell supports a range of execution backends such as cloud services and batch schedulers. 
These can be configured in files such as those present in the `cromwell_configuration` directory.
For more examples you can also visit https://github.com/broadinstitute/cromwell/tree/develop/cromwell.example.backends

### Local

REAT finds and executes a binary named `cromwell` which can be configured by placing a file similar to:

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

## Workflow

### Transcriptome Workflow
![Alt text](/docs/_static/REAT_Transcriptome.png)

### Homology Workflow
![Alt text](/docs/_static/REAT_Homology.png)

### Prediction Workflow
![Alt text](/docs/_static/REAT_Prediction.png)

## Releasing new versions

To create a new release, make sure your repository is in a clean state and please use:

```shell
bumpversion minor --verbose --tag --message "Minor version bump" --tag-message "Minor version bump

List of changes (note that this message will be used in the Release)
Finally you may want to update the release message by editing it directly
```

The above will bump a minor version (options are major, minor and patch), create a tag with the message above, trigger the github actions for checking the commit and generating a new github release.
