#!/usr/bin/env python3

# __main__.py can result in the following outcomes
#     Success:
#         Reat is started in a 'run' backend (local, hpc)
#         Reat is started in a 'server' hosted on the network
#     Failure:
#         Reat fails to start or be submitted and the reason is provided to the user

# __main__.py checks the user environment to ensure the required software is available if this is not the case
# the user is informed of which software is and which isn't available.

# __main__.py parses all the arguments required, generates an input file for cromwell and submits a job to the requested
# backend (run or server mode). The arguments are validated using the json input schema defined in the validation
# directory.

# __main__.py supports multiple input json files which are subsequently merged into a single file, allows reusability
# of some parts of the input such as resource requirements and default extra parameters.

# __main__.py accepts many customisation parameters in the command-line which are translated into the inputs.json before
# inputs.json is validated.

# If the user wishes to 'cancel' the workflow, SIGTERM or SIGINT will be managed by cascading them to the cromwell
# process, SIGINT should allow for 'happy' process termination.

# Generates all required inputs and parses mikado's multiple runs extra configurations and places them in a reusable
# location where they can be edited for the user's convenience
import datetime
import json
import os
import subprocess
import sys
import textwrap
import time
from textwrap import wrap

from annotation import UTR_SELECTION_OPTIONS, LONG_READ_ALIGNER_CHOICES
from annotation import VERSION
from annotation.homology import homology_module
from annotation.prediction import prediction_module
from annotation.prediction_module import add_classification_parser_parameters
from annotation.transcriptome import (
    transcriptome_cli_validation,
    genetic_code_str_to_int,
)
from annotation.transcriptome import transcriptome_module

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from argparse import (
    ArgumentDefaultsHelpFormatter,
    ArgumentParser,
    FileType,
    RawTextHelpFormatter,
)


class ReatHelpFormatter(ArgumentDefaultsHelpFormatter, RawTextHelpFormatter):
    pass


def check_environment(force_quit=True):
    """
    Check user's environment for required packages and their versions. Writing an to stderr if any of the software
    dependencies version doesn't match with the expected ones.
    :param force_quit: If any of the software versions doesn't match REAT will exit
    :return:
    """
    software_available = {
        "spaln": {"command": ["spaln"], "result": "SPALN version 2.4.0"},
        "sortgrcd": {
            "command": ["sortgrcd", "--help"],
            "result": "sortgrcd version 2.2",
        },
        "mikado": {
            "command": "mikado --version".split(" "),
            "result": "Mikado v2.0rc2",
        },
        "diamond": {
            "command": "diamond version".split(" "),
            "result": "diamond version 0.9.31",
        },
        "blastn": {"command": "blastn -version".split(" "), "result": "blastn: 2.7.1+"},
        "blastx": {"command": "blastx -version".split(" "), "result": "blastx: 2.7.1+"},
        "samtools": {
            "command": "samtools --version".split(" "),
            "result": "samtools 1.9",
        },
        "gffread": {"command": "gffread --version".split(" "), "result": "0.12.2"},
        "gmap": {
            "command": "gmap --version".split(" "),
            "result": "version 2019-02-15",
        },
        "minimap2": {"command": "minimap2 --version".split(" "), "result": "2.17-r941"},
        "hisat2": {"command": "hisat2 --version".split(" "), "result": "version 2.1.0"},
        "star": {"command": "STAR --version".split(" "), "result": "2.7.3a"},
        "seqtk": {"command": ["seqtk"], "result": "Version: 1.3-r116-dirty"},
        "stringtie": {"command": "stringtie --version".split(" "), "result": "2.1.1"},
        "scallop": {"command": "scallop --version".split(" "), "result": "v0.10.4"},
        "scallop-lr": {
            "command": "scallop-lr --version".split(" "),
            "result": "v0.9.2",
        },
        "prodigal": {
            "command": "prodigal -v".split(" "),
            "result": "Prodigal V2.6.3: February, 2016",
        },
        "transdecoder": {
            "command": "TransDecoder.LongOrfs --version".split(" "),
            "result": "TransDecoder.LongOrfs 5.5.0",
        },
        "portcullis": {
            "command": "portcullis --version".split(" "),
            "result": "portcullis 1.2.0",
        },
        "junctools": {"command": "junctools --version".split(" "), "result": "1.2.0"},
    }

    programs_not_found = set()
    for key, item in software_available.items():
        try:
            result = subprocess.run(item["command"], capture_output=True)
        except FileNotFoundError:
            programs_not_found.add(key)
            continue
        except RuntimeError as e:
            print(
                f"When executing {key}, found the following error:\n{e}\n"
                f"Please ensure your environment and hardware support REAT's requirements",
                file=sys.stderr,
            )
        except KeyboardInterrupt:
            exit(10)
        output = result.stdout.decode()
        output += result.stderr.decode()
        item["rc"] = result.returncode
        if key == "seqtk":
            if "Version" in output:
                item["rc"] = 0
        elif key == "spaln":
            if "No input seq file !" in output:
                item["rc"] = 0
        elif key == "sortgrcd":
            if output != "":
                item["rc"] = 0
        if item["result"] not in output:
            print('"', key, '"', " version information:", sep="", file=sys.stderr)
            print('"""', file=sys.stderr)
            print(output.strip(), file=sys.stderr)
            print('"""', file=sys.stderr)
            print(file=sys.stderr)
            print(
                "Does not contain indication of the required version:", file=sys.stderr
            )
            print('"""', file=sys.stderr)
            print(item["result"], file=sys.stderr)
            print('"""', file=sys.stderr)
            print(file=sys.stderr)
            print(file=sys.stderr)
        # Command not in path, wrong version or failed to execute
        if item["rc"] != 0:
            raise FileNotFoundError(
                "Command {0} is missing, please check your PATH".format(key)
            )
    if len(programs_not_found) > 0:
        print(
            f"When checking the environment, the software {', '.join([p for p in programs_not_found])}, not found.\n"
            f"Please make sure it is in the PATH environment variable of the shell executing REAT.\n\n"
            f"Currently PATH contains the following:"
        )
        print("\n".join(wrap(", ".join(os.environ["PATH"].split(os.pathsep)))))
        if force_quit:
            sys.exit(2)
    return software_available


def parse_arguments():
    """
    Parses and validates REAT's CLI arguments, the values inside defined on the CLI files have not yet been validated.

    :return: Object containing the validated CLI input arguments.
    """
    reat_ap = ArgumentParser(add_help=True, formatter_class=ReatHelpFormatter)

    reat_ap.add_argument(
        "-j", "--jar_cromwell", type=FileType("r"), help="Cromwell server jar file"
    )
    reat_ap.add_argument(
        "-r",
        "--runtime_configuration",
        type=FileType("r"),
        help="Configuration file for the backend, please follow "
        "https://cromwell.readthedocs.io/en/stable/backends/HPC/ for more information.\n"
        "An example of this file can be found at ",
    )

    reat_ap.add_argument(
        "-c",
        "--computational_resources",
        type=FileType("r"),
        help="Computational resources for REAT, please look at the template for more information",
        required=True,
    )
    reat_ap.add_argument(
        "-o",
        "--output_parameters_file",
        type=str,
        help="REAT parameters file, this file will be used as the input for REAT. "
        "It provides the arguments for the workflow runtime.",
        default="reat_input.json",
    )
    reat_ap.add_argument(
        "--workflow_options_file",
        type=FileType("r"),
        help="Workflow execution options, includes cache usage and result directories "
        "structure and location",
    )

    subparsers = reat_ap.add_subparsers(help="sub-command help", dest="reat_module")

    transcriptome_ap = subparsers.add_parser(
        "transcriptome", help="Transcriptome module", formatter_class=ReatHelpFormatter
    )
    # General inputs
    transcriptome_ap.add_argument(
        "--reference",
        type=FileType("r"),
        help="Reference FASTA to annotate",
        required=True,
    )
    transcriptome_ap.add_argument(
        "--samples",
        nargs="+",
        type=FileType("r"),
        help="Reads organised in the input specification for REAT, for more information "
        "please look at https://github.com/ei-corebioinformatics/reat\nfor an example",
    )
    transcriptome_ap.add_argument(
        "--csv_paired_samples",
        type=FileType("r"),
        help=textwrap.dedent(
            """\
CSV formatted input paired read samples. Without headers.

The CSV fields are as follows name, strand, files (because this is an array that can contain one or more pairs, this fields' values are separated by semi-colon and space. Files in a pair are separated by semi-colon pairs are separated by a single space), merge, score, is_ref, exclude_redundant.

sample_strand takes values \'fr-firststrand\', \'fr-unstranded\', \'fr-secondstrand\'

merge, is_ref and exclude_redundant are boolean and take values 'true', 'false'

Example:
PR1,fr-secondstrand,A_R1.fq;A_R2.fq /samples/paired/B1.fq;/samples/paired/B2.fq,false,2
"""
        ),
    )
    transcriptome_ap.add_argument(
        "--csv_long_samples",
        type=FileType("r"),
        help=textwrap.dedent(
            """\
CSV formatted input long read samples. Without headers."
The CSV fields are as follows name, strand, files (space separated if there is more than one), quality, score, is_ref, exclude_redundant

sample_strand takes values \'fr-firststrand\', \'fr-unstranded\', \'fr-secondstrand\'
quality takes values 'low', 'high'
is_ref and exclude_redundant are booleans and take values 'true', 'false'

Example:

Sample1,fr-firststrand,A.fq /samples/long/B.fq ./inputs/C.fq,low,2"""
        ),
    )
    transcriptome_ap.add_argument(
        "--annotation",
        type=FileType("r"),
        help=textwrap.dedent(
            """\
Annotation of the reference, this file will be used as the base for the new annotation which will incorporate from the available evidence new gene models or update existing ones"""
        ),
    )
    transcriptome_ap.add_argument(
        "--annotation_score",
        type=int,
        default=1,
        help="Score for models in the reference annotation file",
    )
    transcriptome_ap.add_argument(
        "--check_reference",
        action="store_true",
        default=False,
        help="At mikado stage, annotation models will be evaluated in the same manner as "
        "RNA-seq based models, removing any models\ndeemed incorrect",
    )
    transcriptome_ap.add_argument(
        "--mode",
        choices=["basic", "update", "only_update"],
        default="basic",
        help="basic: Annotation models are treated the same as the RNA-Seq models at the pick"
        " stage.\n"
        "update: Annotation models are prioritised but also novel loci are reported.\n"
        "only_update: Annotation models are prioritised and non-reference loci are "
        "excluded.",
    )
    transcriptome_ap.add_argument(
        "--extra_junctions",
        type=FileType("r"),
        help="Extra junctions provided by the user, this file will be used as a set of valid"
        " junctions for alignment of short and\nlong read samples, in the case of long"
        " reads, these junctions are combined with the results of portcullis whenever\n"
        " short read samples have been provided as part of the input datasets",
    )
    transcriptome_ap.add_argument(
        "--skip_mikado_long",
        action="store_true",
        default=False,
        help="Disables generation of the long read only mikado run",
    )
    transcriptome_ap.add_argument(
        "--filter_HQ_assemblies",
        action="store_true",
        default=False,
        help="Use all the junctions available to filter the HQ_assemblies before mikado",
    )
    transcriptome_ap.add_argument(
        "--filter_LQ_assemblies",
        action="store_true",
        default=False,
        help="Use all the junctions available to filter the LQ_assemblies before mikado",
    )
    transcriptome_ap.add_argument(
        "--parameters_file",
        type=FileType("r"),
        help="Base parameters file, this file can be the output of a previous REAT run "
        "which will be used as the base for a new\nparameters file written to the"
        " output_parameters_file argument",
    )
    transcriptome_ap.add_argument(
        "--genetic_code",
        help="\n".join(
            textwrap.wrap(
                "Parameter for the translation table used in Mikado for translating CDS "
                "sequences, and for ORF calling, can take values in the genetic code range of "
                "NCBI as an integer. E.g 1, 6, 10 or when using TransDecoder as ORF caller, "
                "one of: {}. 0 is equivalent to Standard, NCBI #1, but only ATG is "
                "considered a valid start codon.".format(
                    ", ".join(genetic_code_str_to_int.keys())
                ),
                110,
            )
        ),
        default="0",
    )

    # Mikado arguments
    mikado_parameters = transcriptome_ap.add_argument_group(
        "Mikado", "Parameters for Mikado runs"
    )
    mikado_parameters.add_argument(
        "--all_extra_config",
        type=FileType("r"),
        help="External configuration file for Paired and Long reads mikado",
    )
    mikado_parameters.add_argument(
        "--long_extra_config",
        type=FileType("r"),
        help="External configuration file for Long reads mikado run",
    )
    mikado_parameters.add_argument(
        "--lq_extra_config",
        type=FileType("r"),
        help="External configuration file for Low-quality long reads only mikado run "
        "(this is only applied when \n'separate_mikado_LQ' is enabled)",
    )
    mikado_parameters.add_argument(
        "--all_scoring_file",
        type=FileType("r"),
        help="Mikado long and short scoring file",
        required=True,
    )
    mikado_parameters.add_argument(
        "--long_scoring_file", type=FileType("r"), help="Mikado long scoring file"
    )
    mikado_parameters.add_argument(
        "--long_lq_scoring_file",
        type=FileType("r"),
        help="Mikado low-quality long scoring file",
    )
    mikado_parameters.add_argument(
        "--homology_proteins",
        type=FileType("r"),
        help="Homology proteins database, used to score transcripts by Mikado",
    )
    mikado_parameters.add_argument(
        "--separate_mikado_LQ",
        type=bool,
        help="Specify whether or not to analyse low-quality long reads separately from "
        "high-quality, this option generates an\nextra set of mikado analyses "
        "including low-quality data",
    )
    mikado_parameters.add_argument(
        "--exclude_LQ_junctions",
        action="store_true",
        default=False,
        help="When this parameter is defined, junctions derived from low-quality long reads "
        "will not be included in the set of\nvalid junctions for the mikado analyses",
    )

    # Aligner choices
    alignment_parameters = transcriptome_ap.add_argument_group(
        "Alignment", "Parameters for alignment of short and long reads"
    )
    alignment_parameters.add_argument(
        "--short_reads_aligner",
        choices=["hisat", "star"],
        help="Choice of short read aligner",
        default="hisat",
    )
    alignment_parameters.add_argument(
        "--skip_2pass_alignment",
        action="store_true",
        default=False,
        help="If not required, the second round of alignments for 2passtools can be "
        "skipped when this parameter\nis active",
    )
    alignment_parameters.add_argument(
        "--HQ_aligner",
        choices=LONG_READ_ALIGNER_CHOICES,
        help="Choice of aligner for high-quality long reads",
        default="minimap2",
    )
    alignment_parameters.add_argument(
        "--LQ_aligner",
        choices=LONG_READ_ALIGNER_CHOICES,
        help="Choice of aligner for low-quality long reads",
        default="minimap2",
    )
    alignment_parameters.add_argument(
        "--min_identity",
        type=int,
        choices=range(0, 101),
        metavar="[0-100]",
        help="Minimum alignment identity (passed only to gmap)",
        default=90,
    )
    alignment_parameters.add_argument(
        "--min_intron_len",
        type=int,
        help="Where available, the minimum intron length allowed will be specified for "
        "the aligners",
        default=20,
    )
    alignment_parameters.add_argument(
        "--max_intron_len",
        type=int,
        help="Where available, the maximum intron length allowed will be specified for "
        "the aligners",
        default=200000,
    )
    alignment_parameters.add_argument(
        "--max_intron_len_ends",
        type=int,
        help="Where available, the maximum *boundary* intron length allowed will be "
        "specified for the aligner, when specified\nthis implies max_intron_len "
        "only applies to the *internal* introns and this parameter to the *boundary*"
        " introns",
        default=100000,
    )

    alignment_parameters.add_argument(
        "--PR_hisat_extra_parameters",
        type=str,
        help="Extra command-line parameters for the selected short read aligner, please "
        "note that extra parameters are not\nvalidated and will have to match the "
        "parameters available for the selected read aligner",
    )
    alignment_parameters.add_argument(
        "--PR_star_extra_parameters",
        type=str,
        help="Extra command-line parameters for the selected short read aligner, please "
        "note that extra parameters are not\nvalidated and will have to match the "
        "parameters available for the selected read aligner",
    )
    alignment_parameters.add_argument(
        "--HQ_aligner_extra_parameters",
        type=str,
        help="Extra command-line parameters for the selected long read aligner, please "
        "note that extra parameters are not\nvalidated and will have to match the "
        "parameters available for the selected read aligner",
    )
    alignment_parameters.add_argument(
        "--LQ_aligner_extra_parameters",
        type=str,
        help="Extra command-line parameters for the selected long read aligner, please "
        "note that extra parameters are not\nvalidated and will have to match the "
        "parameters available for the selected read aligner",
    )

    # Assembler choices
    assembly_parameters = transcriptome_ap.add_argument_group(
        "Assembly", "Parameters for assembly of short and long reads"
    )
    assembly_parameters.add_argument(
        "--skip_scallop", action="store_true", default=False
    )
    assembly_parameters.add_argument(
        "--HQ_assembler",
        choices=["filter", "merge", "stringtie", "stringtie_collapse"],
        help="Choice of long read assembler."
        "\n- filter: Simply filters the reads based on identity and coverage"
        "\n- merge: cluster the input transcripts into loci, discarding "
        '"duplicated" transcripts (those with the same exact\n\tintrons and fully '
        "contained or equal boundaries). This option also discards contained "
        "transcripts"
        "\n- stringtie: Assembles the long reads alignments into transcripts"
        "\n- stringtie_collapse: Cleans and collapses long reads but does not "
        "assemble them",
        default="filter",
    )
    assembly_parameters.add_argument(
        "--LQ_assembler",
        choices=["filter", "merge", "stringtie", "stringtie_collapse"],
        help="Choice of long read assembler."
        "\n- filter: Simply filters the reads based on identity and coverage"
        "\n- merge: cluster the input transcripts into loci, discarding "
        '"duplicated" transcripts (those with the same exact\n\tintrons and fully '
        "contained or equal boundaries). This option also discards contained "
        "transcripts"
        "\n- stringtie: Assembles the long reads alignments into transcripts"
        "\n- stringtie_collapse: Cleans and collapses long reads but does not "
        "assembles them",
        default="stringtie_collapse",
    )
    assembly_parameters.add_argument(
        "--HQ_min_identity",
        type=int,
        choices=range(0, 101),
        metavar="[0-100]",
        help="When the 'filter' option is selected, this parameter defines the minimum "
        "identity used to filtering",
    )
    assembly_parameters.add_argument(
        "--HQ_min_coverage",
        type=int,
        choices=range(0, 101),
        metavar="[0-100]",
        help="When the 'filter' option is selected, this parameter defines the minimum "
        "coverage used for filtering",
    )
    assembly_parameters.add_argument(
        "--HQ_assembler_extra_parameters",
        help="Extra parameters for the long reads assembler, please note that extra "
        "parameters are not validated and will have to\nmatch the parameters "
        "available for the selected assembler",
    )
    assembly_parameters.add_argument(
        "--LQ_min_identity",
        type=int,
        choices=range(0, 101),
        metavar="[0-100]",
        help="When the 'filter' option is selected, this parameter defines the minimum "
        "identity used to filtering",
    )
    assembly_parameters.add_argument(
        "--LQ_min_coverage",
        type=int,
        choices=range(0, 101),
        metavar="[0-100]",
        help="When the 'filter' option is selected, this parameter defines the minimum "
        "coverage used for filtering",
    )
    assembly_parameters.add_argument(
        "--LQ_assembler_extra_parameters",
        help="Extra parameters for the long reads assembler, please note that extra "
        "parameters are not validated and will have to\nmatch the parameters "
        "available for the selected assembler",
    )
    assembly_parameters.add_argument(
        "--PR_stringtie_extra_parameters",
        help="Extra parameters for stringtie, please note that extra "
        "parameters are not validated and will have to\nmatch the parameters "
        "available for stringtie",
    )
    assembly_parameters.add_argument(
        "--PR_scallop_extra_parameters",
        help="Extra parameters for scallop, please note that extra "
        "parameters are not validated and will have to\nmatch the parameters "
        "available for scallop",
    )

    # Portcullis extra parameters
    portcullis_parameters = transcriptome_ap.add_argument_group(
        "Portcullis", "Parameters specific to portcullis"
    )
    portcullis_parameters.add_argument(
        "--extra_parameters", type=str, help="Extra parameters for portcullis execution"
    )

    # Orf calling
    orf_calling_parameters = transcriptome_ap.add_argument_group(
        "ORF Caller", "Parameters for ORF calling programs"
    )
    orf_calling_parameters.add_argument(
        "--orf_caller",
        choices=["prodigal", "transdecoder", "none"],
        help="Choice of available orf calling softwares",
        default="prodigal",
    )
    orf_calling_parameters.add_argument(
        "--orf_calling_proteins",
        type=FileType("r"),
        help="Set of proteins to be aligned to the genome for orf prediction by "
        "Transdecoder",
    )

    homology_ap = subparsers.add_parser(
        "homology", help="Homology module", formatter_class=ReatHelpFormatter
    )

    homology_ap.add_argument(
        "--genome",
        type=FileType("r"),
        help="Fasta file of the genome to annotate",
        required=True,
    )
    homology_ap.add_argument(
        "-p",
        "--output_prefix",
        type=str,
        default="xspecies",
        help="Prefix for the final output files",
    )
    homology_ap.add_argument(
        "--alignment_species",
        type=str,
        help="Species specific parameters, select a value from the first or second column of "
        "https://raw.githubusercontent.com/ogotoh/spaln/master/table/gnm2tab",
        required=True,
    )
    homology_ap.add_argument(
        "--codon_table", type=int, default=1, help="NCBI based codon translation table"
    )
    homology_ap.add_argument(
        "--annotations_csv",
        type=FileType("r"),
        help="CSV file with reference annotations to extract proteins/cdnas for spliced alignments"
        ". The CSV fields are: genome_fasta,annotation_gff"
        "\nExample:\nAthaliana.fa,Athaliana.gff",
    )
    homology_ap.add_argument(
        "--protein_sequences",
        type=str,
        nargs="*",
        help="List of files containing protein sequences to use as evidence",
    )
    homology_ap.add_argument(
        "--annotation_filters",
        choices=[
            "all",
            "none",
            "exon_len",
            "intron_len",
            "internal_stop",
            "aa_len",
            "splicing",
        ],
        nargs="+",
        help="Filter annotation coding genes by the filter types specified",
        default=["none"],
    )
    homology_ap.add_argument(
        "--mikado_config",
        type=FileType("r"),
        help="Base configuration for Mikado consolidation stage.",
        required=True,
    )
    homology_ap.add_argument(
        "--mikado_scoring",
        type=FileType("r"),
        help="Scoring file for Mikado pick at consolidation stage.",
        required=True,
    )
    homology_ap.add_argument(
        "--junctions",
        type=FileType("r"),
        help="Validated junctions BED file for use in Mikado consolidation stage.",
    )
    homology_ap.add_argument(
        "--utrs",
        type=FileType("r"),
        help="Gene models that may provide UTR extensions to the homology based models at the "
        "mikado stage",
    )
    homology_ap.add_argument(
        "--pick_extra_config",
        type=FileType("r"),
        help="Extra configuration for Mikado pick stage",
    )
    homology_ap.add_argument(
        "--min_cdna_length",
        type=int,
        default=100,
        help="Minimum cdna length for models to consider in Mikado consolidation stage",
    )
    homology_ap.add_argument(
        "--max_intron_length",
        type=int,
        default=1000000,
        help="Maximum intron length for models to consider in Mikado consolidation stage",
    )
    homology_ap.add_argument(
        "--filter_min_cds",
        type=int,
        help="If 'aa_len' filter is enabled for annotation coding features, any CDS smaller than"
        "this parameter will be filtered out",
        default=20,
    )
    homology_ap.add_argument(
        "--filter_max_intron",
        type=int,
        help="If 'intron_len' filter is enabled, any features with introns longer than this "
        "parameter will be filtered out",
        default=200000,
    )
    homology_ap.add_argument(
        "--filter_min_exon",
        type=int,
        help="If 'exon_len' filter is enabled, any features with exons shorter than this "
        "parameter will be filtered out",
        default=20,
    )
    homology_ap.add_argument(
        "--alignment_min_exon_len",
        type=int,
        help="Minimum exon length, alignment parameter",
        default=20,
    )
    homology_ap.add_argument(
        "--alignment_filters",
        choices=[
            "all",
            "none",
            "exon_len",
            "intron_len",
            "internal_stop",
            "aa_len",
            "splicing",
        ],
        help="Filter alignment results by the filter types specified",
        nargs="+",
        default=["none"],
    )
    homology_ap.add_argument(
        "--alignment_min_identity",
        type=int,
        help="Minimum identity filter for alignments",
        default=50,
    )
    homology_ap.add_argument(
        "--alignment_min_coverage",
        type=int,
        help="Minimum coverage filter for alignments",
        default=80,
    )
    homology_ap.add_argument(
        "--alignment_max_per_query",
        type=int,
        default=4,
        help="Maximum number of alignments per input query protein",
    )
    homology_ap.add_argument(
        "--alignment_recursion_level",
        type=int,
        default=6,
        help="SPALN's Q value, indicating the level of recursion for the Hirschberg algorithm",
    )
    homology_ap.add_argument(
        "--alignment_show_intron_length",
        action="store_true",
        help="Add an attribute to the alignment gff with the maximum intron len for each mRNA",
    )
    homology_ap.add_argument(
        "--exon_f1_filter",
        type=int,
        help="Filter alignments scored against its original structure with a CDS exon f1 "
        "lower than this value",
    )
    homology_ap.add_argument(
        "--junction_f1_filter",
        type=int,
        help="Filter alignments scored against its original structure with a CDS junction f1 "
        "lower than this value",
    )

    prediction_ap = subparsers.add_parser(
        "prediction", help="Prediction module", formatter_class=ReatHelpFormatter
    )

    prediction_ap.add_argument(
        "--genome", type=FileType("r"), required=True, help="Genome fasta file"
    )
    prediction_ap.add_argument(
        "--augustus_config_path",
        type=str,
        required=True,
        help="Template path for augustus config, this path will not be modified as a copy will "
        "be created internally for the workflow's use",
    )
    prediction_ap.add_argument(
        "--extrinsic_config",
        type=FileType("r"),
        help="Augustus extrinsic configuration file, defines the boni/mali for each type of "
        "feature-evidence combination",
    )
    prediction_ap.add_argument(
        "--species",
        type=str,
        required=True,
        help="Name of the species to train models for, if it does not exist in the augustus "
        "config path it will be created.",
    )
    prediction_ap.add_argument(
        "--codon_table", type=int, default=1, help="NCBI based codon translation table"
    )
    prediction_ap.add_argument(
        "--chunk_size",
        type=int,
        default=5000000,
        help="Maximum length of sequence to be processed by Augustus or EVM",
    )
    prediction_ap.add_argument(
        "--overlap_size",
        type=int,
        default=500000,
        help="Overlap length for sequences longer than chunk_size for EVM and Augustus",
    )
    prediction_ap.add_argument(
        "--transcriptome_models",
        type=FileType("r"),
        nargs="*",
        help="Models derived from transcriptomic data",
    )
    prediction_ap.add_argument(
        "--homology_models",
        type=FileType("r"),
        nargs="*",
        help="Models derived from protein alignments",
    )
    prediction_ap.add_argument(
        "--introns", type=FileType("r"), help="Introns to be used as hints for Augustus"
    )
    prediction_ap.add_argument(
        "--firststrand_expression",
        type=FileType("r"),
        help="Sorted by position first-strand RNASeq alignments used for coverage hints",
    )
    prediction_ap.add_argument(
        "--secondstrand_expression",
        type=FileType("r"),
        help="Sorted by position second-strand RNAseq alignments used for coverage hints",
    )
    prediction_ap.add_argument(
        "--unstranded_expression",
        type=FileType("r"),
        help="Sorted by position unstranded RNAseq alignments used for coverage hints",
    )
    prediction_ap.add_argument(
        "--repeats", type=FileType("r"), help="Repeat annotation GFF file."
    )
    prediction_ap.add_argument(
        "--homology_proteins",
        type=FileType("r"),
        required=True,
        help="Protein DB of sequences used for determining whether the evidence provided is "
        "full-length or not",
    )
    prediction_ap.add_argument(
        "--optimise_augustus",
        action="store_true",
        help="Enable augustus metaparameter optimisation",
    )
    prediction_ap.add_argument(
        "--kfold",
        type=int,
        default=8,
        help="Number of batches for augustus optimisation",
    )
    prediction_ap.add_argument(
        "--force_train",
        action="store_true",
        help="Re-train augustus even if the species is found in the 'augustus_config_path'",
    )
    prediction_ap.add_argument(
        "--augustus_runs",
        type=FileType("r"),
        nargs="*",
        help="File composed of 13 lines with SOURCE PRIORITY pairs for each of the types of "
        "evidence that can be used in\nan Augustus run. These evidence types are: "
        "gold models, silver models, bronze models, all models,\n"
        "gold introns, silver introns, protein models, coverage hints, repeat hints, "
        "high quality assemblies,\nlow quality assemblies, high quality proteins, and low "
        "quality proteins.",
    )
    prediction_ap.add_argument(
        "--EVM_weights",
        type=FileType("r"),
        required=True,
        help="Evidence modeler requires a weighting to be provided for each source of evidence,"
        " this file is the means to do so.",
    )
    prediction_ap.add_argument(
        "--hq_protein_alignments",
        type=FileType("r"),
        nargs="*",
        help="High confidence protein alignments to be used as hints for Augustus runs",
    )
    prediction_ap.add_argument(
        "--lq_protein_alignments",
        type=FileType("r"),
        nargs="*",
        help="Low confidence protein alignments to be used as hints for Augustus runs",
    )
    prediction_ap.add_argument(
        "--hq_assembly",
        type=FileType("r"),
        nargs="*",
        help="High confidence assemblies (for example from HiFi source) to be used as hints for "
        "Augustus runs",
    )
    prediction_ap.add_argument(
        "--lq_assembly",
        type=FileType("r"),
        nargs="*",
        help="Low confidence assemblies (short reads or low quality long reads) to be used as "
        "hints for Augustus runs",
    )
    prediction_ap.add_argument(
        "--mikado_utr_files",
        choices=UTR_SELECTION_OPTIONS,
        nargs="*",
        default=["augustus", "gold", "silver"],
        help=f"Choose any combination of space separated values from: "
        f"{' '.join(UTR_SELECTION_OPTIONS)}",
    )
    prediction_ap.add_argument(
        "--mikado_config",
        type=FileType("r"),
        help="Base configuration for Mikado consolidation stage.",
    )
    prediction_ap.add_argument(
        "--mikado_scoring",
        type=FileType("r"),
        help="Scoring file for Mikado pick at consolidation stage.",
    )
    prediction_ap.add_argument(
        "--do_glimmer",
        nargs="?",
        const=True,
        help="Enables GlimmerHmm predictions, optionally accepts a training directory",
    )
    prediction_ap.add_argument(
        "--do_snap",
        nargs="?",
        const=True,
        help="Enables SNAP predictions, optionally accepts a training directory",
    )
    prediction_ap.add_argument(
        "--do_codingquarry",
        nargs="?",
        const=True,
        help="Enables CodingQuarry predictions, optionally accepts a training directory",
    )
    prediction_ap.add_argument("--no_augustus", action="store_false")
    prediction_ap.add_argument(
        "--filter_top_n",
        type=int,
        default=0,
        help="Only output the top N transcripts that pass the self blast filter (0 outputs all)",
    )
    prediction_ap.add_argument(
        "--filter_max_identity",
        type=int,
        default=80,
        help="Maximum identity between models for redundancy classification",
    )
    prediction_ap.add_argument(
        "--filter_max_coverage",
        type=int,
        default=80,
        help="Maximum coverage between models for redundancy classification",
    )
    prediction_ap.add_argument(
        "--codingquarry_extra_params",
        help="Extra parameters for CodingQuarry predictions",
    )
    prediction_ap.add_argument(
        "--glimmer_extra_params", help="Extra parameters for glimmer predictions"
    )
    prediction_ap.add_argument(
        "--snap_extra_params", help="Extra parameters for snap predictions"
    )
    prediction_ap.add_argument(
        "--augustus_extra_params", help="Extra parameters for all Augustus predictions"
    )
    prediction_ap.add_argument(
        "--evm_extra_params",
        help="Extra parameters for EVM gene predictions consolidation",
    )
    prediction_ap.add_argument(
        "--min_train_models",
        type=int,
        default=400,
        help="Minimum number of training models",
    )
    prediction_ap.add_argument(
        "--max_train_models",
        type=int,
        default=1000,
        help="Maximum number of training models",
    )
    prediction_ap.add_argument(
        "--max_test_models", type=int, default=200, help="Maximum number of test models"
    )
    prediction_ap.add_argument(
        "--target_mono_exonic_percentage",
        type=int,
        default=20,
        help="Target percentage of mono-exonic models in the training set",
    )
    prediction_ap.add_argument(
        "--force_train_few_models",
        action="store_true",
        help="Train Augustus regardless monoexonic model ratio and number of models",
    )
    add_classification_parser_parameters(prediction_ap)

    args = reat_ap.parse_args()

    genetic_code = 0
    mikado_genetic_code = 0
    if args.reat_module == "transcriptome":
        genetic_code, mikado_genetic_code = transcriptome_cli_validation(args, reat_ap)
    return args, genetic_code, mikado_genetic_code


def main():
    print("Welcome to REAT")
    print("version -", VERSION)
    print("\nCommand-line call:")
    print(" ".join(sys.argv))
    print("\n")

    start_time = time.time()
    cli_arguments, genetic_code, mikado_genetic_code = parse_arguments()
    check_environment()

    # Check options.json
    try:
        json.load(cli_arguments.workflow_options_file)
    except json.JSONDecodeError as err:
        start, stop = max(0, err.pos - 20), err.pos + 20
        snippet = err.doc[start:stop]
        print(err)
        print(
            "... " if start else "",
            snippet,
            " ..." if stop < len(err.doc) else "",
            sep="",
        )

    if cli_arguments.reat_module == "transcriptome":
        rc = transcriptome_module(cli_arguments, genetic_code, mikado_genetic_code)
    elif cli_arguments.reat_module == "homology":
        rc = homology_module(cli_arguments)
    elif cli_arguments.reat_module == "prediction":
        rc = prediction_module(cli_arguments)
    else:
        rc = 1

    print(f"Done in {str(datetime.timedelta(seconds=time.time() - start_time))}")
    return rc


if __name__ == "__main__":
    main()
