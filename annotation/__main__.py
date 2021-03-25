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

import argparse
# Generates all required inputs and parses mikado's multiple runs extra configurations and places them in a reusable
# location where they can be edited for the user's convenience
import datetime
import io
import json
import os
import signal
import subprocess
import sys
import time
from textwrap import wrap

from annotation import VERSION
from annotation.homology import combine_arguments_homology, validate_homology_inputs
from annotation.transcriptome import combine_arguments, validate_transcriptome_inputs

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def check_environment(force_quit=True):
    """
    Check user's environment for required packages and their versions. Writing an to stderr if any of the software
    dependencies version doesn't match with the expected ones.
    :param force_quit: If any of the software versions doesn't match REAT will exit
    :return:
    """
    software_available = {
        "spaln": {"command": ["spaln"],
                  "result": "SPALN version 2.4.0"},
        "sortgrcd": {"command": ["sortgrcd", "--help"],
                     "result": "sortgrcd version 2.2"},
        "mikado": {"command": "mikado --version".split(' '), "result": "Mikado v2.0rc2"},
        "diamond": {"command": "diamond version".split(' '), "result": "diamond version 0.9.31"},
        "blastn": {"command": "blastn -version".split(' '), "result": "blastn: 2.7.1+"},
        "blastx": {"command": "blastx -version".split(' '), "result": "blastx: 2.7.1+"},
        "samtools": {"command": "samtools --version".split(' '), "result": "samtools 1.9"},
        "bamtools": {"command": "bamtools --version".split(' '), "result": "bamtools 2.5.1"},
        "gffread": {"command": "gffread --version".split(' '), "result": "0.12.2"},
        "gmap": {"command": "gmap --version".split(' '), "result": "version 2019-02-15"},
        "minimap2": {"command": "minimap2 --version".split(' '), "result": "2.17-r941"},
        "genometools": {"command": "gt --version".split(' '), "result": "gt (GenomeTools) 1.5.10"},
        "hisat2": {"command": "hisat2 --version".split(' '), "result": "version 2.1.0"},
        "star": {"command": "STAR --version".split(' '), "result": "2.7.3a"},
        "seqtk": {"command": ["seqtk"], "result": "Version: 1.3-r116-dirty"},
        "stringtie": {"command": "stringtie --version".split(' '), "result": "2.1.1"},
        "scallop": {"command": "scallop --version".split(' '), "result": "v0.10.4"},
        "scallop-lr": {"command": "scallop-lr --version".split(' '), "result": "v0.9.2"},
        "prodigal": {"command": "prodigal -v".split(' '), "result": "Prodigal V2.6.3: February, 2016"},
        "transdecoder": {"command": "TransDecoder.LongOrfs --version".split(' '),
                         "result": "TransDecoder.LongOrfs 5.5.0"},
        "portcullis": {"command": "portcullis --version".split(' '), "result": "portcullis 1.2.0"},
        "junctools": {"command": "junctools --version".split(' '), "result": "1.2.0"},
    }

    programs_not_found = set()
    for key, item in software_available.items():
        try:
            result = subprocess.run(item["command"], capture_output=True)
        except FileNotFoundError:
            programs_not_found.add(key)
            continue
        except RuntimeError as e:
            print(f"When executing {key}, found the following error:\n{e}\n"
                  f"Please ensure your environment and hardware support REAT's requirements", file=sys.stderr)
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
            print("\"", key, "\"", " version information:", sep="", file=sys.stderr)
            print('"""', file=sys.stderr)
            print(output.strip(), file=sys.stderr)
            print('"""', file=sys.stderr)
            print(file=sys.stderr)
            print("Does not contain indication of the required version:", file=sys.stderr)
            print('"""', file=sys.stderr)
            print(item["result"], file=sys.stderr)
            print('"""', file=sys.stderr)
            print(file=sys.stderr)
            print(file=sys.stderr)
        # Command not in path, wrong version or failed to execute
        if item["rc"] != 0:
            raise FileNotFoundError("Command {0} is missing, please check your PATH".format(key))
    if len(programs_not_found) > 0:
        print(f"When checking the environment, the software {', '.join([p for p in programs_not_found])}, not found.\n"
              f"Please make sure it is in the PATH environment variable of the shell executing REAT.\n\n"
              f"Currently PATH contains the following:")
        print('\n'.join(wrap(', '.join(os.environ['PATH'].split(os.pathsep)))))
        if force_quit:
            sys.exit(2)
    return software_available


def parse_arguments():
    """
    Parses and validates REAT's CLI arguments, the values inside defined on the CLI files have not yet been validated.

    :return: Object containing the validated CLI input arguments.
    """
    reat_ap = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reat_ap.add_argument("-j", "--jar_cromwell", type=argparse.FileType('r'),
                         help="Cromwell server jar file", required=True)
    reat_ap.add_argument("-r", "--runtime_configuration", type=argparse.FileType('r'),
                         help="Configuration file for the backend, please follow "
                              "https://cromwell.readthedocs.io/en/stable/backends/HPC/ for more information.\n"
                              "An example of this file can be found at ",
                         required=True)

    reat_ap.add_argument("-c", "--computational_resources", type=argparse.FileType('r'),
                         help="Computational resources for REAT, please look at the template for more information",
                         required=True)
    reat_ap.add_argument("-o", "--output_parameters_file", type=str,
                         help="REAT parameters file, this file will be used as the input for REAT. "
                              "It provides the arguments for the workflow runtime.",
                         default="reat_input.json")
    reat_ap.add_argument("--workflow_options_file", type=argparse.FileType('r'),
                         help="Workflow execution options, includes cache usage and result directories "
                              "structure and location")

    subparsers = reat_ap.add_subparsers(help="sub-command help", dest="reat_module")

    transcriptome_ap = subparsers.add_parser('transcriptome',
                                             help="Transcriptome module",
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # General inputs
    transcriptome_ap.add_argument("--reference", type=argparse.FileType('r'),
                                  help="Reference FASTA to annotate", required=True)
    transcriptome_ap.add_argument("--samples", nargs='+', type=argparse.FileType('r'),
                                  help="Reads organised in the input specification for REAT, for more information "
                                       "please look at https://github.com/ei-corebioinformatics/reat for an example")
    transcriptome_ap.add_argument("--csv_paired_samples", type=argparse.FileType('r'),
                                  help="CSV formatted input paired read samples. Without headers.\n"
                                       "The CSV fields are as follows name, strand, files (because"
                                       " this is an array that can contain one or more pairs, this fields' values are "
                                       "separated by semi-colon and space. Files in a pair are separated by semi-colon"
                                       "pairs are separated by a single space), merge, score, is_ref, "
                                       "exclude_redundant\n\n"
                                       "sample_strand takes values \'fr-firststrand\', \'fr-unstranded\', "
                                       "\'fr-secondstrand\'\n"
                                       "merge, is_ref and exclude_redundant are boolean and take values 'true', 'false'\n\n"
                                       "Example:\n"
                                       "PR1,fr-secondstrand,A_R1.fq;A_R2.fq /samples/paired/B1.fq;/samples/paired/B2.fq"
                                       ",false,2")
    transcriptome_ap.add_argument("--csv_long_samples", type=argparse.FileType('r'),
                                  help="CSV formatted input long read samples. Without headers.\n"
                                       "The CSV fields are as follows name, strand, files (space "
                                       "separated if there is more than one), quality, score, is_ref, "
                                       "exclude_redundant\n\n"
                                       "sample_strand takes values \'fr-firststrand\', \'fr-unstranded\', "
                                       "\'fr-secondstrand\'\n"
                                       "quality takes values 'low', 'high'\n"
                                       "is_ref and exclude_redundant are booleans and take values 'true', 'false'\n\n"
                                       "Example:\n"
                                       "Sample1,fr-firststrand,A.fq /samples/long/B.fq ./inputs/C.fq,low,2")
    transcriptome_ap.add_argument("--annotation", type=argparse.FileType('r'),
                                  help="Annotation of the reference, this file will be used as the base for the new"
                                       " annotation which will incorporate from the available evidence new gene models"
                                       " or update existing ones")
    transcriptome_ap.add_argument("--annotation_score", type=int, default=1,
                                  help="Score for models in the reference annotation file")
    transcriptome_ap.add_argument("--check_reference", action="store_true", default=False,
                                  help="At mikado stage, annotation models will be evaluated in the same manner as "
                                       "RNA-seq based models, removing any models deemed incorrect")
    transcriptome_ap.add_argument("--mode", choices=['basic', 'update', 'only_update'], default='basic',
                                  help="basic: Annotation models are treated the same as the RNA-Seq models at the pick"
                                       " stage."
                                       "update: Annotation models are prioritised but also novel loci are reported."
                                       "only_update: Annotation models are prioritised and non-reference loci are "
                                       "excluded.")
    transcriptome_ap.add_argument("--extra_junctions", type=argparse.FileType('r'),
                                  help="Extra junctions provided by the user, this file will be used as a set of valid"
                                       " junctions for alignment of short and long read samples, in the case of long"
                                       " reads, these junctions are combined with the results of portcullis whenever"
                                       " short read samples have been provided as part of the input datasets")
    transcriptome_ap.add_argument("--skip_mikado_long", action='store_true', default=False,
                                  help="Disables generation of the long read only mikado run")
    transcriptome_ap.add_argument("--parameters_file", type=argparse.FileType('r'),
                                  help="Base parameters file, this file can be the output of a previous REAT run "
                                       "which will be used as the base for a new parameters file written to the"
                                       " output_parameters_file argument")

    # Mikado arguments
    mikado_parameters = transcriptome_ap.add_argument_group("Mikado", "Parameters for Mikado runs")
    mikado_parameters.add_argument("--all_extra_config", type=argparse.FileType('r'),
                                   help="External configuration file for Paired and Long reads mikado")
    mikado_parameters.add_argument("--long_extra_config", type=argparse.FileType('r'),
                                   help="External configuration file for Long reads mikado run")
    mikado_parameters.add_argument("--lq_extra_config", type=argparse.FileType('r'),
                                   help="External configuration file for Low-quality long reads only mikado run "
                                        "(this is only applied when 'separate_mikado_LQ' is enabled)")
    mikado_parameters.add_argument("--all_scoring_file", type=argparse.FileType('r'),
                                   help="Mikado long and short scoring file", required=True)
    mikado_parameters.add_argument("--long_scoring_file", type=argparse.FileType('r'),
                                   help="Mikado long scoring file")
    mikado_parameters.add_argument("--long_lq_scoring_file", type=argparse.FileType('r'),
                                   help="Mikado low-quality long scoring file")
    mikado_parameters.add_argument("--homology_proteins", type=argparse.FileType('r'),
                                   help="Homology proteins database, used to score transcripts by Mikado")
    mikado_parameters.add_argument("--separate_mikado_LQ", type=bool,
                                   help="Specify whether or not to analyse low-quality long reads separately from "
                                        "high-quality, this option generates an extra set of mikado analyses "
                                        "including low-quality data")

    # Aligner choices
    alignment_parameters = transcriptome_ap.add_argument_group("Alignment",
                                                               "Parameters for alignment of short and long reads")
    alignment_parameters.add_argument("--short_reads_aligner", choices=['hisat', 'star'],
                                      help="Choice of short read aligner", default='hisat')
    alignment_parameters.add_argument("--HQ_aligner", choices=['minimap2', 'gmap'],
                                      help="Choice of aligner for high-quality long reads", default='gmap')
    alignment_parameters.add_argument("--LQ_aligner", choices=['minimap2', 'gmap'],
                                      help="Choice of aligner for low-quality long reads", default='minimap2')
    alignment_parameters.add_argument("--min_identity", type=float,
                                      help="Minimum alignment identity to retain transcript", default=0.9)
    alignment_parameters.add_argument("--min_intron_len", type=int,
                                      help="Where available, the minimum intron length allowed will be specified for "
                                           "the aligners", default=20)
    alignment_parameters.add_argument("--max_intron_len", type=int,
                                      help="Where available, the maximum intron length allowed will be specified for "
                                           "the aligners", default=200000)
    alignment_parameters.add_argument("--max_intron_len_ends", type=int,
                                      help="Where available, the maximum *boundary* intron length allowed will be "
                                           "specified for the aligner, when specified this implies max_intron_len "
                                           "only applies to the *internal* introns and this parameter to the *boundary*"
                                           " introns",
                                      default=100000)

    alignment_parameters.add_argument("--PR_hisat_extra_parameters", type=str,
                                      help="Extra command-line parameters for the selected short read aligner, please "
                                           "note that extra parameters are not validated and will have to match the "
                                           "parameters available for the selected read aligner")
    alignment_parameters.add_argument("--PR_star_extra_parameters", type=str,
                                      help="Extra command-line parameters for the selected short read aligner, please "
                                           "note that extra parameters are not validated and will have to match the "
                                           "parameters available for the selected read aligner")
    alignment_parameters.add_argument("--HQ_aligner_extra_parameters", type=str,
                                      help="Extra command-line parameters for the selected long read aligner, please "
                                           "note that extra parameters are not validated and will have to match the "
                                           "parameters available for the selected read aligner")
    alignment_parameters.add_argument("--LQ_aligner_extra_parameters", type=str,
                                      help="Extra command-line parameters for the selected long read aligner, please "
                                           "note that extra parameters are not validated and will have to match the "
                                           "parameters available for the selected read aligner")

    # Assembler choices
    assembly_parameters = transcriptome_ap.add_argument_group("Assembly",
                                                              "Parameters for assembly of short and long reads")
    assembly_parameters.add_argument("--skip_scallop", action='store_true', default=False)
    assembly_parameters.add_argument("--HQ_assembler",
                                     choices=["filter", "merge", "stringtie", "stringtie_collapse"],
                                     help="Choice of long read assembler."
                                          "\n- filter: Simply filters the reads based on identity and coverage"
                                          "- merge: cluster the input transcripts into loci, discarding "
                                          "\"duplicated\" transcripts (those with the same exact introns and fully "
                                          "contained or equal boundaries). This option also discards contained "
                                          "transcripts"
                                          "- stringtie: Assembles the long reads alignments into transcripts"
                                          "- stringtie_collapse: Cleans and collapses long reads but does not "
                                          "assemble them", default='filter')
    assembly_parameters.add_argument("--LQ_assembler",
                                     choices=["filter", "merge", "stringtie", "stringtie_collapse"],
                                     help="Choice of long read assembler."
                                          "\n- filter: Simply filters the reads based on identity and coverage"
                                          "- merge: cluster the input transcripts into loci, discarding "
                                          "\"duplicated\" transcripts (those with the same exact introns and fully "
                                          "contained or equal boundaries). This option also discards contained "
                                          "transcripts"
                                          "- stringtie: Assembles the long reads alignments into transcripts"
                                          "- stringtie_collapse: Cleans and collapses long reads but does not "
                                          "assembles them", default='stringtie_collapse')
    assembly_parameters.add_argument("--HQ_assembler_extra_parameters",
                                     help="Extra parameters for the long reads assembler, please note that extra "
                                          "parameters are not validated and will have to match the parameters "
                                          "available for the selected assembler")
    assembly_parameters.add_argument("--LQ_assembler_extra_parameters",
                                     help="Extra parameters for the long reads assembler, please note that extra "
                                          "parameters are not validated and will have to match the parameters "
                                          "available for the selected assembler")
    assembly_parameters.add_argument("--PR_stringtie_extra_parameters",
                                     help="Extra parameters for stringtie, please note that extra "
                                          "parameters are not validated and will have to match the parameters "
                                          "available for stringtie")
    assembly_parameters.add_argument("--PR_scallop_extra_parameters",
                                     help="Extra parameters for scallop, please note that extra "
                                          "parameters are not validated and will have to match the parameters "
                                          "available for scallop")

    # Portcullis extra parameters
    portcullis_parameters = transcriptome_ap.add_argument_group("Portcullis", "Parameters specific to portcullis")
    portcullis_parameters.add_argument("--extra_parameters", type=str, help="Extra parameters for portcullis execution")

    # Orf calling
    orf_calling_parameters = transcriptome_ap.add_argument_group("ORF Caller", "Parameters for ORF calling programs")
    orf_calling_parameters.add_argument("--orf_caller", choices=['prodigal', 'transdecoder', 'none'],
                                        help="Choice of available orf calling softwares", default='none')
    orf_calling_parameters.add_argument("--orf_calling_proteins", type=argparse.FileType('r'),
                                        help="Set of proteins to be aligned to the genome for orf prediction by "
                                             "Transdecoder")

    homology_ap = subparsers.add_parser('homology', help="Homology module",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    homology_ap.add_argument("--genome", type=argparse.FileType('r'),
                             help="Fasta file of the genome to annotate",
                             required=True)
    homology_ap.add_argument("--alignment_species", type=str,
                             help="Species specific parameters, select a value from the first or second column of "
                                  "https://raw.githubusercontent.com/ogotoh/spaln/master/table/gnm2tab",
                             required=True)
    homology_ap.add_argument("--annotations_csv", type=argparse.FileType('r'),
                             help="CSV file with reference annotations to extract proteins/cdnas for spliced alignments"
                                  " in csv format. The CSV fields are as follows genome_fasta,annotation_gff  "
                                  "e.g Athaliana.fa,Athaliana.gff",
                             required=True)
    homology_ap.add_argument("--annotation_filters",
                             choices=['all', 'none', 'exon_len', 'intron_len', 'internal_stop', 'aa_len',
                                      'splicing'],
                             nargs='+',
                             help="Filter annotation coding genes by the filter types specified",
                             default=['none'])
    homology_ap.add_argument("--mikado_config", type=argparse.FileType('r'),
                             help="Base configuration for Mikado consolidation stage.",
                             required=True)
    homology_ap.add_argument("--mikado_scoring", type=argparse.FileType('r'),
                             help="Scoring file for Mikado pick at consolidation stage.",
                             required=True)
    homology_ap.add_argument("--junctions", type=argparse.FileType('r'),
                             help="Validated junctions BED file for use in Mikado consolidation stage.")
    homology_ap.add_argument("--utrs", type=argparse.FileType('r'),
                             help="Gene models that may provide UTR extensions to the homology based models at the "
                                  "mikado stage")
    homology_ap.add_argument("--pick_extra_config", type=argparse.FileType('r'),
                             help="Extra configuration for Mikado pick stage")
    homology_ap.add_argument("--min_cdna_length", type=int, default=100,
                             help="Minimum cdna length for models to consider in Mikado consolidation stage")
    homology_ap.add_argument("--max_intron_length", type=int, default=1000000,
                             help="Maximum intron length for models to consider in Mikado consolidation stage")
    homology_ap.add_argument("--filter_min_cds", type=int,
                             help="If 'aa_len' filter is enabled for annotation coding features, any CDS smaller than"
                                  "this parameter will be filtered out",
                             default=20)
    homology_ap.add_argument("--filter_max_intron", type=int,
                             help="If 'intron_len' filter is enabled, any features "
                                  "with introns longer than this parameter will be filtered out",
                             default=200000)
    homology_ap.add_argument("--filter_min_exon", type=int,
                             help="If 'exon_len' filter is enabled, any features "
                                  "with exons shorter than this parameter will be filtered out",
                             default=20)
    homology_ap.add_argument("--alignment_min_exon_len", type=int, help="Minimum exon length, alignment parameter",
                             default=20)
    homology_ap.add_argument("--alignment_filters",
                             choices=['all', 'none', 'exon_len', 'intron_len', 'internal_stop', 'aa_len',
                                      'splicing'],
                             help="Filter alignment results by the filter types specified", nargs='+', default=['none'])
    homology_ap.add_argument("--alignment_min_identity", type=int, help="Minimum identity filter for alignments",
                             default=50)
    homology_ap.add_argument("--alignment_min_coverage", type=int, help="Minimum coverage filter for alignments",
                             default=80)
    homology_ap.add_argument("--alignment_max_per_query", type=int, default=4,
                             help="Maximum number of alignments per input query protein")
    homology_ap.add_argument("--alignment_recursion_level", type=int, default=6,
                             help="SPALN's Q value, indicating the level of recursion for the Hirschberg algorithm")
    homology_ap.add_argument("--alignment_show_intron_length", action='store_true',
                             help="Add an attribute to the alignment gff with the maximum intron len for each mRNA")
    homology_ap.add_argument("--exon_f1_filter", type=int,
                             help="Filter alignments scored against its original structure with a CDS exon f1 "
                                  "lower than this value")
    homology_ap.add_argument("--junction_f1_filter", type=int,
                             help="Filter alignments scored against its original structure with a CDS junction f1 "
                                  "lower than this value")

    args = reat_ap.parse_args()

    if args.reat_module == 'transcriptome':
        if args.separate_mikado_LQ:
            if not args.long_lq_scoring_file:
                reat_ap.error("When '--separate_mikado_LQ' is enabled, --long_lq_scoring_file is required, please "
                              "provide it.")

        if args.samples and (args.csv_paired_samples or args.csv_long_samples):
            reat_ap.error("Conflicting arguments '--samples' and ['--csv_paired_samples' or '--csv_long_samples'] "
                          "provided, please choose one of csv or json sample input format")
        if not args.samples and not args.csv_paired_samples and not args.csv_long_samples:
            reat_ap.error("Please provide at least one of --samples, --csv_paired_samples, --csv_long_samples")
    return args


def transcriptome_module(cli_arguments):
    """
    Collects the CLI arguments and combines them with CLI defined input files.
    The resulting object is validated and the final inputs are written to a json Cromwell input file.
    :param cli_arguments: Validated CLI values with un-inspected files
    :return: Return code from Cromwell
    """
    # Print input file for cromwell
    cromwell_inputs = combine_arguments(cli_arguments)
    # Validate input against schema
    validate_transcriptome_inputs(cromwell_inputs)
    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.transcriptome_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        return execute_cromwell(cli_arguments.runtime_configuration, cli_arguments.jar_cromwell,
                                cli_arguments.output_parameters_file, workflow_options_file, wdl_file)


def homology_module(cli_arguments):
    cromwell_inputs = combine_arguments_homology(cli_arguments)
    validate_homology_inputs(cromwell_inputs)
    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.homology_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        return execute_cromwell(cli_arguments.runtime_configuration, cli_arguments.jar_cromwell,
                                cli_arguments.output_parameters_file, workflow_options_file, wdl_file)


def execute_cromwell(workflow_configuration_file, jar_cromwell, input_parameters_filepath, workflow_options_file, wdl_file):
    return cromwell_run(workflow_configuration_file, jar_cromwell,
                        input_parameters_filepath, workflow_options_file, wdl_file)


def cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    # FIXME
    #  Package the pipeline dependencies from the installed resources folder into a zip file for submitting

    subprocess.run(
        ["cromwell", "submit", "-h", cli_arguments.server, "-i",
         input_parameters_filepath, "-o", workflow_options_file, wdl_file]
    )
    # FIXME return the code of the request or some mapping to useful error codes
    return 0


def kill_cromwell(sig, frame):
    raise KeyboardInterrupt


def cromwell_run(workflow_configuration_file, jar_cromwell, input_parameters_filepath, workflow_options_file, wdl_file,
                 log_level="INFO"):
    formatted_command_line = ["java", f"-Dconfig.file={str(workflow_configuration_file.name)}",
                              f"-DLOG_LEVEL={log_level}",
                              "-jar", str(jar_cromwell.name),
                              "run",
                              "-i", str(input_parameters_filepath)]
    if workflow_options_file:
        formatted_command_line.extend(["-o", str(workflow_options_file)])

    formatted_command_line.extend(["-m", "run_details.json", str(wdl_file)])

    print("Starting:")
    print(' '.join(formatted_command_line))
    cromwell_sp_output = io.StringIO()
    try:
        sp_cromwell = subprocess.Popen(
            formatted_command_line,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        signal.signal(signal.SIGINT, kill_cromwell)
        signal.signal(signal.SIGTERM, kill_cromwell)
        while True:
            output = sp_cromwell.stdout.readline()
            if sp_cromwell.poll() is not None:
                break
            if output:
                cromwell_sp_output.write(output)
                print(output.strip())
                sys.stdout.flush()
    except KeyboardInterrupt:
        sp_cromwell.send_signal(signal.SIGINT)
        sp_cromwell.wait()
        for line in sp_cromwell.stdout:
            print(line.strip())
    sys.stdout.flush()
    rc = sp_cromwell.poll()
    if rc != 0:
        if rc == 130:
            print("REAT stopped by user request")
        else:
            sentinel = "Check the content of stderr for potential additional information: "
            cromwell_sp_output_str = cromwell_sp_output.getvalue()
            error_file_start_pos = cromwell_sp_output_str.find(sentinel)
            for line in sp_cromwell.stderr:
                print(line)
            if error_file_start_pos < 0:
                # FIXME Unhandled error
                print("Unhandled errors, please report this as an issue to add support for improved messages and "
                      "suggestions for actions on how to resolve it")
                print(cromwell_sp_output_str)
            else:
                error_file = cromwell_sp_output_str[error_file_start_pos + len(sentinel):].split("\n")[0][:-1]
                print("\n\n\nREAT Failed, the following file might contain information with the reasons behind"
                      " the failure")
                print(error_file)
                if os.path.exists(error_file):
                    with open(error_file, 'r') as failed_job_stderr_file:
                        print(failed_job_stderr_file.read())
                else:
                    print("The stderr file for the process that failed was not found, this can happen when the job was "
                          "killed outside REAT i.e when there was an out-of-memory issue, please check the logs for "
                          "failed jobs and provide the necessary resources for these to complete")
    return rc


def main():
    print("Welcome to REAT")
    print("version -", VERSION)
    print("\nCommand-line call:")
    print(' '.join(sys.argv))
    print("\n")

    start_time = time.time()
    cli_arguments = parse_arguments()
    check_environment()

    if cli_arguments.reat_module == "transcriptome":
        rc = transcriptome_module(cli_arguments)
    elif cli_arguments.reat_module == "homology":
        rc = homology_module(cli_arguments)
    else:
        rc = 1

    print(f"Done in {str(datetime.timedelta(seconds=time.time() - start_time))}")
    return rc


if __name__ == '__main__':
    main()
