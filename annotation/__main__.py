#!/usr/bin/env python3

# __main__.py can result in the following outcomes
#     Success:
#         Reat is started in a 'run' backend (local, hpc)
#         Reat is started in a 'server' hosted on the network
#     Failure:
#         Reat fails to start or be submitted and the reason is provided to the user

# __main__.py checks the user environment to ensure the required software is available if this is not the case
# the user is informed of which software is and which isn't available

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
import json
import subprocess
import sys

from jsonschema import ValidationError, validators, Draft7Validator

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources


def is_valid_name(validator, value, instance, schema):
    if not isinstance(instance, str):
        yield ValidationError("%r is not a string" % instance)

    if value and set(instance).intersection("][/?\\\'\" .*$)(}{"):
        yield ValidationError("%r is not alphanumeric" % (instance,))


def check_environment():
    # Check the following software is installed and available in the user's environment
    software_available = {
        "spaln": {"command": ["spaln"],
                  "result": "SPALN version 2.4.03"},
        "mikado": {"command": "mikado --version".split(' '), "result": "Mikado v2.0rc2"},
        "diamond": {"command": "diamond version".split(' '), "result": "diamond version 0.9.31"},
        "blastn": {"command": "blastn -version".split(' '), "result": "blastn: 2.7.1+"},
        "blastx": {"command": "blastx -version".split(' '), "result": "blastx: 2.7.1+"},
        "samtools": {"command": "samtools --version".split(' '), "result": "samtools 1.9"},
        "bamtools": {"command": "bamtools --version".split(' '), "result": "bamtools 2.5.1"},
        "gffread": {"command": "gffread --version".split(' '), "result": "0.10.1"},
        "gmap": {"command": "gmap --version".split(' '), "result": "version 2019-02-15"},
        "minimap2": {"command": "minimap2 --version".split(' '), "result": "2.17-r941"},
        "genometools": {"command": "gt --version".split(' '), "result": "gt (GenomeTools) 1.5.10"},
        "hisat2": {"command": "hisat2 --version".split(' '), "result": "version 2.1.0"},
        "star": {"command": "STAR --version".split(' '), "result": "2.7.3a"},
        "seqtk": {"command": ["seqtk"], "result": "Version: 1.3-r114-dirty"},
        "stringtie": {"command": "stringtie --version".split(' '), "result": "2.1.1"},
        "scallop": {"command": "scallop --version".split(' '), "result": "v0.10.4"},
        "scallop-lr": {"command": "scallop-lr --version".split(' '), "result": "v0.9.2"},
        "prodigal": {"command": "prodigal -v".split(' '), "result": "Prodigal V2.6.3: February, 2016"},
        "transdecoder": {"command": "TransDecoder.LongOrfs --version".split(' '),
                         "result": "TransDecoder.LongOrfs 5.5.0"},
        "portcullis": {"command": "portcullis --version".split(' '), "result": "portcullis 1.2.0"},
        "BioPerl": {"command": ["perl", "-MBio::Root::Version", "-e", """print $Bio::Root::Version::VERSION\n"""],
                    "result": "1.7.7"}
    }

    for key, item in software_available.items():
        result = subprocess.run(item["command"], capture_output=True)
        output = result.stdout.decode()
        output += result.stderr.decode()
        item["rc"] = result.returncode
        if key == "seqtk":
            if "Version" in output:
                item["rc"] = 0
        if key == "spaln":
            if "No input seq file !" in output:
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

    return software_available


def parse_arguments():
    reat_ap = argparse.ArgumentParser(add_help=True)

    reat_ap.add_argument("--computational_resources", type=argparse.FileType('r'),
                         help="Computational resources for REAT, please look at the template for more information",
                         required=True)
    reat_ap.add_argument("--output_parameters_file", type=str,
                         help="REAT parameters file, this file will be used as the input for REAT. "
                              "It provides the arguments for the workflow runtime.",
                         default="reat_input.json")
    reat_ap.add_argument("--workflow_options_file", type=argparse.FileType('r'),
                         help="Workflow execution options, includes cache usage and result directories "
                              "structure and location")

    runtime = reat_ap.add_mutually_exclusive_group(required=True)

    runtime.add_argument("--server", type=str,
                         help="Run the workflow in a cromwell server. Address and port of the cromwell server to "
                              "submit the workflow")
    runtime.add_argument("--run", type=argparse.FileType('r'),
                         help="Configuration file for the backend, please follow "
                              "https://cromwell.readthedocs.io/en/stable/backends/HPC/ for more information")

    subparsers = reat_ap.add_subparsers(help="sub-command help", dest="reat_module")

    transcriptome_ap = subparsers.add_parser('transcriptome',
                                             help="Transcriptome module")
    # General inputs
    transcriptome_ap.add_argument("--samples", nargs='+', type=argparse.FileType('r'),
                                  help="Reads organised in the input specification for REAT, for more information "
                                       "please look at the template file", required=True)
    transcriptome_ap.add_argument("--reference", type=argparse.FileType('r'),
                                  help="Reference FASTA to annotate", required=True)
    transcriptome_ap.add_argument("--annotation", type=argparse.FileType('r'),
                                  help="Annotation of the reference, this file will be used as the base for the new"
                                       " annotation which will incorporate from the available evidence new gene models"
                                       " or update existing ones")
    transcriptome_ap.add_argument("--parameters_file", type=argparse.FileType('r'),
                                  help="Base parameters file, this file can be the output of a previous REAT run "
                                       "which will be used as the base for a new parameters file written to the"
                                       " output_parameters_file argument")

    # Mikado arguments
    mikado_parameters = transcriptome_ap.add_argument_group("mikado", "Parameters for Mikado runs")
    mikado_parameters.add_argument("--run_mikado_homology", type=bool,
                                   help="Use the homology proteins provided for scoring transcripts")
    mikado_parameters.add_argument("--mikado_scoring_file", type=argparse.FileType('r'),
                                   help="Mikado scoring file", required=True)
    mikado_parameters.add_argument("--homology_proteins", type=argparse.FileType('r'),
                                   help="Homology proteins database, used to score transcripts by Mikado")
    mikado_parameters.add_argument("--do_mikado_homology", type=bool,
                                   help="Specify whether or not to use homology data for Mikado scoring of gene models")
    mikado_parameters.add_argument("--separate_mikado_LQ", type=bool,
                                   help="Specify whether or not to analyse low-quality long reads separately from "
                                        "high-quality, this option generates an extra set of mikado analyses "
                                        "including low-quality data")

    # Aligner choices
    alignment_parameters = transcriptome_ap.add_argument_group("alignment",
                                                               "Parameters for alignment of short and long reads")
    alignment_parameters.add_argument("--short_reads_aligner", choices=['hisat', 'star'],
                                      help="Choice of short read aligner", default='hisat')
    alignment_parameters.add_argument("--HQ_aligner", choices=['minimap2', 'gmap'],
                                      help="Choice of aligner for high-quality long reads", default='gmap')
    alignment_parameters.add_argument("--LQ_aligner", choices=['minimap2', 'gmap'],
                                      help="Choice of aligner for low-quality long reads", default='gmap')
    alignment_parameters.add_argument("--min_identity", type=float,
                                      help="Minimum alignment identity to retain transcript", default=0.9)
    alignment_parameters.add_argument("--min_intron_len", type=int,
                                      help="Where available, the minimum intron length allowed will be specified for "
                                           "the aligners", default=200)
    alignment_parameters.add_argument("--max_intron_len", type=int,
                                      help="Where available, the maximum intron length allowed will be specified for "
                                           "the aligners", default=2000)
    alignment_parameters.add_argument("--max_intron_len_middle", type=int,
                                      help="Where available, the maximum *internal* intron length allowed will be "
                                           "specified for the aligner, when specified this implies max_intron_length "
                                           "only applies to the *ends* and this parameter to the *internal* introns",
                                      default=4000)

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
    assembly_parameters = transcriptome_ap.add_argument_group("assembly",
                                                              "Parameters for assembly of short and long reads")
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
                                          "assemble them", default='merge')
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
                                          "assembles them", default='stringtie')
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
    portcullis_parameters = transcriptome_ap.add_argument_group("portcullis", "Parameters specific to portcullis")
    portcullis_parameters.add_argument("--extra_parameters", type=str, help="Extra parameters for portcullis execution")

    # Orf calling
    orf_calling_parameters = transcriptome_ap.add_argument_group("ORF Caller", "Parameters for ORF calling programs")
    orf_calling_parameters.add_argument("--orf_caller", choices=['prodigal', 'transdecoder', 'none'],
                                        help="Choice of available orf calling softwares", default='none')
    orf_calling_parameters.add_argument("--orf_calling_proteins", type=argparse.FileType('r'),
                                        help="Set of proteins to be aligned to the genome for orf prediction by "
                                             "Transdecoder")

    homology_ap = subparsers.add_parser('homology', help="Homology module")

    homology_ap.add_argument("--genome", type=argparse.FileType('r'),
                             help="Fasta file of the genome to annotate")
    homology_ap.add_argument("--annotations", nargs='+', type=argparse.FileType('r'),
                             help="Reference annotations to extract proteins/cdnas for spliced alignments")
    homology_ap.add_argument("--annotation_filters",
                             choices=['all', 'none', 'intron_len', 'internal_stop', 'aa_len', 'splicing'], nargs='+',
                             help="Filter annotation coding genes by the filter types specified")
    homology_ap.add_argument("--annotation_min_cds", type=int,
                             help="If 'aa_len' filter is enabled for annotation coding features, any CDS smaller than"
                                  "this parameter will be filtered out")
    homology_ap.add_argument("--alignment_species", type=str,
                             help="Available aligner species, for more information, please look at URL")
    homology_ap.add_argument("--alignment_min_exon_len", type=int, help="Minimum exon length, alignment parameter")
    homology_ap.add_argument("--alignment_filters",
                             choices=['all', 'none', 'intron_len', 'internal_stop', 'aa_len', 'splicing'],
                             help="Filter alignment results by the filter types specified", nargs='+')
    homology_ap.add_argument("--alignment_min_identity", type=int, help="Minimum identity filter for alignments")
    homology_ap.add_argument("--alignment_min_coverage", type=int, help="Minimum coverage filter for alignments")

    args = reat_ap.parse_args()
    return args


def combine_arguments(cli_arguments):
    computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources
    for s in cli_arguments.samples:
        sample = json.load(s)
        cromwell_inputs.update(sample)

    cromwell_inputs["ei_annotation.reference_genome"] = cli_arguments.reference.name
    cromwell_inputs["ei_annotation.mikado_scoring_file"] = cli_arguments.mikado_scoring_file.name

    cromwell_inputs["ei_annotation.wf_align.PR_hisat_extra_parameters"] = cli_arguments.PR_hisat_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.PR_star_extra_parameters"] = cli_arguments.PR_star_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.HQ_aligner"] = cli_arguments.HQ_aligner
    cromwell_inputs["ei_annotation.wf_align.LQ_aligner"] = cli_arguments.LQ_aligner

    cromwell_inputs["ei_annotation.wf_align.HQ_aligner_extra_parameters"] = cli_arguments.HQ_aligner_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.LQ_aligner_extra_parameters"] = cli_arguments.LQ_aligner_extra_parameters

    cromwell_inputs[
        "ei_annotation.wf_align.PR_stringtie_extra_parameters"] = cli_arguments.PR_stringtie_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.PR_scallop_extra_parameters"] = cli_arguments.PR_scallop_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.HQ_assembler"] = cli_arguments.HQ_assembler
    cromwell_inputs["ei_annotation.wf_align.LQ_assembler"] = cli_arguments.LQ_assembler

    cromwell_inputs[
        "ei_annotation.wf_align.HQ_assembler_extra_parameters"] = cli_arguments.HQ_assembler_extra_parameters
    cromwell_inputs[
        "ei_annotation.wf_align.LQ_assembler_extra_parameters"] = cli_arguments.LQ_assembler_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.min_intron_len"] = cli_arguments.min_intron_len
    cromwell_inputs["ei_annotation.wf_align.max_intron_len"] = cli_arguments.max_intron_len
    cromwell_inputs["ei_annotation.wf_align.max_intron_len_middle"] = cli_arguments.max_intron_len_middle
    cromwell_inputs["ei_annotation.wf_align.min_identity"] = cli_arguments.min_identity

    if cli_arguments.run_mikado_homology:
        cromwell_inputs["ei_annotation.wf_main_mikado.run_mikado_homology"] = "true"
    else:
        cromwell_inputs["ei_annotation.wf_main_mikado.run_mikado_homology"] = "false"

    if cli_arguments.homology_proteins is not None:
        cromwell_inputs["ei_annotation.homology_proteins"] = cli_arguments.homology_proteins

    if cli_arguments.orf_calling_proteins is not None:
        cromwell_inputs["ei_annotation.orf_calling_proteins"] = cli_arguments.orf_calling_proteins

    if cli_arguments.orf_caller != "none":
        cromwell_inputs["ei_annotation.wf_main_mikado.orf_calling_program"] = cli_arguments.orf_caller

    if cli_arguments.separate_mikado_LQ:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = "true"
    else:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = "false"

    return cromwell_inputs


def validate_transcriptome_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "transcriptome.schema") as schema_file:
        with open(schema_file, 'r') as schema:
            reat_schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    all_validators["is_name"] = is_valid_name
    reat_validator = validators.create(meta_schema=reat_schema, validators=all_validators)
    reat_validator(reat_schema).validate(cromwell_inputs)


def execute_cromwell(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    if cli_arguments.run is None:
        return cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file)
    else:
        return cromwell_run(input_parameters_filepath, workflow_options_file, wdl_file)


def cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    # FIXME
    #  Package the pipeline dependencies from the installed resources folder into a zip file for submitting

    subprocess.run(
        ["cromwell", "submit", "-h", cli_arguments.server, "-i",
         input_parameters_filepath, "-o", workflow_options_file, wdl_file]
    )
    # FIXME return the code of the request or some mapping to useful error codes
    return 0


def cromwell_run(input_parameters_filepath, workflow_options_file, wdl_file):
    workflow_options = "-o " + workflow_options_file if workflow_options_file is not None else ""
    print("Starting:")
    print(' '.join(["cromwell", "run", "-i", input_parameters_filepath, workflow_options, str(wdl_file)]))
    sp_cromwell = subprocess.run(
        ["cromwell", "run", "-i", input_parameters_filepath, wdl_file],
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    if sp_cromwell.returncode != 0:
        sentinel = "Check the content of stderr for potential additional information: "
        error_file_start_pos = sp_cromwell.stdout.find(sentinel)
        if error_file_start_pos < 0:
            # FIXME Unhandled error
            print("Unhandled errors, please create an issue in the github to add support for improved messages and "
                  "actions on how to resolve it")
            print(sp_cromwell.stdout)
        else:
            error_file = sp_cromwell.stdout[error_file_start_pos + len(sentinel):].split("\n")[0][:-1]
            print(error_file)
            with open(error_file, 'r') as failed_job_stderr_file:
                print(failed_job_stderr_file.read())
    return sp_cromwell.returncode


def main():
    try:
        check_environment()
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        sys.exit(2)
    cli_arguments = parse_arguments()
    if cli_arguments.reat_module == "transcriptome":
        return transcriptome_module(cli_arguments)
    elif cli_arguments.reat_module == "homology":
        return homology_module(cli_arguments)


def transcriptome_module(cli_arguments):
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
        return execute_cromwell(cli_arguments, cli_arguments.output_parameters_file,
                                workflow_options_file, wdl_file)


def combine_arguments_homology(cli_arguments):
    computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources
    for s in cli_arguments.annotations:
        annotation = json.load(s)
        cromwell_inputs.update(annotation)

    cromwell_inputs["ei_homology.genome_to_annotate"] = cli_arguments.genome.name
    cromwell_inputs["ei_homology.AlignProteins.species"] = cli_arguments.alignment_species

    # Optional extra parameters
    if cli_arguments.annotation_filters:
        cromwell_inputs["ei_homology.PrepareAnnotations.filters"] = cli_arguments.annotation_filters
    if cli_arguments.annotation_min_cds:
        cromwell_inputs["ei_homology.PrepareAnnotations.min_cds_len"] = cli_arguments.annotation_min_cds
    if cli_arguments.alignment_min_exon_len:
        cromwell_inputs["ei_homology.AlignProteins.min_exon_len"] = cli_arguments.alignment_min_exon_len
    if cli_arguments.alignment_filters:
        cromwell_inputs["ei_homology.AlignProteins.filters"] = cli_arguments.alignment_filters
    if cli_arguments.alignment_min_identity:
        cromwell_inputs["ei_homology.AlignProteins.min_identity"] = cli_arguments.alignment_min_identity
    if cli_arguments.alignment_min_coverage:
        cromwell_inputs["ei_homology.AlignProteins.min_coverage"] = cli_arguments.alignment_min_coverage

    return cromwell_inputs


def validate_homology_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "homology.schema") as schema_file:
        with open(schema_file, 'r') as schema:
            schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    reat_validator = validators.create(meta_schema=schema, validators=all_validators)
    reat_validator(schema).validate(cromwell_inputs)


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
        return execute_cromwell(cli_arguments, cli_arguments.output_parameters_file,
                                workflow_options_file, wdl_file)


if __name__ == '__main__':
    main()
