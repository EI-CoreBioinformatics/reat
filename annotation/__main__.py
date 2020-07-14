#!/usr/bin/env python3

# __main__.py can result in the following outcomes
#     Success:
#         Reat is started in a 'run' backend (local, hpc)
#         Reat is started in a 'server' hosted on the network
#     Failure:
#         Reat fails to start or be submitted and the reason is provided to the user

# __main__.py parses all the arguments required, generates an input file for cromwell and submits a job to the requested
# backend. The arguments are validated using the json input schema defined in the validation directory.

# __main__.py supports multiple input json files which are subsequently merged into a single file, allows reusability
# of some parts of the input such as resource requirements and default extra parameters

# __main__.py accepts many customisation parameters in the command-line which are translated into the inputs.json before
# inputs.json is validated

import argparse
import json
from jsonschema import ValidationError, validators, Draft7Validator
import subprocess

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


def collect_arguments():
    ap = argparse.ArgumentParser(add_help=True)

    # General inputs
    ap.add_argument("--samples", nargs='+', type=argparse.FileType('r'),
                    help="Reads organised in the input specification for REAT, for more information please look at "
                         "the template file", required=True)
    ap.add_argument("--reference", type=argparse.FileType('r'),
                    help="Reference FASTA to annotate", required=True)
    ap.add_argument("--annotation", type=argparse.FileType('r'),
                    help="Annotation of the reference, this file will be used as the base for the new annotation "
                         "which will incorporate from the available evidence new gene models or update existing "
                         "ones")
    ap.add_argument("--computational_resources", type=argparse.FileType('r'),
                    help="Computational resources for REAT, please look at the template for more information",
                    required=True)
    ap.add_argument("--parameters_file", type=argparse.FileType('r'),
                    help="Base parameters file, this file can be the output of a previous REAT run which will be used "
                         "as the base for a new parameters file written to the output_parameters_file argument")
    ap.add_argument("--output_parameters_file", type=str,
                    help="REAT parameters file, this file will be used as the input for REAT. It provides the "
                         "arguments for the workflow.", default="reat_input.json")
    ap.add_argument("--workflow_options_file", type=argparse.FileType('r'),
                    help="Workflow execution options, includes cache usage and result directories structure and "
                         "location")

    # Mikado arguments
    mikado_parameters = ap.add_argument_group("mikado", "Parameters for Mikado runs")
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

    runtime = ap.add_mutually_exclusive_group(required=True)
    runtime.add_argument("--server", type=str,
                         help="Run the workflow in a cromwell server. Address and port of the cromwell server to "
                              "submit the workflow")
    runtime.add_argument("--run", type=argparse.FileType('r'),
                         help="Configuration file for the backend, please follow "
                              "https://cromwell.readthedocs.io/en/stable/backends/HPC/ for more information")

    # Aligner choices
    alignment_parameters = ap.add_argument_group("alignment", "Parameters for alignment of short and long reads")
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
    assembly_parameters = ap.add_argument_group("assembly", "Parameters for assembly of short and long reads")
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
    portcullis_parameters = ap.add_argument_group("portcullis", "Parameters specific to portcullis")
    portcullis_parameters.add_argument("--extra_parameters", type=str, help="Extra parameters for portcullis execution")

    # Orf calling
    orf_calling_parameters = ap.add_argument_group("ORF Caller", "Parameters for ORF calling programs")
    orf_calling_parameters.add_argument("--orf_caller", choices=['prodigal', 'transdecoder', 'none'],
                                        help="Choice of available orf calling softwares", default='none')
    orf_calling_parameters.add_argument("--orf_calling_proteins", type=argparse.FileType('r'),
                                        help="Set of proteins to be aligned to the genome for orf prediction by "
                                             "Transdecoder")

    args = ap.parse_args()
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

    cromwell_inputs["ei_annotation.wf_align.PR_stringtie_extra_parameters"] = cli_arguments.PR_stringtie_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.PR_scallop_extra_parameters"] = cli_arguments.PR_scallop_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.HQ_assembler"] = cli_arguments.HQ_assembler
    cromwell_inputs["ei_annotation.wf_align.LQ_assembler"] = cli_arguments.LQ_assembler

    cromwell_inputs["ei_annotation.wf_align.HQ_assembler_extra_parameters"] = cli_arguments.HQ_assembler_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.LQ_assembler_extra_parameters"] = cli_arguments.LQ_assembler_extra_parameters

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

    if cli_arguments.orf_caller is not "none":
        cromwell_inputs["ei_annotation.wf_main_mikado.orf_calling_program"] = cli_arguments.orf_caller

    if cli_arguments.separate_mikado_LQ:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = "true"
    else:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = "false"

    return cromwell_inputs


def main():
    cli_arguments = collect_arguments()
    # Print input file for cromwell
    cromwell_inputs = combine_arguments(cli_arguments)
    # Validate input against schema
    validate_cromwell_inputs(cromwell_inputs)

    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)

    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.transcriptome_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        return execute_cromwell(cli_arguments, cli_arguments.output_parameters_file,
                                workflow_options_file, wdl_file)


def validate_cromwell_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "reat.schema") as schema_file:
        with open(schema_file, 'r') as schema:
            reat_schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    all_validators["is_name"] = is_valid_name
    reat_validator = validators.create(meta_schema=reat_schema, validators=all_validators)
    reat_validator(reat_schema).validate(cromwell_inputs)


def execute_cromwell(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    if cli_arguments.server is None:
        return cromwell_run(input_parameters_filepath, workflow_options_file, wdl_file)
    else:
        return cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file)


def cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    # FIXME
    # To submit a workflow to a server subworkflows need to be included as a zip file, this could be a file that is pre-
    # packaged during installation and pointed to, or generated for each run. Pre-packaging seems to be a more sensible
    # option

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
    if sp_cromwell.returncode is not 0:
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


if __name__ == '__main__':
    main()
