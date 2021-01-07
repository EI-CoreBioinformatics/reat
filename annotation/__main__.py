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
import io
import argparse
import json
import os
import signal
import subprocess
import sys
import time
from collections import defaultdict
from json.decoder import JSONDecodeError
from textwrap import wrap

from jsonschema import ValidationError, validators, Draft7Validator

from annotation import VERSION
from pathlib import Path

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


def eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def is_valid_name(validator, value, instance, schema):
    if not isinstance(instance, str):
        yield ValidationError("%r is not a string" % instance)

    if value and set(instance).intersection("][/?\\\'\" .*$)(}{"):
        yield ValidationError("%r is not alphanumeric" % (instance,))


def separate_mikado_config(mikado_config, mikado_run):
    # Check empty/None creates no files!
    config = load(open(mikado_config, 'r'), Loader=Loader)
    prepare = config.get('prepare', None)
    serialise = config.get('serialise', None)
    pick = config.get('pick', None)

    prepare_path = Path(mikado_run).joinpath("prepare.yaml") if prepare else None
    serialise_path = Path(mikado_run).joinpath("serialise.yaml") if serialise else None
    pick_path = Path(mikado_run).joinpath("pick.yaml") if pick else None

    if any((prepare_path, serialise_path, pick_path)):
        Path(mikado_run).mkdir(exist_ok=True)
        if prepare:
            print(dump({'prepare': prepare}, default_flow_style=False), file=open(prepare_path, 'w'))
        if serialise:
            print(dump({'serialise': serialise}, default_flow_style=False), file=open(serialise_path, 'w'))
        if pick:
            print(dump({'pick': pick}, default_flow_style=False), file=open(pick_path, 'w'))
    return prepare_path, serialise_path, pick_path


def check_environment():
    # Check the following software is installed and available in the user's environment
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
        # "BioPerl": {"command": ["perl", "-MBio::Root::Version", "-e", """print $Bio::Root::Version::VERSION\n"""],
        #             "result": "1.7.7"}
    }

    programs_not_found = set()
    for key, item in software_available.items():
        try:
            result = subprocess.run(item["command"], capture_output=True)
        except FileNotFoundError:
            programs_not_found.add(key)
            continue
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
        sys.exit(2)
    return software_available


def parse_arguments():
    reat_ap = argparse.ArgumentParser(add_help=True, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # runtime = reat_ap.add_mutually_exclusive_group(required=True)
    #
    # runtime.add_argument("--server", type=str,
    #                      help="Run the workflow in a cromwell server. Address and port of the cromwell server to "
    #                           "submit the workflow")
    reat_ap.add_argument("-r", "--runtime_configuration", type=argparse.FileType('r'),
                         help="Configuration file for the backend, please follow "
                              "https://cromwell.readthedocs.io/en/stable/backends/HPC/ for more information.\n"
                              "An example of this file can be found at ",
                         required=True)

    reat_ap.add_argument("-c", "--computational_resources", type=argparse.FileType('r'),
                         help="Computational resources for REAT, please look at the template for more information",
                         )  #required=True)
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
                                   help="Mikado long scoring file", required=True)
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
                             help="Available aligner species, for more information, please look at URL",
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


def validate_long_samples(samples):
    result = {}
    names = []
    errors = defaultdict(list)
    strands = ("fr-firststrand", "fr-secondstrand", "fr-unstranded")
    quality_choices = ('high', 'low')
    hq_samples = []
    lq_samples = []

    for line in samples:
        out_files = []
        fields = line.rstrip().split(",")
        try:
            name, strand, files, quality = fields[:4]
        except ValueError as e:
            errors[line].append(f"Unexpected input '{fields}'\n\t\tPlease make sure this is a csv file with at minimum"
                                f"the following fields name, quality, strand, files")
        if name in names:
            errors[line].append(
                ("Non-unique name '{}' specified, please make sure sample names are unique".format(name)))

        if quality.lower() not in quality_choices:
            errors[line].append(
                "Incorrect quality '{}' specification, please choose one of {}".format(quality, quality_choices))

        for file in files.split(' '):
            if not os.path.exists(file):
                errors[line].append(("File not found: '{}'".format(file)))
            out_files.append(file)
        if strand.lower() not in strands:
            errors[line].append(
                ("Incorrect strand '{}' specification, please choose one of {}".format(strand, strands)))
        names.append(name)

        try:
            score = float(fields[4])
        except ValueError as e:
            errors[line].append("Cannot parse score value '{}' as float".format(fields[4]))
        except IndexError as e:
            score = 0

        try:
            is_ref = fields[5].lower()
            if is_ref in ("true", "false"):
                is_ref = True if is_ref == 'true' else False
            else:
                errors[line].append("is_ref field with value '{}' should be either 'true' or 'false'".format(is_ref))
        except IndexError:
            is_ref = False

        try:
            exclude_redundant = fields[6].lower()
            if exclude_redundant in ("true", "false"):
                exclude_redundant = True if exclude_redundant == 'true' else False
            else:
                errors[line].append(
                    "exclude_redundant field with value '{}' should be either 'true' or 'false'".format(exclude_redundant))
        except IndexError:
            exclude_redundant = False

        if not errors:
            if quality == 'high':
                hq_samples.append({'name': name, 'strand': strand, 'LR': out_files,
                                   'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant})
            if quality == 'low':
                lq_samples.append({'name': name, 'strand': strand, 'LR': out_files,
                                   'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant})

    if any([len(error_list) for error_list in errors.values()]):
        print(f"File {samples.name} parsing failed, errors found:\n", file=sys.stderr)
        for line, error_list in errors.items():
            if not error_list:
                continue
            print("Line:", line, sep='\n\t', file=sys.stderr)
            print("was not parsed successfully, the following errors were found:", file=sys.stderr)
            [print("\t-", e, file=sys.stderr) for e in error_list]
        raise ValueError(f"Could not parse file {samples.name}")

    if hq_samples:
        result['ei_annotation.HQ_long_read_samples'] = hq_samples
    if lq_samples:
        result['ei_annotation.LQ_long_read_samples'] = lq_samples

    return result


def validate_paired_samples(samples):
    names = []
    errors = defaultdict(list)
    strands = ("fr-firststrand", "fr-secondstrand", "fr-unstranded")
    result = {'ei_annotation.paired_samples': []}
    first_line = True
    for line in samples:
        out_files = []
        fields = line.rstrip().split(",")
        try:
            name, strand, files, merge = fields[:4]
        except ValueError as e:
            errors[line].append(f"Unexpected input '{fields}'\n\t\tPlease make sure this is a csv file with at minimum"
                                f" the following fields name, strand, files, merge")
            break
        for file in files.split(' '):
            try:
                r1, r2 = file.split(';')
            except ValueError as e:
                if first_line:
                    errors[line].append("Unexpected input '{}'\n\t\tPlease remove this line if it is a header"
                                        ", otherwise, make sure paired reads are separated using ';'".format(file))
                else:
                    errors[line].append(
                        "Unexpected input '{}'\n\t\tPlease make sure paired reads are correctly separated using "
                        "';'".format(file))
                break
            if not os.path.exists(r1):
                errors[line].append(("File not found: '{}'".format(r1)))
            if not os.path.exists(r2):
                errors[line].append(("File not found: '{}'".format(r2)))
            out_files.append({'R1': r1, 'R2': r2})
        if name in names:
            errors[line].append(("Non-unique label specified: '{}'".format(name)))
        if strand.lower() not in strands:
            errors[line].append(("Incorrect strand '{}' specification, please choose one of {}".format(strand, strands)))

        merge = merge.lower()
        if merge in ("true", "false"):
            merge = True if merge == 'true' else False
        else:
            errors[line].append("merge field with value '{}' should be either 'true' or 'false'".format(merge))

        names.append(name)

        try:
            score = float(fields[4])
        except ValueError as e:
            errors[line].append("Cannot parse score value '{}' as float".format(fields[4]))
        except IndexError as e:
            score = 0

        try:
            is_ref = fields[5].lower()
            if is_ref in ("true", "false"):
                is_ref = True if is_ref == 'true' else False
            else:
                errors[line].append("is_ref field with value '{}' should be either 'true' or 'false'".format(is_ref))
        except IndexError:
            is_ref = False

        try:
            exclude_redundant = fields[6].lower()
            if exclude_redundant in ("true", "false"):
                exclude_redundant = True if exclude_redundant == 'true' else False
            else:
                errors[line].append(
                    "exclude_redundant field with value '{}' should be either 'true' or 'false'".format(exclude_redundant))
        except IndexError:
            exclude_redundant = False

        if not errors[line]:
            result['ei_annotation.paired_samples'].append(
                {'name': name, 'strand': strand, 'read_pair': out_files,
                 'merge': merge,
                 'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant}
            )
        first_line = False

    if any([len(error_list) for error_list in errors.values()]):
        print(f"File {samples.name} parsing failed, errors found:\n", file=sys.stderr)
        for line, error_list in errors.items():
            if not error_list:
                continue
            print("Line:", line.strip(), sep='\n\t', file=sys.stderr)
            print("was not parsed successfully, the following errors were found:", file=sys.stderr)
            [print("\t-", e, file=sys.stderr) for e in error_list]
        raise ValueError(f"Could not parse file {samples.name}")

    return result


def combine_arguments(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = parse_json_input(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources

    if cli_arguments.samples:
        for s in cli_arguments.samples:
            sample = parse_json_input(s)
            cromwell_inputs.update(sample)

    if cli_arguments.csv_paired_samples:
        cromwell_inputs.update(validate_paired_samples(cli_arguments.csv_paired_samples))

    if cli_arguments.csv_long_samples:
        cromwell_inputs.update(validate_long_samples(cli_arguments.csv_long_samples))

    cromwell_inputs["ei_annotation.reference_genome"] = cli_arguments.reference.name
    cromwell_inputs["ei_annotation.all_scoring_file"] = cli_arguments.all_scoring_file.name
    cromwell_inputs["ei_annotation.long_scoring_file"] = cli_arguments.long_scoring_file.name
    if cli_arguments.long_lq_scoring_file:
        cromwell_inputs["ei_annotation.long_lq_scoring_file"] = cli_arguments.long_lq_scoring_file.name

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
    cromwell_inputs["ei_annotation.wf_align.max_intron_len_ends"] = cli_arguments.max_intron_len_ends
    cromwell_inputs["ei_annotation.wf_align.min_identity"] = cli_arguments.min_identity

    # Separate these config files onto multiple files one for each mikado step and point the workflow to the files
    if cli_arguments.all_extra_config is not None:
        prepare_cfg, serialise_cfg, pick_cfg = separate_mikado_config(cli_arguments.all_extra_config.name, "all")
        if prepare_cfg:
            cromwell_inputs["ei_annotation.all_prepare_cfg"] = str(prepare_cfg)
        if serialise_cfg:
            cromwell_inputs["ei_annotation.all_serialise_cfg"] = str(serialise_cfg)
        if pick_cfg:
            cromwell_inputs["ei_annotation.all_pick_cfg"] = str(pick_cfg)
    if cli_arguments.long_extra_config is not None:
        prepare_cfg, serialise_cfg, pick_cfg = separate_mikado_config(cli_arguments.all_extra_config.name, "long")
        if prepare_cfg:
            cromwell_inputs["ei_annotation.long_prepare_cfg"] = str(prepare_cfg)
        if serialise_cfg:
            cromwell_inputs["ei_annotation.long_serialise_cfg"] = str(serialise_cfg)
        if pick_cfg:
            cromwell_inputs["ei_annotation.long_pick_cfg"] = str(pick_cfg)
    if cli_arguments.lq_extra_config is not None:
        prepare_cfg, serialise_cfg, pick_cfg = separate_mikado_config(cli_arguments.all_extra_config.name, "long_lq")
        if prepare_cfg:
            cromwell_inputs["ei_annotation.long_lq_prepare_cfg"] = str(prepare_cfg)
        if serialise_cfg:
            cromwell_inputs["ei_annotation.long_lq_serialise_cfg"] = str(serialise_cfg)
        if pick_cfg:
            cromwell_inputs["ei_annotation.long_lq_pick_cfg"] = str(pick_cfg)

    if cli_arguments.homology_proteins is not None:
        cromwell_inputs["ei_annotation.homology_proteins"] = cli_arguments.homology_proteins.name

    if cli_arguments.orf_calling_proteins is not None:
        cromwell_inputs["ei_annotation.orf_calling_proteins"] = cli_arguments.orf_calling_proteins.name

    if cli_arguments.orf_caller != "none":
        cromwell_inputs["ei_annotation.wf_main_mikado.orf_calling_program"] = cli_arguments.orf_caller

    if cli_arguments.separate_mikado_LQ:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = True
    else:
        cromwell_inputs["ei_annotation.wf_main_mikado.separate_LQ"] = False

    sample_validation(cromwell_inputs)

    return cromwell_inputs


def parse_json_input(s):
    try:
        result = json.load(s)
    except JSONDecodeError as e:
        lines = e.doc.split('\n')
        eprint(e)
        eprint(f"Please check '{s.name}' is a valid json file around this context:")
        error_line = e.lineno - 1
        if 'Expecting \',\' delimiter' == e.msg:
            error_line = e.lineno
        [eprint(l) for l in lines[max(0, error_line - 3):error_line]]
        eprint(' ' * (len(lines[error_line - 1]) - 1), '^', sep='')
        [eprint(l) for l in lines[error_line:min(error_line + 3, len(lines))]]
        exit(1)
    return result


def sample_validation(cromwell_inputs):
    paired_sample_names = set()
    if cromwell_inputs['ei_annotation.paired_samples']:
        for sample in cromwell_inputs['ei_annotation.paired_samples']:
            l = len(paired_sample_names)
            paired_sample_names.add(sample['name'])
            if len(paired_sample_names) == l:
                raise ValueError(f"Sample {sample['name']} is repeated, please make sure sample names are unique")
    if cromwell_inputs.get("ei_annotation.wf_align.group_to_samples", None):
        samples_in_groups = defaultdict(list)
        group_names = set()
        for group_name, group_samples in cromwell_inputs["ei_annotation.wf_align.group_to_samples"].items():
            l = len(group_names)
            group_names.add(group_name)
            if len(group_names) == l:
                raise ValueError(f"Group name {group_name} has already been used, please make sure group names are "
                                 f"unique")
            for sample_name in group_samples:
                if sample_name not in paired_sample_names:
                    raise ValueError(f"The name '{sample_name}' is not a paired_samples name {paired_sample_names}, "
                                     f"please make sure the samples in the groups have been defined previously as "
                                     f"paired samples")
                samples_in_groups[sample_name].append(group_name)
        for sample, groups in samples_in_groups.items():
            if len(groups) > 1:
                raise ValueError(f"The sample '{sample}' appears in more than one group ({groups}), please make sure "
                                 f"samples are only present in a single group")


def validate_transcriptome_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "transcriptome.schema") as schema_file:
        with open(schema_file, 'r') as schema:
            reat_schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    all_validators["is_name"] = is_valid_name
    reat_validator = validators.create(meta_schema=reat_schema, validators=all_validators)
    reat_validator(reat_schema).validate(cromwell_inputs)


def execute_cromwell(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file):
    if cli_arguments.runtime_configuration is None:
        return cromwell_submit(cli_arguments, input_parameters_filepath, workflow_options_file, wdl_file)
    else:
        return cromwell_run(input_parameters_filepath, cli_arguments.runtime_configuration.name,
                            workflow_options_file, wdl_file)


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


def cromwell_run(input_parameters_filepath, cromwell_configuration, workflow_options_file, wdl_file):
    if workflow_options_file:
        formatted_command_line = ["java", f"-Dconfig.file={cromwell_configuration}", "-jar", "cromwell.jar", "run",
                                  "-i", str(input_parameters_filepath), "-o", str(workflow_options_file),
                                  "-m", "run_details.json", str(wdl_file)]
    else:
        formatted_command_line = ["java", f"-Dconfig.file={cromwell_configuration}", "-jar", "cromwell.jar", "run",
                                  "-i", str(input_parameters_filepath), "-m", "run_details.json", str(wdl_file)]

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
                print("\n\n\nREAT Failed, the following file contains might contain information with the reasons behind"
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

    print(f"Done in {str(datetime.timedelta(seconds=time.time() - start_time))}")
    return rc


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


def validate_annotations(csv_annotation_file):
    result = {'ei_homology.annotations': []}
    errors = defaultdict(list)
    for line in csv_annotation_file:
        try:
            fasta, gff = line.strip().split(',')
        except ValueError as e:
            errors[line].append(f"Unexpected input '{line}', please make sure this is a comma-separated file with"
                                f" 2 values per line. A fasta file followed by a gff file for the associated genome")
            continue
        if not os.path.exists(fasta):
            errors[line].append(f"The fasta file '{fasta}' does not exist, please make sure this is the correct path")
        if not os.path.exists(gff):
            errors[line].append(f"The gff file '{gff}' does not exist, please make sure this is the correct path")

        if not errors[line]:
            result['ei_homology.annotations'].append({'genome': fasta, 'annotation_gff': gff})

    if any([len(error_list) for error_list in errors.values()]):
        print(f"File {csv_annotation_file.name} parsing failed, errors found:\n", file=sys.stderr)
        for line, error_list in errors.items():
            if not error_list:
                continue
            print("Line:", line.strip(), sep='\n\t', file=sys.stderr)
            print("was not parsed successfully, the following errors were found:", file=sys.stderr)
            [print("\t-", e, file=sys.stderr) for e in error_list]
        raise ValueError(f"Could not parse file {csv_annotation_file.name}")

    return result


def combine_arguments_homology(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources
    annotation = validate_annotations(cli_arguments.annotations_csv)
    cromwell_inputs.update(annotation)

    cromwell_inputs["ei_homology.genome_to_annotate"] = cli_arguments.genome.name
    cromwell_inputs["ei_homology.AlignProteins.species"] = cli_arguments.alignment_species

    # Optional extra parameters
    if cli_arguments.junction_f1_filter:
        cromwell_inputs["ei_homology.CombineResults.junction_f1_filter"] = cli_arguments.junction_f1_filter
    if cli_arguments.exon_f1_filter:
        cromwell_inputs["ei_homology.CombineResults.exon_f1_filter"] = cli_arguments.exon_f1_filter
    if cli_arguments.annotation_filters:
        cromwell_inputs["ei_homology.PrepareAnnotations.filters"] = cli_arguments.annotation_filters
    if cli_arguments.filter_min_cds:
        cromwell_inputs["ei_homology.PrepareAnnotations.min_cds_len"] = cli_arguments.filter_min_cds
    if cli_arguments.filter_max_intron:
        cromwell_inputs["ei_homology.PrepareAnnotations.max_intron_len"] = cli_arguments.filter_max_intron
    if cli_arguments.filter_min_exon:
        cromwell_inputs["ei_homology.PrepareAnnotations.min_exon_len"] = cli_arguments.filter_min_exon
    if cli_arguments.alignment_min_exon_len:
        cromwell_inputs["ei_homology.AlignProteins.min_spaln_exon_len"] = cli_arguments.alignment_min_exon_len
    if cli_arguments.alignment_filters:
        cromwell_inputs["ei_homology.AlignProteins.filters"] = cli_arguments.alignment_filters
    if cli_arguments.filter_min_cds:
        cromwell_inputs["ei_homology.AlignProteins.min_cds_len"] = cli_arguments.filter_min_cds
    if cli_arguments.filter_max_intron:
        cromwell_inputs["ei_homology.AlignProteins.max_intron_len"] = cli_arguments.filter_max_intron
    if cli_arguments.filter_min_exon:
        cromwell_inputs["ei_homology.AlignProteins.min_filter_exon_len"] = cli_arguments.filter_min_exon
    if cli_arguments.alignment_min_identity:
        cromwell_inputs["ei_homology.AlignProteins.min_identity"] = cli_arguments.alignment_min_identity
    if cli_arguments.alignment_min_coverage:
        cromwell_inputs["ei_homology.AlignProteins.min_coverage"] = cli_arguments.alignment_min_coverage
    if cli_arguments.alignment_max_per_query:
        cromwell_inputs["ei_homology.AlignProteins.max_per_query"] = cli_arguments.alignment_max_per_query
    if cli_arguments.alignment_show_intron_length:
        cromwell_inputs["ei_homology.AlignProteins.show_intron_len"] = cli_arguments.alignment_show_intron_length

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
