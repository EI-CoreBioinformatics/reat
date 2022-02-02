import glob
from datetime import datetime
import json
import os
import sys
from collections import defaultdict
from importlib import resources as pkg_resources
from json.decoder import JSONDecodeError
from pathlib import Path

from jsonschema import ValidationError, Draft7Validator, validators

from annotation.utils import symlink, link_mikado, link_assemblies, link_bams
from annotation import RUN_METADATA, report_errors, prepare_cromwell_arguments, execute_cromwell

UNSUPPORTED_GENETIC_CODES_INT = [2, 3, 4, 5, 9, 11, 13, 14, 16, 21, 22, 23, 24]

genetic_code_str_to_int = {'Universal': 1,
                           'Tetrahymena': 6,
                           'Acetabularia': 6,
                           'Ciliate': 6,
                           'Dasycladacean': 6,
                           'Hexamita': 6,
                           'Candida': 1,
                           'Euplotid': 10,
                           'SR1_Gracilibacteria': 25,
                           'Pachysolen_tannophilus': 1,
                           'Peritrich': 6}

try:
    from yaml import CLoader as Loader, CDumper as Dumper, load, dump
except ImportError:
    from yaml import Loader, Dumper, load, dump


def eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def is_valid_name(validator, value, instance, schema):
    if not isinstance(instance, str):
        yield ValidationError("%r is not a string" % instance)

    if value and set(instance).intersection("][/?\\\'\" .*$)(}{"):
        yield ValidationError("%r is not alphanumeric" % (instance,))


def separate_mikado_config(mikado_config, mikado_run):
    # Creates a new file with a timestamp if the current configuration is different from the previous.
    #  - No previous run in the same directory -> Folders created, new files with timestamps
    #  - A run with new parameters -> New files with timestamps, previous files removed
    #  - A run with the same parameters -> No changes to the existing files
    config = load(open(mikado_config, 'r'), Loader=Loader)
    prepare = config.get('prepare', None)
    serialise = config.get('serialise', None)
    pick = config.get('pick', None)

    # Load the latest version of the files to compare
    previous_prepare, prepare_path = get_previous_config_part(mikado_run, "prepare", prepare)
    previous_serialise, serialise_path = get_previous_config_part(mikado_run, "serialise", serialise)
    previous_pick, pick_path = get_previous_config_part(mikado_run, "pick", pick)

    ctimestamp = str(datetime.timestamp(datetime.now()))
    if any((prepare, serialise, pick)):
        Path(mikado_run).mkdir(exist_ok=True)
        if prepare != previous_prepare:
            prepare_path = Path(mikado_run).joinpath(ctimestamp + "-prepare.yaml")
            print(dump({'prepare': prepare}, default_flow_style=False),
                  file=open(prepare_path, 'w'))
        if serialise != previous_serialise:
            serialise_path = Path(mikado_run).joinpath(ctimestamp + "-serialise.yaml")
            print(dump({'serialise': serialise}, default_flow_style=False),
                  file=open(serialise_path, 'w'))
        if pick != previous_pick:
            pick_path = Path(mikado_run).joinpath(ctimestamp + "-pick.yaml")
            print(dump({'pick': pick}, default_flow_style=False),
                  file=open(pick_path, 'w'))

    return prepare_path, serialise_path, pick_path


def scoring_setup(cli_arguments, cromwell_inputs):
    all_scoring_filepath = cli_arguments.all_scoring_file.name
    all_scoring_path = get_prev_scoring(all_scoring_filepath)
    cromwell_inputs["ei_annotation.all_scoring_file"] = all_scoring_path

    if cli_arguments.long_scoring_file:
        long_scoring_filepath = cli_arguments.long_scoring_file.name
        long_scoring_path = get_prev_scoring(long_scoring_filepath)
        cromwell_inputs["ei_annotation.long_scoring_file"] = long_scoring_path
    else:
        cromwell_inputs["ei_annotation.long_scoring_file"] = all_scoring_path

    if cli_arguments.long_lq_scoring_file:
        long_lq_scoring_filepath = cli_arguments.long_lq_scoring_file.name
        long_lq_scoring_path = get_prev_scoring(long_lq_scoring_filepath)
        cromwell_inputs["ei_annotation.long_lq_scoring_file"] = long_lq_scoring_path


def get_prev_scoring(scoring_filepath):
    file_found = False
    prev_filepath = None
    prev_scoring = None
    scoring = load(open(scoring_filepath), Loader=Loader)
    list_of_files = glob.glob(f"scoring/*-scoring.yaml")
    for prev_filepath in sorted(list_of_files, key=os.path.getctime):
        prev_scoring = load(open(prev_filepath), Loader=Loader)
        if scoring == prev_scoring:
            file_found = True
            break
    latest_filepath = Path(prev_filepath) if prev_filepath else None
    previous_config = prev_scoring if file_found else None
    scoring_path = latest_filepath if previous_config else None
    ctimestamp = str(datetime.timestamp(datetime.now()))
    if scoring:
        Path('scoring').mkdir(exist_ok=True)
        if scoring != prev_scoring:
            scoring_path = Path('scoring').joinpath(ctimestamp + "-scoring.yaml")
            print(dump(scoring, default_flow_style=False),
                  file=open(scoring_path, 'w'))

    return str(scoring_path)


def get_previous_config_part(mikado_run, config_part, cur_config_part):
    file_found = False
    prev_filepath = None
    prev_config_part = None
    list_of_files = glob.glob(f"{mikado_run}/*-{config_part}.yaml")
    for prev_filepath in sorted(list_of_files, key=os.path.getctime):
        prev_config_part = load(open(prev_filepath), Loader=Loader).get(f"{config_part}", None)
        if cur_config_part and prev_config_part and cur_config_part == prev_config_part:
            file_found = True
            break
    latest_filepath = Path(prev_filepath) if prev_filepath else None
    previous_config = prev_config_part if file_found else None
    config_part_filepath = latest_filepath if previous_config else None

    return previous_config, config_part_filepath


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
            continue

        if name == "reference":
            errors[line].append(f"The 'reference' name is reserved for internal use, please rename the sample")

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

        exclude_redundant, is_ref, score = parse_sample_extra_fields(errors, fields, line, name, names)
        if not errors:
            if quality == 'high':
                hq_samples.append({'name': name, 'strand': strand, 'LR': out_files,
                                   'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant})
            if quality == 'low':
                lq_samples.append({'name': name, 'strand': strand, 'LR': out_files,
                                   'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant})

    report_errors(errors, samples)

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
        if name == "reference":
            errors[line].append(f"The 'reference' name is reserved for internal use, please rename the sample")
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
            errors[line].append(
                ("Incorrect strand '{}' specification, please choose one of {}".format(strand, strands)))

        merge = merge.lower()
        if merge in ("true", "false"):
            merge = True if merge == 'true' else False
        else:
            errors[line].append("merge field with value '{}' should be either 'true' or 'false'".format(merge))

        exclude_redundant, is_ref, score = parse_sample_extra_fields(errors, fields, line, name, names)

        if not errors[line]:
            result['ei_annotation.paired_samples'].append(
                {'name': name, 'strand': strand, 'read_pair': out_files,
                 'merge': merge,
                 'score': score, 'is_ref': is_ref, 'exclude_redundant': exclude_redundant}
            )
        first_line = False

    report_errors(errors, samples)

    return result


def parse_sample_extra_fields(errors, fields, line, name, names):
    names.append(name)
    score = 0
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
    return exclude_redundant, is_ref, score


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
    scoring_setup(cli_arguments, cromwell_inputs)

    cromwell_inputs["ei_annotation.wf_align.PR_hisat_extra_parameters"] = cli_arguments.PR_hisat_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.PR_star_extra_parameters"] = cli_arguments.PR_star_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.HQ_aligner"] = cli_arguments.HQ_aligner
    cromwell_inputs["ei_annotation.wf_align.LQ_aligner"] = cli_arguments.LQ_aligner

    cromwell_inputs["ei_annotation.wf_align.HQ_aligner_extra_parameters"] = cli_arguments.HQ_aligner_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.LQ_aligner_extra_parameters"] = cli_arguments.LQ_aligner_extra_parameters

    cromwell_inputs[
        "ei_annotation.wf_align.PR_stringtie_extra_parameters"] = cli_arguments.PR_stringtie_extra_parameters
    cromwell_inputs["ei_annotation.wf_align.PR_scallop_extra_parameters"] = cli_arguments.PR_scallop_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.LQ_min_identity"] = cli_arguments.LQ_min_identity
    cromwell_inputs["ei_annotation.wf_align.LQ_min_coverage"] = cli_arguments.LQ_min_coverage

    cromwell_inputs["ei_annotation.wf_align.HQ_min_identity"] = cli_arguments.HQ_min_identity
    cromwell_inputs["ei_annotation.wf_align.HQ_min_coverage"] = cli_arguments.HQ_min_coverage

    cromwell_inputs["ei_annotation.wf_align.skip_scallop"] = cli_arguments.skip_scallop
    cromwell_inputs["ei_annotation.wf_align.skip_2pass_alignment"] = cli_arguments.skip_2pass_alignment

    # Reference annotation parameters
    cromwell_inputs["ei_annotation.mode"] = cli_arguments.mode
    cromwell_inputs["ei_annotation.check_reference"] = cli_arguments.check_reference
    if cli_arguments.annotation:
        cromwell_inputs["ei_annotation.annotation"] = cli_arguments.annotation.name
        cromwell_inputs["ei_annotation.annotation_score"] = cli_arguments.annotation_score
    #####

    if cli_arguments.extra_junctions:
        cromwell_inputs["ei_annotation.wf_align.extra_junctions"] = cli_arguments.extra_junctions.name

    if cli_arguments.skip_mikado_long:
        cromwell_inputs["ei_annotation.wf_main_mikado.skip_mikado_long"] = cli_arguments.skip_mikado_long

    if cli_arguments.filter_HQ_assemblies:
        cromwell_inputs["ei_annotation.wf_main_mikado.filter_HQ_assemblies"] = cli_arguments.filter_HQ_assemblies
    if cli_arguments.filter_LQ_assemblies:
        cromwell_inputs["ei_annotation.wf_main_mikado.filter_LQ_assemblies"] = cli_arguments.filter_LQ_assemblies

    cromwell_inputs["ei_annotation.wf_align.HQ_assembler"] = cli_arguments.HQ_assembler
    cromwell_inputs["ei_annotation.wf_align.LQ_assembler"] = cli_arguments.LQ_assembler

    cromwell_inputs[
        "ei_annotation.wf_align.HQ_assembler_extra_parameters"] = cli_arguments.HQ_assembler_extra_parameters
    cromwell_inputs[
        "ei_annotation.wf_align.LQ_assembler_extra_parameters"] = cli_arguments.LQ_assembler_extra_parameters

    cromwell_inputs["ei_annotation.wf_align.min_intron_len"] = cli_arguments.min_intron_len
    cromwell_inputs["ei_annotation.wf_align.max_intron_len"] = cli_arguments.max_intron_len
    cromwell_inputs["ei_annotation.wf_align.max_intron_len_ends"] = cli_arguments.max_intron_len_ends
    cromwell_inputs["ei_annotation.wf_align.min_identity"] = cli_arguments.min_identity / 100.0

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

    if cli_arguments.exclude_LQ_junctions:
        cromwell_inputs["ei_annotation.wf_main_mikado.exclude_LQ_junctions"] = True
    else:
        cromwell_inputs["ei_annotation.wf_main_mikado.exclude_LQ_junctions"] = False

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
    sample_strand = dict()
    if cromwell_inputs.get('ei_annotation.paired_samples', None):
        for sample in cromwell_inputs['ei_annotation.paired_samples']:
            l = len(paired_sample_names)
            paired_sample_names.add(sample['name'])
            sample_strand[sample['name']] = sample['strand']
            if len(paired_sample_names) == l:
                raise ValueError(f"Sample {sample['name']} is repeated, please make sure sample names are unique")

    if cromwell_inputs.get("ei_annotation.wf_align.group_to_samples", None) and not sample_strand:
        raise ValueError(f"Sorry! It appears there are no paired samples to group. Did you mean to include the "
                         f"csv_paired_samples argument?")
    if cromwell_inputs.get("ei_annotation.wf_align.group_to_samples", None):
        seen_samples = set(paired_sample_names)
        samples_in_groups = defaultdict(list)
        group_names = set()
        for group_name, group_samples in cromwell_inputs["ei_annotation.wf_align.group_to_samples"].items():
            l = len(group_names)
            group_names.add(group_name)
            if len(group_names) == l:
                raise ValueError(f"Group name {group_name} has already been used, please make sure group names are "
                                 f"unique")
            group_strand = sample_strand[group_samples[0]]
            for sample_name in group_samples:
                if sample_name not in paired_sample_names:
                    raise ValueError(f"The name '{sample_name}' is not a paired_samples name {paired_sample_names}, "
                                     f"please make sure the samples in the groups have been defined previously as "
                                     f"paired samples")
                if group_strand != sample_strand[sample_name]:
                    raise ValueError(f"The strandness of '{sample_name}' ({sample_strand[sample_name]}) does not match "
                                     f"the strandness ({group_strand}) of other samples in this group, please ensure "
                                     f"all samples in a group have the same strandness")
                samples_in_groups[sample_name].append(group_name)
        for sample, groups in samples_in_groups.items():
            if len(groups) > 1:
                raise ValueError(f"The sample '{sample}' appears in more than one group ({groups}), please make sure "
                                 f"samples are only present in a single group")
        for sample in samples_in_groups.keys():
            seen_samples.remove(sample)
        if len(seen_samples) != 0:
            raise ValueError(f"The samples {seen_samples} do not belong to any groups, please include them in a group")


def validate_transcriptome_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "transcriptome.schema.json") as schema_file:
        with open(schema_file, 'r') as schema:
            reat_schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    all_validators["is_name"] = is_valid_name
    reat_validator = validators.create(meta_schema=reat_schema, validators=all_validators)
    reat_validator(reat_schema).validate(cromwell_inputs)


def transcriptome_cli_validation(args, reat_ap):
    if args.separate_mikado_LQ:
        if not args.long_lq_scoring_file:
            reat_ap.error("When '--separate_mikado_LQ' is enabled, --long_lq_scoring_file is required, please "
                          "provide it.")
    if args.samples and (args.csv_paired_samples or args.csv_long_samples):
        reat_ap.error("Conflicting arguments '--samples' and ['--csv_paired_samples' or '--csv_long_samples'] "
                      "provided, please choose one of csv or json sample input format")
    if not args.samples and not args.csv_paired_samples and not args.csv_long_samples:
        reat_ap.error("Please provide at least one of --samples, --csv_paired_samples, --csv_long_samples")

    genetic_code, mikado_genetic_code = validate_genetic_code(args, reat_ap)
    return genetic_code, mikado_genetic_code


def validate_genetic_code(args, reat_ap):
    genetic_code = 0
    mikado_genetic_code = 0
    if args.orf_caller == 'transdecoder':
        if args.genetic_code != '0':
            try:
                mikado_genetic_code = genetic_code_str_to_int[args.genetic_code]
                genetic_code = args.genetic_code
            except KeyError:
                reat_ap.error(
                    f"Sorry, genetic_code={args.genetic_code} is not supported when orf caller is TransDecoder, please "
                    f"use one of {', '.join(genetic_code_str_to_int.keys())}")
        else:
            genetic_code = 'Universal'
            mikado_genetic_code = 0
    elif args.orf_caller == 'prodigal':
        if args.genetic_code != '0' and not args.genetic_code.isdigit():
            try:
                mikado_genetic_code = genetic_code_str_to_int[args.genetic_code]
                genetic_code = args.genetic_code
            except KeyError:
                reat_ap.error(
                    f"Sorry, genetic_code={args.genetic_code} is not supported, please "
                    f"use one of {', '.join(genetic_code_str_to_int.keys())}")
        else:
            genetic_code = int(args.genetic_code)
            if genetic_code in UNSUPPORTED_GENETIC_CODES_INT or genetic_code > 25:
                reat_ap.error(
                    f"Sorry, genetic_code={args.genetic_code} is not supported, please use one of "
                    f"{set(range(0, 25)) - set(UNSUPPORTED_GENETIC_CODES_INT)}")
            mikado_genetic_code = genetic_code
        if int(genetic_code) == 0:
            genetic_code = 1
    return genetic_code, mikado_genetic_code


def transcriptome_module(cli_arguments, genetic_code='0', mikado_genetic_code=0):
    """
    Collects the CLI arguments and combines them with CLI defined input files.
    The resulting object is validated and the final inputs are written to a json Cromwell input file.
    :param cli_arguments: Validated CLI values with un-inspected files
    :return: Return code from Cromwell
    """
    # Print input file for cromwell
    cromwell_inputs = combine_arguments(cli_arguments)
    cromwell_inputs['ei_annotation.mikado_genetic_code'] = mikado_genetic_code
    cromwell_inputs['ei_annotation.genetic_code'] = str(genetic_code)
    cromwell_jar, runtime_config = prepare_cromwell_arguments(cli_arguments)

    # Validate input against schema
    validate_transcriptome_inputs(cromwell_inputs)
    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.transcriptome_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        rc = execute_cromwell(runtime_config, cromwell_jar,
                              cli_arguments.output_parameters_file, workflow_options_file, wdl_file)
        if rc == 0:
            collect_transcriptome_output(RUN_METADATA)
        return rc


def collect_transcriptome_output(RUN_METADATA, output_path="outputs"):
    # Get the outputs and symlink them to the output folder
    run_metadata = json.load(open(RUN_METADATA))
    outputs = run_metadata["outputs"]
    prfx = "ei_annotation."
    outputs_path = output_path
    if not os.path.exists(outputs_path):
        os.mkdir(outputs_path)

    # reat/outputs
    # ├── align_long.HQ_read_samples.summary.stats.tsv
    # ├── alignments/
    # ├── align_short.SR_read_samples.summary.stats.tsv
    # ├── assembly_long/
    # ├── assembly_short/
    # ├── assembly_short.scallop.summary.stats.tsv
    # ├── assembly_short.stringtie.summary.stats.tsv
    # ├── plots/
    # ├── Mikado_long-permissive.loci.gff3
    # ├── Mikado_long-permissive.loci.gff3.stats
    # ├── Mikado_long-permissive.loci.metrics.tsv
    # ├── Mikado_long-permissive.loci.scores.tsv
    # ├── Mikado_short_and_long-permissive.loci.gff3
    # ├── Mikado_short_and_long-permissive.loci.gff3.stats
    # ├── Mikado_short_and_long-permissive.loci.metrics.tsv
    # ├── Mikado_short_and_long-permissive.loci.scores.tsv
    # ├── mikado.summary.stats.tsv
    # └── portcullis/

    # Alignments
    # Array[AlignedSample]? SR_bams = wf_align.SR_bams
    # Array[File]? SR_alignment_summary_stats = wf_align.SR_summary_stats
    # File? SR_alignment_summary_stats_table = wf_align.SR_summary_stats_table
    link_bams(outputs, outputs_path, prfx + 'SR_bams', prfx + 'SR_alignment_summary_stats',
              prfx + 'SR_alignment_summary_stats_table')
    # Array[AlignedSample]? LQ_bams = wf_align.LQ_bams
    # Array[File]? LQ_alignment_summary_stats = wf_align.LQ_summary_stats
    # File? LQ_alignment_summary_stats_table = wf_align.LQ_summary_stats_table
    link_bams(outputs, outputs_path, prfx + 'LQ_bams', prfx + 'LQ_alignment_summary_stats',
              prfx + 'LQ_alignment_summary_stats_table')
    # Array[AlignedSample]? HQ_bams = wf_align.HQ_bams
    # Array[File]? HQ_alignment_summary_stats = wf_align.HQ_summary_stats
    # File? HQ_alignment_summary_stats_table = wf_align.HQ_summary_stats_table
    link_bams(outputs, outputs_path, prfx + 'HQ_bams', prfx + 'HQ_alignment_summary_stats',
              prfx + 'HQ_alignment_summary_stats_table')

    # assembly_short
    # Array[AssembledSample]? SR_asms = wf_align.SR_gff
    # Array[File]? SR_assembly_stats = wf_align.SR_assembly_stats
    link_assemblies(prfx + 'SR_asms', os.path.join(outputs_path, 'assembly_short'), prfx + 'SR_assembly_stats', outputs)

    # assembly_long
    # Array[AssembledSample]? LQ_asms = wf_align.LQ_gff
    # Array[File]? LQ_assembly_stats = wf_align.LQ_assembly_stats
    link_assemblies(prfx + 'LQ_asms', os.path.join(outputs_path, 'assembly_long'), prfx + 'LQ_assembly_stats', outputs)
    # Array[AssembledSample]? HQ_asms = wf_align.HQ_gff
    # Array[File]? HQ_assembly_stats = wf_align.HQ_assembly_stats
    link_assemblies(prfx + 'HQ_asms', os.path.join(outputs_path, 'assembly_long'), prfx + 'HQ_assembly_stats', outputs)

    # portcullis
    if any((outputs[prfx + 'portcullis_pass_tab'],
            outputs[prfx + 'portcullis_pass_bed'],
            outputs[prfx + 'portcullis_pass_gff3'],
            outputs[prfx + 'portcullis_fail_tab'],
            outputs[prfx + 'portcullis_fail_bed'],
            outputs[prfx + 'portcullis_fail_gff3'])):
        portcullis_path = os.path.join(outputs_path, 'portcullis')
        os.mkdir(portcullis_path) if not os.path.exists(portcullis_path) else ""
        # File? portcullis_pass_tab = wf_align.pass_filtered_tab
        if outputs[prfx + 'portcullis_pass_tab']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_pass_tab'])
        # File? portcullis_pass_bed = wf_align.pass_filtered_bed
        if outputs[prfx + 'portcullis_pass_bed']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_pass_bed'])
        # File? portcullis_pass_gff3 = wf_align.pass_filtered_gff3
        if outputs[prfx + 'portcullis_pass_gff3']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_pass_gff3'])
        # File? portcullis_fail_tab = wf_align.fail_filtered_tab
        if outputs[prfx + 'portcullis_fail_tab']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_fail_tab'])
        # File? portcullis_fail_bed = wf_align.fail_filtered_bed
        if outputs[prfx + 'portcullis_fail_bed']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_fail_bed'])
        # File? portcullis_fail_gff3 = wf_align.fail_filtered_gff3
        if outputs[prfx + 'portcullis_fail_gff3']:
            symlink(portcullis_path, outputs[prfx + 'portcullis_fail_gff3'])

    # TODO
    #  Plots
    #  Array[Array[File]]? stats = wf_align.stats
    #  Array[Array[File]]? actg_cycles_plots = wf_align.actg_cycles_plots
    #  Array[Array[File]]? coverage_plots = wf_align.coverage_plots
    #  Array[Array[File]]? gc_content_plots = wf_align.gc_content_plots
    #  Array[Array[File]]? gc_depth_plots = wf_align.gc_depth_plots
    #  Array[Array[File]]? htmls = wf_align.htmls
    #  Array[Array[File]]? indel_cycles_plots = wf_align.indel_cycles_plots
    #  Array[Array[File]]? indel_dist_plots = wf_align.indel_dist_plots
    #  Array[Array[File]]? insert_size_plots = wf_align.insert_size_plots
    #  Array[Array[File]]? quals_plots = wf_align.quals_plots
    #  Array[Array[File]]? quals2_plots = wf_align.quals2_plots
    #  Array[Array[File]]? quals3_plots = wf_align.quals3_plots
    #  Array[Array[File]]? quals_hm_plots = wf_align.quals_hm_plots

    # File? SR_stringtie_summary_stats = wf_align.SR_stringtie_summary_stats
    if outputs[prfx + 'SR_stringtie_summary_stats']:
        symlink(outputs_path, outputs[prfx + 'SR_stringtie_summary_stats'])
    # File? SR_scallop_summary_stats = wf_align.SR_scallop_summary_stats
    if outputs[prfx + 'SR_scallop_summary_stats']:
        symlink(outputs_path, outputs[prfx + 'SR_scallop_summary_stats'])
    # File? LQ_assembly_summary_stats = wf_align.LQ_assembly_summary_stats
    if outputs[prfx + 'LQ_assembly_summary_stats']:
        symlink(outputs_path, outputs[prfx + 'LQ_assembly_summary_stats'])
    # File? HQ_assembly_summary_stats = wf_align.HQ_assembly_summary_stats
    if outputs[prfx + 'HQ_assembly_summary_stats']:
        symlink(outputs_path, outputs[prfx + 'HQ_assembly_summary_stats'])

    # File? mikado_long_loci = wf_main_mikado.long_loci
    # File? mikado_long_scores = wf_main_mikado.long_scores
    # File? mikado_long_metrics = wf_main_mikado.long_metrics
    # File? mikado_long_stats = wf_main_mikado.long_stats
    link_mikado(outputs, outputs_path, prfx + 'mikado_long_')
    #
    # File? mikado_short_loci = wf_main_mikado.short_loci
    # File? mikado_short_scores = wf_main_mikado.short_scores
    # File? mikado_short_metrics = wf_main_mikado.short_metrics
    # File? mikado_short_stats = wf_main_mikado.short_stats
    link_mikado(outputs, outputs_path, prfx + 'mikado_short_')
    #
    # File? mikado_short_and_long_noLQ_loci = wf_main_mikado.short_and_long_noLQ_loci
    # File? mikado_short_and_long_noLQ_scores = wf_main_mikado.short_and_long_noLQ_scores
    # File? mikado_short_and_long_noLQ_metrics = wf_main_mikado.short_and_long_noLQ_metrics
    # File? mikado_short_and_long_noLQ_stats = wf_main_mikado.short_and_long_noLQ_stats
    link_mikado(outputs, outputs_path, prfx + 'mikado_short_and_long_noLQ_')
    #
    # File? mikado_longHQ_loci = wf_main_mikado.longHQ_loci
    # File? mikado_longHQ_scores = wf_main_mikado.longHQ_scores
    # File? mikado_longHQ_metrics = wf_main_mikado.longHQ_metrics
    # File? mikado_longHQ_stats = wf_main_mikado.longHQ_stats
    link_mikado(outputs, outputs_path, prfx + 'mikado_longHQ_')
    #
    # File? mikado_longLQ_loci = wf_main_mikado.longLQ_loci
    # File? mikado_longLQ_scores = wf_main_mikado.longLQ_scores
    # File? mikado_longLQ_metrics = wf_main_mikado.longLQ_metrics
    # File? mikado_longLQ_stats = wf_main_mikado.longLQ_stats
    link_mikado(outputs, outputs_path, prfx + 'mikado_longLQ_')

    # File mikado_summary_stats = wf_main_mikado.mikado_stats_summary
    if outputs[prfx + 'mikado_summary_stats']:
        symlink(outputs_path, outputs[prfx + 'mikado_summary_stats'])