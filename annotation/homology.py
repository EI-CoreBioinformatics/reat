import json
import os
import sys
from collections import defaultdict
from importlib import resources as pkg_resources

from jsonschema import Draft7Validator, validators


def combine_arguments_homology(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources
    annotation = validate_annotations(cli_arguments.annotations_csv)
    cromwell_inputs.update(annotation)

    cromwell_inputs["ei_homology.genome_to_annotate"] = cli_arguments.genome.name
    cromwell_inputs["ei_homology.AlignProteins.species"] = cli_arguments.alignment_species
    cromwell_inputs["ei_homology.mikado_config"] = cli_arguments.mikado_config.name
    cromwell_inputs["ei_homology.mikado_scoring"] = cli_arguments.mikado_scoring.name

    # Optional extra parameters
    if cli_arguments.pick_extra_config:
        cromwell_inputs['ei_homology.MikadoPick.extra_config'] = cli_arguments.pick_extra_config.name
    if cli_arguments.junctions:
        cromwell_inputs['ei_homology.Mikado.junctions'] = cli_arguments.junctions.name
    if cli_arguments.utrs:
        cromwell_inputs['ei_homology.Mikado.utrs'] = cli_arguments.utrs.name
    if cli_arguments.min_cdna_length:
        cromwell_inputs["ei_homology.Mikado.min_cdna_length"] = cli_arguments.min_cdna_length
    if cli_arguments.max_intron_length:
        cromwell_inputs["ei_homology.Mikado.max_intron_length"] = cli_arguments.max_intron_length

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
    with pkg_resources.path("validation", "homology.schema.json") as schema_file:
        with open(schema_file, 'r') as schema:
            schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    reat_validator = validators.create(meta_schema=schema, validators=all_validators)
    reat_validator(schema).validate(cromwell_inputs)


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