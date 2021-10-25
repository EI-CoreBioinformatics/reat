import json
import os
from collections import defaultdict
from importlib import resources as pkg_resources

from jsonschema import Draft7Validator, validators

from annotation.utils import symlink
from annotation import RUN_METADATA, report_errors, prepare_cromwell_arguments, execute_cromwell


def combine_arguments_homology(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources
    annotation = validate_annotations(cli_arguments.annotations_csv)
    cromwell_inputs.update(annotation)

    cromwell_inputs["ei_homology.genome_to_annotate"] = cli_arguments.genome.name
    cromwell_inputs["ei_homology.species"] = cli_arguments.alignment_species
    cromwell_inputs["ei_homology.mikado_config"] = cli_arguments.mikado_config.name
    cromwell_inputs["ei_homology.mikado_scoring"] = cli_arguments.mikado_scoring.name
    cromwell_inputs["ei_homology.output_prefix"] = cli_arguments.output_prefix

    if cli_arguments.protein_sequences:
        cromwell_inputs["ei_homology.protein_sequence_files"] = cli_arguments.protein_sequences

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

    if cli_arguments.filter_min_cds:
        cromwell_inputs["ei_homology.PrepareAlignments.min_cds_len"] = cli_arguments.filter_min_cds
    if cli_arguments.alignment_filters:
        cromwell_inputs["ei_homology.PrepareAlignments.filters"] = cli_arguments.alignment_filters
    if cli_arguments.filter_min_exon:
        cromwell_inputs["ei_homology.PrepareAlignments.min_filter_exon_len"] = cli_arguments.filter_min_exon
    if cli_arguments.alignment_show_intron_length:
        cromwell_inputs["ei_homology.PrepareAlignments.show_intron_len"] = cli_arguments.alignment_show_intron_length
    if cli_arguments.filter_max_intron:
        cromwell_inputs["ei_homology.PrepareAlignments.max_intron_len"] = cli_arguments.filter_max_intron

    if cli_arguments.alignment_min_exon_len:
        cromwell_inputs["ei_homology.AlignProteins.min_spaln_exon_len"] = cli_arguments.alignment_min_exon_len
    if cli_arguments.alignment_min_identity:
        cromwell_inputs["ei_homology.AlignProteins.min_identity"] = cli_arguments.alignment_min_identity
    if cli_arguments.alignment_min_coverage:
        cromwell_inputs["ei_homology.AlignProteins.min_coverage"] = cli_arguments.alignment_min_coverage
    if cli_arguments.alignment_max_per_query:
        cromwell_inputs["ei_homology.AlignProteins.max_per_query"] = cli_arguments.alignment_max_per_query
    if cli_arguments.alignment_recursion_level:
        cromwell_inputs["ei_homology.AlignProteins.recursion_level"] = cli_arguments.alignment_recursion_level

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

    report_errors(errors, csv_annotation_file)
    return result


def collect_homology_output(run_metadata):
    # Get the outputs and symlink them to the output folder
    run_metadata = json.load(open(run_metadata))
    outputs = run_metadata['outputs']
    outputs_path = 'outputs'
    if not os.path.exists(outputs_path):
        os.mkdir(outputs_path)

    # ├── Annotations
    # │      ├── Homo_sapiens.coding.annotation.stats
    # │      ├── Homo_sapiens.coding.clean.extra_attr.gff
    # │      ├── Monodelphis_domestica.coding.annotation.stats
    # │      ├── Monodelphis_domestica.coding.clean.extra_attr.gff
    # │      ├── Notamacropus_eugenii.coding.annotation.stats
    # │      ├── Notamacropus_eugenii.coding.clean.extra_attr.gff
    # │      ├── Sarcophilus_harrisii.coding.annotation.stats
    # │      └── Sarcophilus_harrisii.coding.clean.extra_attr.gff
    # ├── Homo_sapiens.alignment.stats
    # ├── Homo_sapiens.alignment.stop_extended.extra_attr.mgc.xspecies_scores.gff
    # ├── Monodelphis_domestica.alignment.stats
    # ├── Monodelphis_domestica.alignment.stop_extended.extra_attr.mgc.xspecies_scores.gff
    # ├── Notamacropus_eugenii.alignment.stats
    # ├── Notamacropus_eugenii.alignment.stop_extended.extra_attr.mgc.xspecies_scores.gff
    # ├── Sarcophilus_harrisii.alignment.stats
    # ├── Sarcophilus_harrisii.alignment.stop_extended.extra_attr.mgc.xspecies_scores.gff
    # ├── ScoreAlignments
    # │      ├── all_avgF1.bin.txt
    # │      ├── comp_Homo_sapiens.alignment.stop_extended.extra_attr_Homo_sapiens_detail.tab
    # │      ├── comp_Homo_sapiens.alignment.stop_extended.extra_attr_Homo_sapiens.tab
    # │      ├── comp_Monodelphis_domestica.alignment.stop_extended.extra_attr_Monodelphis_domestica_detail.tab
    # │      ├── comp_Monodelphis_domestica.alignment.stop_extended.extra_attr_Monodelphis_domestica.tab
    # │      ├── comp_Notamacropus_eugenii.alignment.stop_extended.extra_attr_Notamacropus_eugenii_detail.tab
    # │      ├── comp_Notamacropus_eugenii.alignment.stop_extended.extra_attr_Notamacropus_eugenii.tab
    # │      ├── comp_Sarcophilus_harrisii.alignment.stop_extended.extra_attr_Sarcophilus_harrisii_detail.tab
    # │      └── comp_Sarcophilus_harrisii.alignment.stop_extended.extra_attr_Sarcophilus_harrisii.tab
    # ├── xspecies.loci.gff3
    # ├── xspecies.loci.gff3.stats
    # ├── xspecies.loci.metrics.tsv
    # └── xspecies.loci.scores.tsv

    # Annotations
    annotations_path = os.path.join(outputs_path, 'ei_homology.annotations')
    if not os.path.exists(annotations_path):
        os.mkdir(annotations_path)
    # Array[File] clean_annotations = PrepareAnnotations.cleaned_up_gff
    for clean_annotation in outputs['ei_homology.clean_annotations']:
        symlink(annotations_path, clean_annotation)
    # Array[File] annotation_filter_stats = PrepareAnnotations.stats
    for annotation_stats in outputs['ei_homology.annotation_filter_stats']:
        symlink(annotations_path, annotation_stats)

    # ScoreAlignments
    score_alignments_path = os.path.join(outputs_path, 'ei_homology.score_alignments')
    if not os.path.exists(score_alignments_path):
        os.mkdir(score_alignments_path)
    # Array[File] mgc_evaluation = ScoreAlignments.alignment_compare
    for mgc_eval in outputs['ei_homology.mgc_evaluation']:
        symlink(score_alignments_path, mgc_eval)
    # Array[File] mgc_evaluation_detail = ScoreAlignments.alignment_compare_detail
    for mgc_eval_detail in outputs['ei_homology.mgc_evaluation_detail']:
        symlink(score_alignments_path, mgc_eval_detail)
    # File        mgc_score_summary = ScoreSummary.summary_table
    symlink(score_alignments_path, outputs['ei_homology.mgc_score_summary'])

    # Main output folder
    # Array[File] xspecies_combined_alignments = CombineXspecies.xspecies_scored_alignment
    for xspc_combined_aln in outputs['ei_homology.xspecies_combined_alignments']:
        symlink(outputs_path, xspc_combined_aln)
    # Array[File] alignment_filter_stats = AlignProteins.stats
    for aln_filter_stat in outputs['ei_homology.alignment_filter_stats']:
        symlink(outputs_path, aln_filter_stat)

    # File loci = MikadoPick.loci
    symlink(outputs_path, outputs['ei_homology.loci'])
    # File scores = MikadoPick.scores
    symlink(outputs_path, outputs['ei_homology.scores'])
    # File metrics = MikadoPick.metrics
    symlink(outputs_path, outputs['ei_homology.metrics'])
    # File stats = MikadoPick.stats
    symlink(outputs_path, outputs['ei_homology.stats'])


def homology_module(cli_arguments):
    cromwell_inputs = combine_arguments_homology(cli_arguments)
    validate_homology_inputs(cromwell_inputs)

    cromwell_jar, runtime_config = prepare_cromwell_arguments(cli_arguments)

    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.homology_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        rc = execute_cromwell(runtime_config, cromwell_jar,
                              cli_arguments.output_parameters_file, workflow_options_file, wdl_file)
        if rc == 0:
            collect_homology_output(RUN_METADATA)

        return rc
