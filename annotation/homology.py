import json
import os
from collections import defaultdict
from importlib import resources as pkg_resources

from jsonschema import Draft7Validator, validators

from annotation.utils import symlink
from annotation import (
    RUN_METADATA,
    report_errors,
    prepare_cromwell_arguments,
    execute_cromwell,
)


def combine_arguments_homology(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources

    if cli_arguments.annotations_csv:
        annotation = validate_annotations(cli_arguments.annotations_csv)
        cromwell_inputs.update(annotation)

    cromwell_inputs["ei_homology.genome_to_annotate"] = cli_arguments.genome.name
    cromwell_inputs["ei_homology.species"] = cli_arguments.alignment_species
    cromwell_inputs["ei_homology.codon_table"] = cli_arguments.codon_table
    cromwell_inputs["ei_homology.mikado_config"] = cli_arguments.mikado_config.name
    cromwell_inputs["ei_homology.mikado_scoring"] = cli_arguments.mikado_scoring.name
    cromwell_inputs["ei_homology.output_prefix"] = cli_arguments.output_prefix
    # cromwell_inputs["ei_homology.aligner"] = cli_arguments.aligner
    # if cli_arguments.aligner == "spaln":
    #     cromwell_inputs["ei_homology.do_spaln"] = True
    # if cli_arguments.aligner == "miniprot":
    #     cromwell_inputs["ei_homology.do_miniprot"] = True
    # if cli_arguments.aligner == "both":
    #     cromwell_inputs["ei_homology.do_spaln"] = True
    #     cromwell_inputs["ei_homology.do_miniprot"] = True

    # check if -Q and -A are present in the extra parameters
    if cli_arguments.spaln_extra_parameters:
        if "-Q" not in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please provide the -Q flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-A" not in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please provide the -A flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
    # exit if options -t, -K, -O, -Q, -M, -T, -d, -o, -yL are present in the extra parameters
    if cli_arguments.spaln_extra_parameters:
        if "-t" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -t flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-K" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -K flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-O" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -O flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-M" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -M flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-T" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -T flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-d" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -d flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-o" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -o flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
        if "-yL" in cli_arguments.spaln_extra_parameters:
            raise ValueError(
                "Please do not provide the -yL flag in the SPALN extra parameters (--spaln_extra_parameters)"
            )
    cromwell_inputs["ei_homology.spaln_extra_parameters"] = (
        cli_arguments.spaln_extra_parameters
    )

    if cli_arguments.protein_sequences and len(cli_arguments.protein_sequences) > 1:
        cromwell_inputs["ei_homology.protein_sequence_files"] = (
            cli_arguments.protein_sequences
        )
    else:
        if not cli_arguments.annotations_csv:
            import sys

            print(
                f"Please use more than a single cross species protein sequence file\n"
                f"{cli_arguments.protein_sequences}",
                file=sys.stderr,
            )
            exit(1)

    # Optional extra parameters
    if cli_arguments.pick_extra_config:
        cromwell_inputs["ei_homology.MikadoPick.extra_config"] = (
            cli_arguments.pick_extra_config.name
        )
    if cli_arguments.junctions:
        cromwell_inputs["ei_homology.Mikado.junctions"] = cli_arguments.junctions.name
    if cli_arguments.utrs:
        cromwell_inputs["ei_homology.Mikado.utrs"] = cli_arguments.utrs.name
    if cli_arguments.min_cdna_length:
        cromwell_inputs["ei_homology.Mikado.min_cdna_length"] = (
            cli_arguments.min_cdna_length
        )
    if cli_arguments.max_intron_length:
        cromwell_inputs["ei_homology.Mikado.max_intron_length"] = (
            cli_arguments.max_intron_length
        )

    if cli_arguments.junction_f1_filter:
        cromwell_inputs["ei_homology.CombineResultsSpaln.junction_f1_filter"] = (
            cli_arguments.junction_f1_filter
        )
        cromwell_inputs["ei_homology.CombineResultsMiniprot.junction_f1_filter"] = (
            cli_arguments.junction_f1_filter
        )
    if cli_arguments.exon_f1_filter:
        cromwell_inputs["ei_homology.CombineResultsSpaln.exon_f1_filter"] = (
            cli_arguments.exon_f1_filter
        )
        cromwell_inputs["ei_homology.CombineResultsMiniprot.exon_f1_filter"] = (
            cli_arguments.exon_f1_filter
        )
    if cli_arguments.annotation_filters:
        cromwell_inputs["ei_homology.PrepareAnnotations.filters"] = (
            cli_arguments.annotation_filters
        )
    if cli_arguments.filter_min_cds:
        cromwell_inputs["ei_homology.PrepareAnnotations.min_cds_len"] = (
            cli_arguments.filter_min_cds
        )
    if cli_arguments.filter_max_intron:
        cromwell_inputs["ei_homology.PrepareAnnotations.max_intron_len"] = (
            cli_arguments.filter_max_intron
        )
    if cli_arguments.filter_min_exon:
        cromwell_inputs["ei_homology.PrepareAnnotations.min_exon_len"] = (
            cli_arguments.filter_min_exon
        )

    if cli_arguments.filter_min_cds:
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsSpaln.min_cds_len"] = (
            cli_arguments.filter_min_cds
        )
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsMiniprot.min_cds_len"] = (
            cli_arguments.filter_min_cds
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsSpaln.min_cds_len"] = (
            cli_arguments.filter_min_cds
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsMiniprot.min_cds_len"] = (
            cli_arguments.filter_min_cds
        )
    if cli_arguments.alignment_filters:
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsSpaln.filters"] = (
            cli_arguments.alignment_filters
        )
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsMiniprot.filters"] = (
            cli_arguments.alignment_filters
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsSpaln.filters"] = (
            cli_arguments.alignment_filters
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsMiniprot.filters"] = (
            cli_arguments.alignment_filters
        )
    if cli_arguments.filter_min_exon:
        cromwell_inputs[
            "ei_homology.PrepareProteinAlignmentsSpaln.min_filter_exon_len"
        ] = cli_arguments.filter_min_exon
        cromwell_inputs[
            "ei_homology.PrepareProteinAlignmentsMiniprot.min_filter_exon_len"
        ] = cli_arguments.filter_min_exon
        cromwell_inputs["ei_homology.PrepareAlignmentsSpaln.min_filter_exon_len"] = (
            cli_arguments.filter_min_exon
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsMiniprot.min_filter_exon_len"] = (
            cli_arguments.filter_min_exon
        )
    if cli_arguments.alignment_show_intron_length:
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsSpaln.show_intron_len"] = (
            cli_arguments.alignment_show_intron_length
        )
        cromwell_inputs[
            "ei_homology.PrepareProteinAlignmentsMiniprot.show_intron_len"
        ] = cli_arguments.alignment_show_intron_length
        cromwell_inputs["ei_homology.PrepareAlignmentsSpaln.show_intron_len"] = (
            cli_arguments.alignment_show_intron_length
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsMiniprot.show_intron_len"] = (
            cli_arguments.alignment_show_intron_length
        )
    if cli_arguments.filter_max_intron:
        cromwell_inputs["ei_homology.PrepareProteinAlignmentsSpaln.max_intron_len"] = (
            cli_arguments.filter_max_intron
        )
        cromwell_inputs[
            "ei_homology.PrepareProteinAlignmentsMiniprot.max_intron_len"
        ] = cli_arguments.filter_max_intron
        cromwell_inputs["ei_homology.PrepareAlignmentsSpaln.max_intron_len"] = (
            cli_arguments.filter_max_intron
        )
        cromwell_inputs["ei_homology.PrepareAlignmentsMiniprot.max_intron_len"] = (
            cli_arguments.filter_max_intron
        )
    if cli_arguments.alignment_max_intron_length:
        cromwell_inputs["ei_homology.AlignProteinsMiniprot.max_intron_length"] = (
            cli_arguments.alignment_max_intron_length
        )
    if cli_arguments.alignment_min_intron_length:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.min_spaln_intron_length"] = (
            cli_arguments.alignment_min_intron_length
        )
    if cli_arguments.alignment_min_identity:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.min_identity"] = (
            cli_arguments.alignment_min_identity
        )
        cromwell_inputs["ei_homology.AlignProteinsMiniprot.min_identity"] = (
            cli_arguments.alignment_min_identity
        )
    if cli_arguments.alignment_min_coverage:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.min_coverage"] = (
            cli_arguments.alignment_min_coverage
        )
        cromwell_inputs["ei_homology.AlignProteinsMiniprot.min_coverage"] = (
            cli_arguments.alignment_min_coverage
        )
    if cli_arguments.alignment_max_per_query:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.max_per_query"] = (
            cli_arguments.alignment_max_per_query
        )
        cromwell_inputs["ei_homology.AlignProteinsMiniprot.max_per_query"] = (
            cli_arguments.alignment_max_per_query
        )
    if cli_arguments.post_alignment_clip:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.clip"] = (
            cli_arguments.post_alignment_clip
        )
    if cli_arguments.term5c_len:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.term5c_len"] = (
            cli_arguments.term5c_len
        )
    if cli_arguments.term5i_len:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.term5i_len"] = (
            cli_arguments.term5i_len
        )
    if cli_arguments.term3c_len:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.term3c_len"] = (
            cli_arguments.term3c_len
        )
    if cli_arguments.term3i_len:
        cromwell_inputs["ei_homology.AlignProteinsSpaln.term3i_len"] = (
            cli_arguments.term3i_len
        )
    return cromwell_inputs


def validate_homology_inputs(cromwell_inputs):
    with pkg_resources.path("validation", "homology.schema.json") as schema_file:
        with open(schema_file, "r") as schema:
            schema = json.load(schema)
    all_validators = dict(Draft7Validator.VALIDATORS)
    reat_validator = validators.create(meta_schema=schema, validators=all_validators)
    reat_validator(schema).validate(cromwell_inputs)


def validate_annotations(csv_annotation_file):
    result = {"ei_homology.annotations": []}
    errors = defaultdict(list)
    for line in csv_annotation_file:
        try:
            fasta, gff = line.strip().split(",")
        except ValueError as e:
            errors[line].append(
                f"Unexpected input '{line}', please make sure this is a comma-separated file with"
                f" 2 values per line. A fasta file followed by a gff file for the associated genome"
            )
            continue
        if not os.path.exists(fasta):
            errors[line].append(
                f"The fasta file '{fasta}' does not exist, please make sure this is the correct path"
            )
        if not os.path.exists(gff):
            errors[line].append(
                f"The gff file '{gff}' does not exist, please make sure this is the correct path"
            )

        if not errors[line]:
            result["ei_homology.annotations"].append(
                {"genome": fasta, "annotation_gff": gff}
            )

    report_errors(errors, csv_annotation_file)
    return result


def collect_homology_output(run_metadata):
    # Get the outputs and symlink them to the output folder
    run_metadata = json.load(open(run_metadata))
    outputs = run_metadata["outputs"]
    outputs_path = "outputs"
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
    # │      ├── Junction_precision_binning.tsv
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

    # cleaned up annotations
    annotations_path = os.path.join(outputs_path, "ei_homology.clean_input_annotations")
    if not os.path.exists(annotations_path):
        os.mkdir(annotations_path)
    # Array[File] clean_annotations = PrepareAnnotations.cleaned_up_gff
    if outputs["ei_homology.clean_annotations"]:
        for clean_annotation in outputs["ei_homology.clean_annotations"]:
            symlink(annotations_path, clean_annotation)
    if outputs["ei_homology.clean_annotations_proteins"]:
        for clean_annotations_proteins in outputs[
            "ei_homology.clean_annotations_proteins"
        ]:
            symlink(annotations_path, clean_annotations_proteins)
    # Array[File] annotation_filter_stats = PrepareAnnotations.stats
    if outputs["ei_homology.annotation_filter_stats"]:
        for annotation_stats in outputs["ei_homology.annotation_filter_stats"]:
            symlink(annotations_path, annotation_stats)

    # Intermediates and Alignments
    intermediates_path_spaln = os.path.join(
        outputs_path, "ei_homology.intermediates", "spaln"
    )
    if not os.path.exists(intermediates_path_spaln):
        os.makedirs(intermediates_path_spaln)
    if outputs["ei_homology.spaln_intermediates"]:
        for intermediate in outputs["ei_homology.spaln_intermediates"]:
            symlink(intermediates_path_spaln, intermediate)
    intermediates_path_miniprot = os.path.join(
        outputs_path, "ei_homology.intermediates", "miniprot"
    )
    if not os.path.exists(intermediates_path_miniprot):
        os.makedirs(intermediates_path_miniprot)
    if outputs["ei_homology.miniprot_intermediates"]:
        for intermediate in outputs["ei_homology.miniprot_intermediates"]:
            symlink(intermediates_path_miniprot, intermediate)

    # ScoreAlignments
    score_alignments_path = os.path.join(outputs_path, "ei_homology.score_alignments")
    if not os.path.exists(score_alignments_path):
        os.mkdir(score_alignments_path)
    # Array[File] mgc_evaluation = ScoreAlignments.alignment_compare
    # spaln
    if outputs["ei_homology.mgc_evaluation_spaln"]:
        for mgc_eval in outputs["ei_homology.mgc_evaluation_spaln"]:
            symlink(score_alignments_path, mgc_eval)
    # Array[File] mgc_evaluation_detail = ScoreAlignments.alignment_compare_detail
    if outputs["ei_homology.mgc_evaluation_spaln_detail"]:
        for mgc_eval_detail in outputs["ei_homology.mgc_evaluation_spaln_detail"]:
            symlink(score_alignments_path, mgc_eval_detail)
    # miniprot
    if outputs["ei_homology.mgc_evaluation_miniprot"]:
        for mgc_eval in outputs["ei_homology.mgc_evaluation_miniprot"]:
            symlink(score_alignments_path, mgc_eval)
    if outputs["ei_homology.mgc_evaluation_miniprot_detail"]:
        for mgc_eval_detail in outputs["ei_homology.mgc_evaluation_miniprot_detail"]:
            symlink(score_alignments_path, mgc_eval_detail)
    # File        mgc_score_summary = ScoreSummary.summary_table
    if outputs["ei_homology.mgc_score_summary"]:
        symlink(score_alignments_path, outputs["ei_homology.mgc_score_summary"])

    # Main output folder
    # Array[File] xspecies_combined_alignments = CombineXspecies.xspecies_scored_alignment
    # xspecies spaln
    xspecies_alignments_spaln = os.path.join(
        outputs_path, "ei_homology.xspecies_alignments", "spaln"
    )
    if not os.path.exists(xspecies_alignments_spaln):
        os.makedirs(xspecies_alignments_spaln)
    if outputs["ei_homology.xspecies_combined_alignments_spaln"]:
        for xspc_combined_aln in outputs[
            "ei_homology.xspecies_combined_alignments_spaln"
        ]:
            symlink(xspecies_alignments_spaln, xspc_combined_aln)
    # Array[File] alignment_filter_stats = AlignProteins.stats
    if outputs["ei_homology.spaln_alignment_filter_stats"]:
        for aln_filter_stat in outputs["ei_homology.spaln_alignment_filter_stats"]:
            symlink(xspecies_alignments_spaln, aln_filter_stat)

    # xspecies miniprot
    xspecies_alignments_miniprot = os.path.join(
        outputs_path, "ei_homology.xspecies_alignments", "miniprot"
    )
    if not os.path.exists(xspecies_alignments_miniprot):
        os.makedirs(xspecies_alignments_miniprot)
    if outputs["ei_homology.xspecies_combined_alignments_miniprot"]:
        for xspc_combined_aln in outputs[
            "ei_homology.xspecies_combined_alignments_miniprot"
        ]:
            symlink(xspecies_alignments_miniprot, xspc_combined_aln)
    if outputs["ei_homology.miniprot_alignment_filter_stats"]:
        for aln_filter_stat in outputs["ei_homology.miniprot_alignment_filter_stats"]:
            symlink(xspecies_alignments_miniprot, aln_filter_stat)

    # File loci = MikadoPick.loci
    symlink(outputs_path, outputs["ei_homology.loci"])
    # File scores = MikadoPick.scores
    symlink(outputs_path, outputs["ei_homology.scores"])
    # File metrics = MikadoPick.metrics
    symlink(outputs_path, outputs["ei_homology.metrics"])
    # File stats = MikadoPick.stats
    symlink(outputs_path, outputs["ei_homology.stats"])
    # File summary_stats = MikadoSummaryStats.summary
    symlink(outputs_path, outputs["ei_homology.stats_summary"])


def homology_module(cli_arguments):
    cromwell_inputs = combine_arguments_homology(cli_arguments)
    validate_homology_inputs(cromwell_inputs)

    cromwell_jar, runtime_config = prepare_cromwell_arguments(cli_arguments)

    with open(cli_arguments.output_parameters_file, "w") as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.homology_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        rc = execute_cromwell(
            runtime_config,
            cromwell_jar,
            cli_arguments.output_parameters_file,
            workflow_options_file,
            wdl_file,
        )
        if rc == 0:
            collect_homology_output(RUN_METADATA)

        return rc
