import os.path
from importlib import resources as pkg_resources
import json

from annotation import RUN_METADATA, prepare_cromwell_arguments, execute_cromwell


def combine_arguments_prediction(cli_arguments):
    computational_resources = {}
    if cli_arguments.computational_resources:
        computational_resources = json.load(cli_arguments.computational_resources)
    cromwell_inputs = computational_resources

    cromwell_inputs['ei_prediction.reference_genome'] = {'fasta': cli_arguments.genome.name}
    if os.path.isfile(cli_arguments.genome.name + '.fai'):
        cromwell_inputs['ei_prediction.reference_genome']['index'] = cli_arguments.genome.name + '.fai'
    cromwell_inputs['ei_prediction.augustus_config_path'] = cli_arguments.augustus_config_path
    cromwell_inputs['ei_prediction.species'] = cli_arguments.species
    cromwell_inputs['ei_prediction.kfold'] = cli_arguments.kfold
    cromwell_inputs['ei_prediction.chunk_size'] = cli_arguments.chunk_size
    cromwell_inputs['ei_prediction.overlap_size'] = cli_arguments.overlap_size

    if cli_arguments.firststrand_expression:
        cromwell_inputs['ei_prediction.firststrand_exon_hints'] = {'bam': cli_arguments.firststrand_expression.name}
        if os.path.isfile(cli_arguments.firststrand_expression.name + '.csi'):
            cromwell_inputs['ei_prediction.firststrand_exon_hints']['index'] = cli_arguments.firststrand_expression.name + '.csi'

    if cli_arguments.secondstrand_expression:
        cromwell_inputs['ei_prediction.secondstrand_exon_hints'] = {'bam': cli_arguments.secondstrand_expression.name}
        if os.path.isfile(cli_arguments.secondstrand_expression.name + '.csi'):
            cromwell_inputs['ei_prediction.secondstrand_exon_hints']['index'] = cli_arguments.secondstrand_expression.name + '.csi'

    if cli_arguments.unstranded_expression:
        cromwell_inputs['ei_prediction.unstranded_exon_hints'] = {'bam': cli_arguments.unstranded_expression.name}
        if os.path.isfile(cli_arguments.unstranded_expression.name + '.csi'):
            cromwell_inputs['ei_prediction.unstranded_exon_hints']['index'] = cli_arguments.unstranded_expression.name + '.csi'

    if cli_arguments.EVM_weights:
        cromwell_inputs['ei_prediction.EVM_weights'] = cli_arguments.EVM_weights.name

    with pkg_resources.path("annotation.prediction_module", "extrinsic.ei_augustus_generic.cfg") as extrinsic_path:
        cromwell_inputs['ei_prediction.extrinsic_config'] = str(extrinsic_path)

    # NOTE: The default extrinsic file provided with REAT can be overriden by the CLI
    if cli_arguments.extrinsic_config:
        cromwell_inputs['ei_prediction.extrinsic_config'] = cli_arguments.extrinsic_config.name

    # TODO: Parse the extrinsic_config and validate the augustus_runs SOURCEs against those available in the extrinsic
    if cli_arguments.augustus_runs:
        # TODO: check all files exist and are valid!
        cromwell_inputs['ei_prediction.augustus_runs'] = [f.name for f in cli_arguments.augustus_runs]

    if cli_arguments.introns:
        cromwell_inputs['ei_prediction.intron_hints'] = cli_arguments.introns.name

    if cli_arguments.repeats:
        cromwell_inputs['ei_prediction.repeats_gff'] = cli_arguments.repeats.name

    if cli_arguments.force_train:
        cromwell_inputs['ei_prediction.force_train'] = cli_arguments.force_train

    if cli_arguments.optimise_augustus:
        cromwell_inputs['ei_prediction.optimise_augustus'] = cli_arguments.optimise_augustus

    if cli_arguments.transcriptome_models:
        cromwell_inputs['ei_prediction.transcriptome_models'] = [f.name for f in cli_arguments.transcriptome_models]

    if cli_arguments.homology_models:
        cromwell_inputs['ei_prediction.homology_models'] = [f.name for f in cli_arguments.homology_models]

    if cli_arguments.homology_proteins:
        cromwell_inputs['ei_prediction.protein_validation_database'] = cli_arguments.homology_proteins.name

    if cli_arguments.hq_protein_alignments:
        cromwell_inputs['ei_prediction.HQ_protein_alignments'] = [f.name for f in cli_arguments.hq_protein_alignments]

    if cli_arguments.lq_protein_alignments:
        cromwell_inputs['ei_prediction.LQ_protein_alignments'] = [f.name for f in cli_arguments.lq_protein_alignments]

    if cli_arguments.hq_assembly:
        cromwell_inputs['ei_prediction.HQ_assembly'] = [f.name for f in cli_arguments.hq_assembly]

    if cli_arguments.lq_assembly:
        cromwell_inputs['ei_prediction.LQ_assembly'] = [f.name for f in cli_arguments.lq_assembly]

    if cli_arguments.mikado_utr_files:
        cromwell_inputs['ei_prediction.mikado_utr_files'] = ' '.join(cli_arguments.mikado_utr_files)

    if cli_arguments.do_glimmer:
        cromwell_inputs['ei_prediction.do_glimmer'] = 'true'
        if cli_arguments.do_glimmer is not True and os.access(cli_arguments.do_glimmer, os.R_OK):
            cromwell_inputs['ei_prediction.glimmer_training'] = cli_arguments.do_glimmer

    if cli_arguments.do_snap:
        cromwell_inputs['ei_prediction.do_snap'] = 'true'
        if cli_arguments.do_snap is not True and os.access(cli_arguments.do_snap, os.R_OK):
            cromwell_inputs['ei_prediction.snap_training'] = cli_arguments.do_snap

    if cli_arguments.do_codingquarry:
        cromwell_inputs['ei_prediction.do_codingquarry'] = 'true'
        if cli_arguments.do_codingquarry is not True and os.access(cli_arguments.do_codingquarry, os.R_OK):
            cromwell_inputs['ei_prediction.codingquarry_training'] = cli_arguments.do_codingquarry

    if cli_arguments.no_augustus and cli_arguments.no_augustus is False:
        cromwell_inputs['ei_prediction.do_augustus'] = 'false'

    if cli_arguments.min_pct_cds_fraction:
        cromwell_inputs['ei_prediction.LengthChecker.min_pct_cds_fraction'] = cli_arguments.min_pct_cds_fraction

    if cli_arguments.max_tp_utr_complete:
        cromwell_inputs['ei_prediction.LengthChecker.max_tp_utr_complete'] = cli_arguments.max_tp_utr_complete

    if cli_arguments.max_tp_utr:
        cromwell_inputs['ei_prediction.LengthChecker.max_tp_utr'] = cli_arguments.max_tp_utr

    if cli_arguments.min_tp_utr:
        cromwell_inputs['ei_prediction.LengthChecker.min_tp_utr'] = cli_arguments.min_tp_utr

    if cli_arguments.max_fp_utr_complete:
        cromwell_inputs['ei_prediction.LengthChecker.max_fp_utr_complete'] = cli_arguments.max_fp_utr_complete

    if cli_arguments.max_fp_utr:
        cromwell_inputs['ei_prediction.LengthChecker.max_fp_utr'] = cli_arguments.max_fp_utr

    if cli_arguments.min_fp_utr:
        cromwell_inputs['ei_prediction.LengthChecker.min_fp_utr'] = cli_arguments.min_fp_utr
    ###
    if cli_arguments.query_start_hard_filter_distance:
        cromwell_inputs['ei_prediction.LengthChecker.query_start_hard_filter_distance'] = cli_arguments.query_start_hard_filter_distance

    if cli_arguments.query_start_score:
        cromwell_inputs['ei_prediction.LengthChecker.query_start_score'] = cli_arguments.query_start_score

    if cli_arguments.query_start_scoring_distance:
        cromwell_inputs['ei_prediction.LengthChecker.query_start_scoring_distance'] = cli_arguments.query_start_scoring_distance
    ####
    if cli_arguments.target_start_hard_filter_distance:
        cromwell_inputs['ei_prediction.LengthChecker.target_start_hard_filter_distance'] = cli_arguments.target_start_hard_filter_distance

    if cli_arguments.target_start_score:
        cromwell_inputs['ei_prediction.LengthChecker.target_start_score'] = cli_arguments.target_start_score

    if cli_arguments.target_start_scoring_distance:
        cromwell_inputs['ei_prediction.LengthChecker.target_start_scoring_distance'] = cli_arguments.target_start_scoring_distance
    ####
    if cli_arguments.query_end_hard_filter_distance:
        cromwell_inputs['ei_prediction.LengthChecker.query_end_hard_filter_distance'] = cli_arguments.query_end_hard_filter_distance

    if cli_arguments.query_end_score:
        cromwell_inputs['ei_prediction.LengthChecker.query_end_score'] = cli_arguments.query_end_score

    if cli_arguments.query_end_scoring_distance:
        cromwell_inputs['ei_prediction.LengthChecker.query_end_scoring_distance'] = cli_arguments.query_end_scoring_distance
    ####
    if cli_arguments.target_end_hard_filter_distance:
        cromwell_inputs['ei_prediction.LengthChecker.target_end_hard_filter_distance'] = cli_arguments.target_end_hard_filter_distance

    if cli_arguments.target_end_score:
        cromwell_inputs['ei_prediction.LengthChecker.target_end_score'] = cli_arguments.target_end_score

    if cli_arguments.target_end_scoring_distance:
        cromwell_inputs['ei_prediction.LengthChecker.target_end_scoring_distance'] = cli_arguments.target_end_scoring_distance
    ###
    if cli_arguments.min_query_coverage_hard_filter:
        cromwell_inputs['ei_prediction.LengthChecker.min_query_coverage_hard_filter'] = cli_arguments.min_query_coverage_hard_filter

    if cli_arguments.min_query_coverage_score:
        cromwell_inputs['ei_prediction.LengthChecker.min_query_coverage_score'] = cli_arguments.min_query_coverage_score

    if cli_arguments.min_query_coverage_scoring_percentage:
        cromwell_inputs['ei_prediction.LengthChecker.min_query_coverage_scoring_percentage'] = cli_arguments.min_query_coverage_scoring_percentage
    ###
    if cli_arguments.min_target_coverage_hard_filter:
        cromwell_inputs['ei_prediction.LengthChecker.min_target_coverage_hard_filter'] = cli_arguments.min_target_coverage_hard_filter

    if cli_arguments.min_target_coverage_score:
        cromwell_inputs['ei_prediction.LengthChecker.min_target_coverage_score'] = cli_arguments.min_target_coverage_score

    if cli_arguments.min_target_coverage_scoring_percentage:
        cromwell_inputs['ei_prediction.LengthChecker.min_target_coverage_scoring_percentage'] = cli_arguments.min_target_coverage_scoring_percentage
    ###
    if cli_arguments.max_single_gap_hard_filter:
        cromwell_inputs['ei_prediction.LengthChecker.max_single_gap_hard_filter'] = cli_arguments.max_single_gap_hard_filter
    if cli_arguments.max_single_gap_score:
        cromwell_inputs['ei_prediction.LengthChecker.max_single_gap_score'] = cli_arguments.max_single_gap_score
    if cli_arguments.max_single_gap_scoring_length:
        cromwell_inputs['ei_prediction.LengthChecker.max_single_gap_scoring_length'] = cli_arguments.max_single_gap_scoring_length

    if cli_arguments.filter_top_n:
        cromwell_inputs['ei_prediction.SelfBlastFilter.top_n'] = cli_arguments.filter_top_n
    if cli_arguments.filter_max_coverage:
        cromwell_inputs['ei_prediction.SelfBlastFilter.coverage'] = cli_arguments.filter_max_coverage
    if cli_arguments.filter_max_identity:
        cromwell_inputs['ei_prediction.SelfBlastFilter.identity'] = cli_arguments.filter_max_identity

    return cromwell_inputs


def collect_prediction_output(run_metadata):
    ...


def prediction_module(cli_arguments):
    cromwell_inputs = combine_arguments_prediction(cli_arguments)

    cromwell_jar, runtime_config = prepare_cromwell_arguments(cli_arguments)

    with open(cli_arguments.output_parameters_file, 'w') as cromwell_input_file:
        json.dump(cromwell_inputs, cromwell_input_file)
    # Submit pipeline to server or run locally depending on the arguments
    with pkg_resources.path("annotation.prediction_module", "main.wdl") as wdl_file:
        workflow_options_file = None
        if cli_arguments.workflow_options_file is not None:
            workflow_options_file = cli_arguments.workflow_options_file.name
        rc = execute_cromwell(runtime_config, cromwell_jar,
                              cli_arguments.output_parameters_file, workflow_options_file, wdl_file)
        if rc == 0:
            collect_prediction_output(RUN_METADATA)
        return rc
