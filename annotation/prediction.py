import os.path
import shutil
from importlib import resources as pkg_resources
import json

from annotation import RUN_METADATA, prepare_cromwell_arguments, execute_cromwell
from annotation.utils import symlink


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

    if cli_arguments.min_train_models:
        cromwell_inputs['ei_prediction.SelectAugustusTestAndTrain.min_train_models'] = cli_arguments.min_train_models

    if cli_arguments.max_train_models:
        cromwell_inputs['ei_prediction.SelectAugustusTestAndTrain.max_train_models'] = cli_arguments.max_train_models

    if cli_arguments.max_test_models:
        cromwell_inputs['ei_prediction.SelectAugustusTestAndTrain.max_test_models'] = cli_arguments.max_test_models

    if cli_arguments.target_mono_exonic_percentage:
        cromwell_inputs['ei_prediction.SelectAugustusTestAndTrain.target_mono_exonic_percentage'] = cli_arguments.target_mono_exonic_percentage

    if cli_arguments.force_train_few_models:
        cromwell_inputs["ei_prediction.SelectAugustusTestAndTrain.force"] = 'true'

    if cli_arguments.evalue_filter:
        cromwell_inputs['ei_prediction.LengthChecker.evalue_filter'] = cli_arguments.evalue_filter

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

    if cli_arguments.codingquarry_extra_params:
        cromwell_inputs['ei_prediction.codingquarry_extra_params'] = cli_arguments.codingquarry_extra_params
    if cli_arguments.glimmer_extra_params:
        cromwell_inputs['ei_prediction.glimmer_extra_params'] = cli_arguments.glimmer_extra_params
    if cli_arguments.snap_extra_params:
        cromwell_inputs['ei_prediction.snap_extra_params'] = cli_arguments.snap_extra_params
    if cli_arguments.augustus_extra_params:
        cromwell_inputs['ei_prediction.augustus_extra_params'] = cli_arguments.augustus_extra_params
    if cli_arguments.evm_extra_params:
        cromwell_inputs['ei_prediction.evm_extra_params'] = cli_arguments.evm_extra_params

    if cli_arguments.mikado_config:
        cromwell_inputs['ei_prediction.mikado_config'] = cli_arguments.mikado_config

    if cli_arguments.mikado_scoring:
        cromwell_inputs['ei_prediction.mikado_scoring'] = cli_arguments.mikado_scoring

    return cromwell_inputs


def collect_prediction_output(run_metadata):
    run_metadata = json.load(open(run_metadata))
    outputs = run_metadata['outputs']
    outputs_path = 'outputs'
    if not os.path.exists(outputs_path):
        os.mkdir(outputs_path)

    if outputs['ei_prediction.augustus_config']:
        symlink(outputs_path, outputs['ei_prediction.augustus_config'], "final_augustus_config")

    if outputs['ei_prediction.glimmer']:
        symlink(outputs_path, outputs['ei_prediction.glimmer'])

    if outputs['ei_prediction.snap']:
        symlink(outputs_path, outputs['ei_prediction.snap'])
    if outputs['ei_prediction.codingquarry']:
        symlink(outputs_path, outputs['ei_prediction.codingquarry'])
    if outputs['ei_prediction.codingquarry_fresh']:
        symlink(outputs_path, outputs['ei_prediction.codingquarry_fresh'])
    if outputs['ei_prediction.augustus_abinitio']:
        symlink(outputs_path, outputs['ei_prediction.augustus_abinitio'])

    if outputs['ei_prediction.augustus']:
        for o in outputs['ei_prediction.augustus']:
            symlink(outputs_path, o)

    if outputs['ei_prediction.evm_predictions']:
        symlink(outputs_path, outputs['ei_prediction.evm_predictions'])
    symlink(outputs_path, outputs['ei_prediction.mikado_loci'])
    symlink(outputs_path, outputs['ei_prediction.mikado_stats'])
    symlink(outputs_path, outputs['ei_prediction.mikado_summary_stats'])

    predictions_path = os.path.join(outputs_path, 'GenePredictors')
    glimmer_prediction_path = os.path.join(predictions_path, "GlimmerHMM")
    snap_prediction_path = os.path.join(predictions_path, "SNAP")
    codingquarry_prediction_path = os.path.join(predictions_path, "CodingQuarry")
    augustus_prediction_path = os.path.join(predictions_path, "Augustus")
    if not os.path.exists(predictions_path):
        os.mkdir(predictions_path)

    if outputs['ei_prediction.predictions_codingquarry']:
        if not os.path.exists(codingquarry_prediction_path):
            os.mkdir(codingquarry_prediction_path)
        symlink(codingquarry_prediction_path, outputs['ei_prediction.predictions_codingquarry'])

    if outputs['ei_prediction.predictions_codingquarry_fresh']:
        if not os.path.exists(codingquarry_prediction_path):
            os.mkdir(codingquarry_prediction_path)
        symlink(codingquarry_prediction_path, outputs['ei_prediction.predictions_codingquarry_fresh'])

    if not (outputs['ei_prediction.predictions_codingquarry'] or outputs['ei_prediction.predictions_codingquarry_fresh']):
        if os.path.exists(codingquarry_prediction_path):
            shutil.rmtree(codingquarry_prediction_path)

    if outputs['ei_prediction.predictions_snap']:
        if not os.path.exists(snap_prediction_path):
            os.mkdir(snap_prediction_path)
        symlink(snap_prediction_path, outputs['ei_prediction.predictions_snap'])
    else:
        if os.path.exists(snap_prediction_path):
            shutil.rmtree(snap_prediction_path)

    if outputs['ei_prediction.predictions_glimmer']:
        if not os.path.exists(glimmer_prediction_path):
            os.mkdir(glimmer_prediction_path)
        symlink(glimmer_prediction_path, outputs['ei_prediction.predictions_glimmer'])
    else:
        if os.path.exists(glimmer_prediction_path):
            shutil.rmtree(glimmer_prediction_path)

    if outputs['ei_prediction.predictions_augustus']:
        if not os.path.exists(augustus_prediction_path):
            os.mkdir(augustus_prediction_path)
        for i, o in enumerate(outputs['ei_prediction.predictions_augustus']):
            symlink(augustus_prediction_path, o, f"augustus_run{i+1}.gff")
    if outputs['ei_prediction.predictions_augustus_abinitio']:
        symlink(augustus_prediction_path, outputs['ei_prediction.predictions_augustus_abinitio'], "augustus_abinitio.gff")

    if not (outputs['ei_prediction.predictions_augustus'] or outputs['ei_prediction.predictions_augustus_abinitio']):
        if os.path.exists(augustus_prediction_path):
            shutil.rmtree(augustus_prediction_path)

    training_path = os.path.join(outputs_path, 'GenePredictorsTraining')
    glimmer_training_path = os.path.join(training_path, "GlimmerHMM")
    snap_training_path = os.path.join(training_path, "SNAP")
    augustus_training_path = os.path.join(training_path, "Augustus")
    if not os.path.exists(training_path):
        os.mkdir(training_path)
        os.mkdir(glimmer_training_path)
        os.mkdir(snap_training_path)
        os.mkdir(augustus_training_path)

    if outputs['ei_prediction.training_selected_models']:
        symlink(training_path, outputs['ei_prediction.training_selected_models'], "training_models.gff")
    if outputs['ei_prediction.training_augustus_etraining_evaluation']:
        if not os.path.exists(augustus_training_path):
            os.mkdir(augustus_training_path)
        symlink(augustus_training_path, outputs['ei_prediction.training_augustus_etraining_evaluation'], "base_training_evaluation.txt")
        symlink(augustus_training_path, outputs['ei_prediction.training_augustus_etraining_training'], "base_augustus_training")
    if outputs['ei_prediction.training_augustus_optimise_augustus_evaluation']:
        symlink(augustus_training_path, outputs['ei_prediction.training_augustus_optimise_augustus_evaluation'], "optimise_training_evaluation.txt")
        symlink(augustus_training_path, outputs['ei_prediction.training_augustus_optimise_augustus_training'], "optimised_augustus_training")
    if not (outputs['ei_prediction.training_augustus_etraining_evaluation'] or outputs['ei_prediction.training_augustus_optimise_augustus_evaluation']):
        if os.path.exists(augustus_training_path):
            shutil.rmtree(augustus_training_path)

    if outputs['ei_prediction.training_glimmer_training']:
        if not os.path.exists(glimmer_training_path):
            os.mkdir(glimmer_training_path)
        symlink(glimmer_training_path, outputs['ei_prediction.training_glimmer_training'])
    else:
        if os.path.exists(glimmer_training_path):
            shutil.rmtree(glimmer_training_path)

    if outputs['ei_prediction.training_snap_training']:
        if not os.path.exists(snap_training_path):
            os.mkdir(snap_training_path)
        symlink(snap_training_path, outputs['ei_prediction.training_snap_training'])
    else:
        if os.path.exists(snap_training_path):
            shutil.rmtree(snap_training_path)

    classification_path = os.path.join(outputs_path, 'Classification')
    if not os.path.exists(classification_path):
        os.mkdir(classification_path)

    symlink(classification_path, outputs['ei_prediction.classification_gold_models'])
    symlink(classification_path, outputs['ei_prediction.classification_silver_models'])
    symlink(classification_path, outputs['ei_prediction.classification_bronze_models'])
    symlink(classification_path, outputs['ei_prediction.classification_all'])
    symlink(classification_path, outputs['ei_prediction.classification_non_redundant'])


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
