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

    with pkg_resources.path("annotation.prediction_module", "extrinsic.ei_augustus_generic.cfg") as extrinsic_path:
        cromwell_inputs['ei_prediction.extrinsic_config'] = str(extrinsic_path)

    if cli_arguments.expression:
        cromwell_inputs['ei_prediction.expressed_exon_hints'] = {'bam': cli_arguments.expression.name}
        if os.path.isfile(cli_arguments.expression.name + '.csi'):
            cromwell_inputs['ei_prediction.expressed_exon_hints']['index'] = cli_arguments.expression.name + '.csi'

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
