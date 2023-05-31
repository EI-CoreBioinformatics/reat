import io
import os
import signal
import subprocess
import sys

from Mikado.parsers import parser_factory
from Mikado.transcripts import Gene, Transcript

VERSION = '0.7.1'
RUN_METADATA = "run_details.json"

UTR_SELECTION_OPTIONS = ('augustus', 'gold', 'silver', 'bronze', 'all', 'hq_assembly', 'lq_assembly')
LONG_READ_ALIGNER_CHOICES = ('minimap2', 'gmap', '2pass', '2pass_merged')

def report_errors(errors, samples):
    if any([len(error_list) for error_list in errors.values()]):
        print(f"File {samples.name} parsing failed, errors found:\n", file=sys.stderr)
        for line, error_list in errors.items():
            if not error_list:
                continue
            print("Line:", line.strip(), sep='\n\t', file=sys.stderr)
            print("was not parsed successfully, the following errors were found:", file=sys.stderr)
            [print("\t-", e, file=sys.stderr) for e in error_list]
        raise ValueError(f"Could not parse file {samples.name}")


def prepare_cromwell_arguments(cli_arguments):
    cromwell_jar = os.environ.get('CROMWELL_JAR', None)
    if cli_arguments.jar_cromwell:
        cromwell_jar = cli_arguments.jar_cromwell.name
    runtime_config = os.environ.get('CROMWELL_RUNTIME_CONFIG', None)
    if cli_arguments.runtime_configuration:
        runtime_config = cli_arguments.runtime_configuration.name
    return cromwell_jar, runtime_config


def execute_cromwell(workflow_configuration_file, jar_cromwell, input_parameters_filepath, workflow_options_file,
                     wdl_file):
    return cromwell_run(workflow_configuration_file, jar_cromwell,
                        input_parameters_filepath, workflow_options_file, wdl_file)


def cromwell_run(workflow_configuration_file, jar_cromwell, input_parameters_filepath, workflow_options_file, wdl_file,
                 log_level="INFO"):
    formatted_command_line = ["java", f"-Dconfig.file={str(workflow_configuration_file)}",
                              f"-DLOG_LEVEL={log_level}",
                              "-jar", str(jar_cromwell),
                              "run",
                              "-i", str(input_parameters_filepath)]
    if workflow_options_file:
        formatted_command_line.extend(["-o", str(workflow_options_file)])

    formatted_command_line.extend(["-m", RUN_METADATA, str(wdl_file)])

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
                print("\n\n\nREAT Failed, the following file might contain information with the reasons behind"
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


def kill_cromwell(sig, frame):
    raise KeyboardInterrupt


class HashableTranscript(Transcript):
    def __repr__(self):
        return self.id

    def __hash__(self):
        return hash(self.id)


def minimal_gxf_parser(file):
    parser = parser_factory(file)
    genes = dict()
    tid2gid = dict()
    for row in parser:
        if row.header is True:
            continue
        elif row.is_gene is True:
            genes[row.id] = Gene(row)
        elif row.is_transcript is True:
            assert len(row.parent) == 1
            parent = row.parent[0]
            tid2gid[row.id] = parent
            genes[parent].add(HashableTranscript(row))
        elif row.is_exon is True:
            if row.gene is None:
                gene = tid2gid[row.parent[0]]
            else:
                gene = row.gene
            genes[gene].add_exon(row)

    for gene in genes:
        genes[gene].finalize()
    return genes, tid2gid
