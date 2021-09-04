version 1.0
import "../../common/rt_struct.wdl"

workflow wf_prodigal {
    input {
        File prepared_transcripts
        Int gencode
        String output_directory
        RuntimeAttr? prodigal_runtime_attr
    }

    call Prodigal {
        input:
        prepared_transcripts = prepared_transcripts,
        gencode_id = gencode,
        output_directory = output_directory,
        runtime_attr_override = prodigal_runtime_attr
    }

    output {
        File orfs = Prodigal.orfs
    }
}

task Prodigal {
    input {
        File prepared_transcripts
        Int gencode_id

        String output_directory
        RuntimeAttr? runtime_attr_override
    }

    output {
        File orfs = output_directory+"/transcripts.fasta.prodigal.gff3"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        prodigal -f gff -g "~{gencode_id}" -i "~{prepared_transcripts}" -o "transcripts.fasta.prodigal.gff3"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}
