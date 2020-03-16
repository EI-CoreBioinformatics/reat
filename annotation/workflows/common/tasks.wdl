version 1.0

import "./structs.wdl"
import "./rt_struct.wdl"

task SplitSequences {
    input {
        File sequences_file
        String prefix = "out"
        Int num_out_files = 10
    }

    output {
        Array[File] seq_files = glob(prefix+"*")
    }

    command <<<
        set -euxo pipefail
        seqtk split -n ~{num_out_files} ~{prefix} ~{sequences_file}
    >>>
}

task MergeFiles {
    input {
        Array[File] files_to_merge
        String output_filename
    }

    output {
        File out = output_filename
    }

    command <<<
    set -euxo pipefail
    cat ~{sep=" " files_to_merge} > ~{output_filename}
    >>>
}

task IndexFasta {
    input {
        File reference_fasta
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        IndexedReference indexed_fasta = {"fasta": basename(reference_fasta), "fai": basename(reference_fasta)+".fai"}
    }

    command <<<
        set -euxo pipefail
        ln -s ~{reference_fasta} .
        samtools faidx ~{basename(reference_fasta)}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task sanitizeFasta {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File sanitised_reference = "reference.san.fasta"
    }

    command <<<
        set -euxo pipefail
        sanitize_sequence_db.py -o "reference.san.fasta" ~{reference}
    >>>
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
