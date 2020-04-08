version 1.0

import "../../common/rt_struct.wdl"

task SanitiseProteinBlastDB {
    input {
        File db
        RuntimeAttr? runtime_attr_override
    }

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

    output {
        File clean_db = "output.db"
    }

    command <<<
        set -euxo pipefail
        sanitize_sequence_db.py -cstop ~{db} | gt seqtransform -addstopaminos -width "60" > "output.db"

    >>>
}

task BlastIndex {
    input {
        File target
        RuntimeAttr? runtime_attr_override
    }
    
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

    output {
        Array[File] index = glob("blast_index.db.*")
    }

    command <<<
        set -euxo pipefail
        makeblastdb -dbtype prot -in ~{target} -out blast_index.db
    >>>
}

task BlastAlign {
    input {
        Array[File] index
        File query
        String? extra
        String blast_type
        String outfmt
        String output_filename
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    runtime {
        cpu: cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File out = output_filename
    }

    command <<<
        set -euxo pipefail
        ~{blast_type} ~{extra} -db ~{sub(index[0], "\.[^.]+$", "")} -num_threads ~{cpus} -query ~{query} -outfmt ~{outfmt} > ~{output_filename}
    >>>
}

task DiamondIndex {
    input {
        File target
        RuntimeAttr? runtime_attr_override
    }
    
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

    output {
        File index = "diamond_index.dmnd"
    }

    command <<<
        set -euxo pipefail
        diamond makedb --in ~{target} --db "diamond_index.dmnd"
    >>>
}

task DiamondAlign {
    input {
        File index
        File query
        String? extra
        String blast_type
        String output_filename
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    runtime {
        cpu: cpus
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File out = output_filename
    }

    command <<<
        set -euxo pipefail
        diamond ~{blast_type} ~{extra} -p "~{cpus}" -d "~{index}" -q "~{query}" > "~{output_filename}"
    >>>
}
