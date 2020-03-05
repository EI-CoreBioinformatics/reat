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
        File index = "blast_index.db"
    }

    command <<<
        set -euxo pipefail
        makeblastdb -in ~{target} > "blast_index.db"
    >>>
}

task BlastAlign {
    input {
        File index
        File query
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
        File out = "output.xml"
    }

    command <<<
        set -euxo pipefail
        blast ~{index} ~{query} > "output.txt"
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
        File index = "diamond_index.db"
    }

    command <<<
        set -euxo pipefail
        diamond index ~{target} > "diamond_index.db"
    >>>
}

task DiamondAlign {
    input {
        File index
        File query
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
        File out = "output.txt"
    }

    command <<<
        set -euxo pipefail
        diamond blast ~{index} ~{query} > "output.xml"
    >>>
}
