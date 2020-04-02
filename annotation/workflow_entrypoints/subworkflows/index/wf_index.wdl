version 1.0

import "../common/rt_struct.wdl"

workflow wf_index {
    input {
        File reference
    }

    call GSnapIndex {
        input: reference = reference
    }

    call Hisat2Index {
        input: reference = reference
    }

    call StarIndex {
        input: reference = reference
    }

    output {
        Array[File] gsnap_index = GSnapIndex.index

        Array[File] hisat_index = Hisat2Index.index

        Array[File] star_index = StarIndex.index
    }
}

task GSnapIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Array[File] index = glob("gsnapIndex/ref/*")
    }

    command <<<
        set -euxo pipefail
    mkdir gsnapIndex
    gmap_build --dir="gsnapIndex" --db=ref ~{reference}
    >>>
}

task Hisat2Index {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Array[File] index = glob("ref*")
    }

    command <<<
        set -euxo pipefail
        hisat2-build ~{reference} "ref"
    >>>
}

task StarIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    
    output {
        Array[File] index = glob('starIndex/*')
    }

    command <<<
        set -euxo pipefail
        mkdir StarIndex
        STAR --runThreadN ~{cpus} --runMode genomeGenerate --genomeDir starIndex \
        --genomeFastaFiles ~{reference}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task tophatIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    
    output {
        Array[File] index = glob('ref*')
    }

    command <<<
        set -euxo pipefail
        bowtie2-build --threads ~{cpus} ~{reference} "ref"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}