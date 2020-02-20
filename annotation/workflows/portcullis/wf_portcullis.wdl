version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow portcullis {
    input {
        File reference
        Array[AlignedSample] aligned_samples
        File? annotation
    }

    if (defined(annotation)) {
        File def_annotation = select_first([annotation])
        call PrepareRef {
            input:
            annotation = def_annotation
        }
    }
    scatter (aligned_sample in aligned_samples) {
        call Prepare {
            input:
            reference = reference,
            sample = aligned_sample
        }
        call Junction {
            input:
            prep_dir = Prepare.prep_dir
        }
        call Filter {
            input:
            reference_bed = PrepareRef.refbed,
            prep_dir = Prepare.prep_dir,
            junc_dir = Junction.junc_dir,
            tab = Junction.tab
        }
    }
    
    call Merge {
        input:
        tabs = Filter.pass
    }

    output {
        File tab = Merge.tab
        File bed = Merge.bed
        File gff3 = Merge.gff3
    }
}

task PrepareRef {
    input {
        File annotation
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
        File refbed = "reference.refbed"
    }

    command <<<
        set -euxo pipefail
        junctools convert -if gtf -of ebed -o "reference.refbed" ~{annotation}
    >>>

}

task Prepare {
    input {
        File? reference
        AlignedSample sample
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        Array[File] prep_dir = glob("portcullis_prep/*")
    }

    command <<<
        set -euxo pipefail
        portcullis prep -c -o portcullis_prep -t ~{cpus} ~{reference} ~{sep=" "sample.bam}
    >>>

   RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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

task Junction {
    input {
        Array[File] prep_dir
        String strand = "firststrand"
        RuntimeAttr? runtime_attr_override
    }
    
    Int cpus = 8

    output {
        Array[File] junc_dir = glob("portcullis_junc/*")
        File tab = "portcullis_junc/portcullis.junctions.tab"
    }

    command <<<
        set -euxo pipefail
        prep_dir_path="$(dirname ~{prep_dir[0]})"
        portcullis junc -c ~{"--strandedness="+strand} -t ~{cpus}  "${prep_dir_path}"
    >>>
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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

task Filter {
    input {
        Array[File] prep_dir
        Array[File] junc_dir
        File? reference_bed
        File tab
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        File pass = "portcullis_filter.pass.junctions.tab"
    }

    command <<<
        set -euxo pipefail
        # junc_dir_path="$(dirname ~{junc_dir[0]})"
        prep_dir_path="$(dirname ~{prep_dir[0]})"

        portcullis filter -o portcullis_filter --canonical=OFF \
        --max_length=2000 ~{"--reference " + reference_bed } \
        --threads=~{cpus} "${prep_dir_path}" ~{tab}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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

task Merge {
    input {
        Array[File] tabs
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
        File tab = "portcullis.merged.tab"
        File bed = "portcullis.merged.bed"
        File gff3 = "portcullis.merged.gff3"
    }

    command <<<
        set -euxo pipefail
        (junctools set --prefix=portcullis_merged --output=portcullis.merged.tab --operator=mean union ~{sep=" " tabs} || touch portcullis.merged.tab)
        junctools convert -if portcullis -of ebed --output=portcullis.merged.bed portcullis.merged.tab
        junctools convert -if portcullis -of igff --output=portcullis.merged.gff3 portcullis.merged.tab
    >>>
}