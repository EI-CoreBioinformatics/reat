version 1.0

import "../common/structs.wdl"
import "../common/tasks.wdl" as tsk

workflow wf_sanitize {
    input {
        File reference_genome
        File? in_annotation
    }

    call tsk.sanitizeFasta {
        input:
            reference = reference_genome
    }

    if (defined(in_annotation)) {
        call sanitizeAnnotation {
            input:
            annotation = in_annotation
        }
        File wf_maybe_clean_annotation = sanitizeAnnotation.sanitised_annotation
    }

    call tsk.IndexFasta {
        input:
        reference_fasta = sanitizeFasta.sanitised_reference
    }
    
    output {
        File? annotation = wf_maybe_clean_annotation
        File reference = sanitizeFasta.sanitised_reference
        IndexedReference indexed_reference = IndexFasta.indexed_fasta
    }
}

task sanitizeAnnotation {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File sanitised_annotation = "reference.san.gtf"
    }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    command <<<
        set -euxo pipefail
        filepath=~{annotation}
        if [ ${filepath##*.} == "gff" ]
        then
            mikado util convert -of gtf ~{annotation} "reference.san.gtf"
        else
            ln ~{annotation} "reference.san.gtf"
        fi
    >>>
}