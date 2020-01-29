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
    }

    output {
        File sanitised_annotation = "reference.san.gtf"
    }

    command <<<
        filepath=~{annotation}
        if [ ${filepath##*.} = "gff" ]
        then
            mikado util convert -of gtf ~{annotation} "reference.san.gtf"
        else
            ln ~{annotation} "reference.san.gtf"
        fi
    >>>
}