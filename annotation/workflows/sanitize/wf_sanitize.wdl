version 1.0
workflow wf_sanitize {
    input {
        File reference_genome
        File? in_annotation
    }

    call sanitizeReference {
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

    call indexReference {
        input:
        reference = sanitizeReference.sanitised_reference
    }
    
    output {
        File? annotation = wf_maybe_clean_annotation
        File reference = sanitizeReference.sanitised_reference
        File index = indexReference.sanitised_reference_index
    }
}

task sanitizeReference {
    input {
        File reference
    }

    output {
        File sanitised_reference = 'reference.san.fasta'
    }

    command <<<
        sanitize_sequence_db.py -o "reference.san.fasta" ~{reference}
    >>>
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

task indexReference {
    input {
        File reference
    }

    output {
        File sanitised_reference_index = "reference.san.fasta.fai"
    }

    command <<<
        ln -s ~{reference} "reference.san.fasta"
        samtools faidx "reference.san.fasta"
    >>>
}
