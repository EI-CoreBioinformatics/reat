version 1.0

import "../common/structs.wdl"

workflow wf_assembly_long {
    input {
        File reference
        Array[AlignedSample] bams
        String assembler = "None"
    }

    if (assembler == "None") {
        call sam2gff {
            input:
            bams = bams
        }
    }
    if (assembler == "merge") {
        call gffread_merge {
            input:
            bams = bams
        }
    }

    if (assembler == "stringtie") {
        call stringtie_long {
            input:
            reference = reference,
            bams = bams
        }
    }

    File def_gff = select_first([sam2gff.gff, gffread_merge.gff, stringtie_long.gff]) 
    output {
        File gtf = def_gff
    }
}

task stringtie_long {
    input {
        File reference
        Array[AlignedSample] bams
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    stringtie ~{"-G " + reference} -L ~{sep=" " bams} -o "result.gff"
    >>>
}

task sam2gff {
    input {
        Array[AlignedSample] bams
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    sam2gff ~{sep=" " bams} > result.gff
    >>>
}

task gffread_merge {
    input {
        Array[AlignedSample] bams
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    sam2gff ~{sep=" " bams} | gffread -T -M -K -o result.gff
    >>>
}