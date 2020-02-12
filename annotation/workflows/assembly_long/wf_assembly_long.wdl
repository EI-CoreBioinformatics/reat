version 1.0

import "../common/structs.wdl"

workflow wf_assembly_long {
    input {
        File reference
        Array[AlignedSample] bams
        String assembler = "None"
    }

    scatter (bam in bams) {
        if (assembler == "None") {
            call sam2gff {
                input:
                bams = [bam]
            }
        }
        if (assembler == "merge") {
            call gffread_merge {
                input:
                bams = [bam]
            }
        }

        if (assembler == "stringtie") {
            call stringtie_long {
                input:
                reference = reference,
                bams = [bam]
            }
        }
        File def_gff = select_first([sam2gff.gff, gffread_merge.gff, stringtie_long.gff])
    }

    output {
        Array[File] gff = def_gff
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