version 1.0

import "../common/structs.wdl"

workflow wf_assembly_long {
    input {
        File reference
        Array[AlignedSample] aligned_samples
        String assembler = "None"
    }

    scatter (sample in aligned_samples) {
        if (assembler == "None") {
            call sam2gff {
                input:
                aligned_sample = sample
            }
        }
        if (assembler == "merge") {
            call gffread_merge {
                input:
                aligned_sample = sample
            }
        }

        if (assembler == "stringtie") {
            call stringtie_long {
                input:
                reference = reference,
                aligned_sample = sample
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
        AlignedSample aligned_sample
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    stringtie ~{"-G " + reference} -L ~{aligned_sample.bam} -o "result.gff"
    >>>
}

task sam2gff {
    input {
        AlignedSample aligned_sample
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    samtools view -F 4 -F 0x900 ~{aligned_sample.bam} | sam2gff -s ~{aligned_sample.name} > result.gff
    >>>
}

task gffread_merge {
    input {
        AlignedSample aligned_sample
    }

    output {
        File gff = "result.gff"
    }

    command <<<
    samtools view -F 4 -F 0x900 ~{aligned_sample.bam} | sam2gff -s ~{aligned_sample.name} | gffread -T -M -K -o result.gff
    >>>
}