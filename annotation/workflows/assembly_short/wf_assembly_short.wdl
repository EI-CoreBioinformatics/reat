version 1.0
workflow wf_assembly_short {
    input {
        Array[Pair[File,File]] bams
        File? annotation
    }

    scatter (bam in bams) {
        call Stringtie{
            input:
            bam = bam.left,
            annotation = annotation
        }

        call Scallop {
            input:
            bam = bam.left
        }
    }

    output {
        Array[File] assemblies = flatten([Stringtie.assembled, Scallop.assembled])
    }
}

task Stringtie {
    input {
        File bam
        File? annotation
        String strand = "--fr" # this needs to be computed from input paramaters
    }

    output {
        File assembled = "assembled.gtf"
    }

    command <<<
        stringtie ~{bam} \
        -p 4 \
        ~{strand} \
        ~{"-G " + annotation} \
        -o assembled.gtf
    >>>
}


# Needs to have the tool available... Not built yet for OSX
task Scallop {
    input {
        File bam
        String strand = "--library_type first"
    }

    output {
        File assembled = "scallop.gtf"
    }

    command <<<
        scallop --verbose 0 -i ~{bam} -o "scallop.gtf" ~{strand}
    >>>
}