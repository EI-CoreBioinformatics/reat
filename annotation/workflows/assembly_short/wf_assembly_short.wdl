workflow wf_assembly_short {
    Array[Pair[File,File]] bams
    File? annotation

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
        Array[Array[File]] assemblies = [Stringtie.assembled, Scallop.assembled]
    }
}

task Stringtie {
    File bam
    File? annotation
    String strand = "--fr" # this needs to be computed from input paramaters
    command {
        stringtie ${bam} \
        -p 4 \
        ${strand} \
        ${"-G " + annotation} \
        -o assembled.gtf
    }

    output {
        File assembled = "assembled.gtf"
    }
}


# Needs to have the tool available... Not built yet for OSX
task Scallop {
    File bam
    String strand = "--library_type first"

    command {
        scallop --verbose 0 -i ${bam} -o "scallop.gtf" ${strand}
    }

    output {
        File assembled = "scallop.gtf"
    }
}