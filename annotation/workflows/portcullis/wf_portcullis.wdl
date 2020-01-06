workflow portcullis {
    Array[Pair[File,File]] bams
    File? annotation

    call PrepareRef {
        input:
        annotation = annotation
    }

    scatter (bam in bams) {
        call Prepare {
            input:
            bam = bam
        }
        call Junction {
            input:
            csi = Prepare.out
        }
        call Filter {
            input:
            tab = Junction.tab,
            bed = Junction.bed
        }
    }

}

task PrepareRef {
    File? annotation

    command {
        junctools convert -if gtf -of ebed -o "reference.refbed" ${annotation}
    }

    output {
        File refbed = "reference.refbed"
    }
}

task Prepare {

}

task Junction {

}

task Filter {

}

task Merge {

}