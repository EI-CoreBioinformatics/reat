workflow portcullis {
    File? annotation
    Array[Pair[File,File]] bams

    call PrepareRef {
        input:
        annotation = annotation
    }

    scatter (bam in bams) {
        call Prepare {
            input:
            reference = PrepareRef.refbed,
            bam = bam.left
        }
        call Junction {
            input:
            prep_dir = Prepare.prep_dir
        }
        call Filter {
            input:
            prep_dir = Prepare.prep_dir,
            junc_dir = Junction.junc_dir,
            tab = Junction.tab
        }
    }
    
    # Add the call to merge after knowing how the output of the previous steps looks like 
    # call Merge {
    #     input:
    #     tabs = Filter.filter_dir
    # }

    output {
        Array[Array[File]] tabs = Filter.filter_dir
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
    File reference
    File bam

    command {
        portcullis prep -c -o portcullis_prep -t 4 ${reference} ${bam}
    }

    output {
        Array[File] prep_dir = glob("portcullis_prep/*")
    }

}

task Junction {
    Array[File] prep_dir
    String dollar = "$"
    String strand = "fr-firststrand"

    command <<<
        prep_dir_path=dirname ${prep_dir[0]}

        portcullis junc -c -o portcullis_junc ${strand} -t 4  ${dollar}{prep_dir_path}
    >>>

    output {
        Array[File] junc_dir = glob("portcullis_junc/*")
        File tab = "portcullus_junc/portcullis.junctions.tab"
    }
}

task Filter {
    Array[File] prep_dir
    Array[File] junc_dir
    File? reference_bed
    String dollar = "$"
    File tab

    command <<<
        junc_dir_path=dirname ${junc_dir[0]}
        prep_dir_path=dirname ${prep_dir[0]}

        portcullis filter -o portcullis_filter --canonical=OFF \
        --max_length=2000 ${"--reference " + reference_bed } \
        --threads=4 ${dollar}{prep_dir_path} ${tab}
    >>>

    output {
        Array[File] filter_dir = glob("portcullis_filter/*")
    }

}

# task Merge {
#     Array[Array[File]] tabs

#     command {
#         (junctools set --prefix=portcullis_merged --output=portcullis.merged.tab --operator=mean union ${sep=" " tabs} || touch portcullis.merged.tab)
#         junctools convert -if portcullis -of ebed --output=portcullis.merged.bed portcullis.merged.tab
#         junctools convert -if portcullis -of igff --output=portcullis.merged.gff3 portcullis.merged.tab
#     }

#     output {
#         File tab = "portcullis.merged.tab"
#         File bed = "portcullis.merged.bed"
#         File gff3 = "portcullis.merged.gff3"
#     }
# }