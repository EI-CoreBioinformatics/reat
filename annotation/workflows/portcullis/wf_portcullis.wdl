version 1.0

import "../structs/structs.wdl"

workflow portcullis {
    input {
        File reference
        Array[IndexedAlignedSample] aligned_samples
        File? annotation
    }

    if (defined(annotation)) {
        File def_annotation = select_first([annotation])
        call PrepareRef {
            input:
            annotation = def_annotation
        }
    }
    scatter (aligned_sample in aligned_samples) {
        call Prepare {
            input:
            reference = reference,
            sample = aligned_sample
        }
        call Junction {
            input:
            prep_dir = Prepare.prep_dir
        }
        call Filter {
            input:
            reference_bed = PrepareRef.refbed,
            prep_dir = Prepare.prep_dir,
            junc_dir = Junction.junc_dir,
            tab = Junction.tab
        }
    }
    
    call Merge {
        input:
        tabs = Filter.pass
    }

    output {
        File tab = Merge.tab
        File bed = Merge.bed
        File gff3 = Merge.gff3
    }
}

task PrepareRef {
    input {
        File annotation
    }

    output {
        File refbed = "reference.refbed"
    }

    command <<<
        junctools convert -if gtf -of ebed -o "reference.refbed" ~{annotation}
    >>>

}

task Prepare {
    input {
        File? reference
        IndexedAlignedSample sample
    }

    output {
        Array[File] prep_dir = glob("portcullis_prep/*")
    }

    command <<<
        portcullis prep -c -o portcullis_prep -t 4 ~{reference} ~{sample.bam}
    >>>

}

task Junction {
    input {
        Array[File] prep_dir
        String strand = "firststrand"
    }

    output {
        Array[File] junc_dir = glob("portcullis_junc/*")
        File tab = "portcullis_junc/portcullis.junctions.tab"
    }

    command <<<
        prep_dir_path="$(dirname ~{prep_dir[0]})"
        portcullis junc -c ~{"--strandedness="+strand} -t 4  "${prep_dir_path}"
    >>>
}

task Filter {
    input {
        Array[File] prep_dir
        Array[File] junc_dir
        File? reference_bed
        File tab
    }

    output {
        File pass = "portcullis_filter.pass.junctions.tab"
    }

    command <<<
        # junc_dir_path="$(dirname ~{junc_dir[0]})"
        prep_dir_path="$(dirname ~{prep_dir[0]})"

        portcullis filter -o portcullis_filter --canonical=OFF \
        --max_length=2000 ~{"--reference " + reference_bed } \
        --threads=4 "${prep_dir_path}" ~{tab}
    >>>
}

task Merge {
    input {
    Array[File] tabs

    }

    output {
        File tab = "portcullis.merged.tab"
        File bed = "portcullis.merged.bed"
        File gff3 = "portcullis.merged.gff3"
    }

    command <<<
        (junctools set --prefix=portcullis_merged --output=portcullis.merged.tab --operator=mean union ~{sep=" " tabs} || touch portcullis.merged.tab)
        junctools convert -if portcullis -of ebed --output=portcullis.merged.bed portcullis.merged.tab
        junctools convert -if portcullis -of igff --output=portcullis.merged.gff3 portcullis.merged.tab
    >>>
}