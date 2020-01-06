workflow wf_index {
    File reference

    call GSnapIndex {
        input: reference = reference
    }

    call hisat2Index {
        input: reference = reference
    }

    call starIndex {
        input: reference = reference
    }

    call tophatIndex {
        input: reference = reference
    }

    output {
        Array[File] gsnap_index = GSnapIndex.index

        Array[File] hisat_index = hisat2Index.index

        Array[File] tophat_index = tophatIndex.index
        Array[File] star_index = starIndex.index
    }
}

task GSnapIndex {
    File reference

    command {
    mkdir gsnapIndex
    gmap_build --dir="gsnapIndex" --db=ref ${reference}
    }
    output {
        Array[File] index = glob("gsnapIndex/ref/*")
    }
}

task hisat2Index {
    File reference

    command {
        hisat2-build ${reference} "ref"
    }

    output {
        Array[File] index = glob("ref*")
    }
}

task starIndex {
    File reference

    command {
        mkdir starIndex
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir starIndex \
        --genomeFastaFiles ${reference}
    }

    output {
        Array[File] index = glob('starIndex/*')
    }
}

task tophatIndex {
    File reference

    command {
        bowtie2-build --threads 4 ${reference} "ref"
    }

    output {
        Array[File] index = glob('ref*')
    }
}