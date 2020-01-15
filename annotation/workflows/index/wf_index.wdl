version 1.0
workflow wf_index {
    input {
        File reference
    }

    call GSnapIndex {
        input: reference = reference
    }

    call hisat2Index {
        input: reference = reference
    }

    call starIndex {
        input: reference = reference
    }

    output {
        Array[File] gsnap_index = GSnapIndex.index

        Array[File] hisat_index = hisat2Index.index

        Array[File] star_index = starIndex.index
    }
}

task GSnapIndex {
    input {
        File reference
    }

    output {
        Array[File] index = glob("gsnapIndex/ref/*")
    }

    command <<<
    mkdir gsnapIndex
    gmap_build --dir="gsnapIndex" --db=ref ~{reference}
    >>>
}

task hisat2Index {
    input {
        File reference
    }

    output {
        Array[File] index = glob("ref*")
    }

    command <<<
        hisat2-build ~{reference} "ref"
    >>>
}

task starIndex {
    input {
        File reference
    }

    output {
        Array[File] index = glob('starIndex/*')
    }

    command <<<
        mkdir starIndex
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir starIndex \
        --genomeFastaFiles ~{reference}
    >>>
}

task tophatIndex {
    input {
        File reference
    }

    output {
        Array[File] index = glob('ref*')
    }

    command <<<
        bowtie2-build --threads 4 ~{reference} "ref"
    >>>
}