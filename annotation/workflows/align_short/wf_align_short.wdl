workflow wf_align_short {
    File R1
    File? R2
    File? annotation
    Array[File] gsnap_index
    Array[File] hisat_index
    Array[File] star_index

    if (defined(annotation)) {
        call hisat2SpliceSites {
            input: annotation = annotation
        }
    }

    call Hisat {
        input:
        sites = hisat2SpliceSites.sites,
        R1 = R1,
        R2 = R2,
        index = hisat_index
    }

    call Star {
        input:
        annotation = annotation,
        R1 = R1,
        R2 = R2,
        index = star_index
    }

    # if (defined(annotation)) {
    #     call GSnapSpliceSites {
    #         input:
    #         annotation = annotation
    #     }
    # }

    # call GSnap {
    #     input:
    #     sites = GSnapSpliceSites.sites,
    #     R1 = R1,
    #     R2 = R2,
    #     index = gsnap_index
    # }

    # call Tophat {
    #     input:
    #     index = tophat_index,
    #     annotation = annotation,
    #     R1 = R1,
    #     R2 = R2,
    #     strand = "fr-firststrand"
    # }

    # Array[File] bams = [Hisat.bam, Star.bam, GSnap.bam, Tophat.bam]

    Array[File] bams = [Hisat.bam, Star.bam]

    scatter (bam in bams) {
        call Sort {
            input:
            bam = bam
        }
    }

    scatter (bam in bams) {
        call Stats {
            input:
            bam = bam
        }
    }

    output {
        Array[Pair[File,File]] indexed_bams = Sort.indexed_bam
        Array[File] stats = Stats.stats
        Array[Array[File]] plots = Stats.plots
    }
}

task Sort{
    File bam
    String name = basename(bam, ".bam")
    command {
        samtools sort ${bam} > ${name + ".sorted.bam"}
        samtools index ${name + ".sorted.bam"}
    }
    output {
        Pair[File,File] indexed_bam = (name + ".sorted.bam", name + ".sorted.bam.bai")
    }
}

task Stats {
    File bam
    String name = basename(bam, ".bam")
    
    command {
        samtools stats ${bam} > ${name + ".stats"} && \
        plot-bamstats -p plot/ ${name + ".stats"}
    }

    output {
        File stats = name + ".stats"
        Array[File] plots = glob("plot/*.png")
    }
}

task GSnapSpliceSites {
    File? annotation
    command {
        gtf_splicesites ${annotation} | iit_store -o gsnapSplicesites.iit
    }
    output {
        File sites = "gsnapSplicesites.iit"
    }
}

task hisat2SpliceSites {
    File? annotation

    command {
        hisat2_extract_splice_sites.py ${annotation} > "hisat2Splicesites.txt"
    }

    output {
        File sites = "hisat2Splicesites.txt"
    }
}

task GSnap {
    Array[File] index
    File? sites
    File R1
    File? R2
    String dollar = "$"

    command <<<
        r1_file=${R1}
        r1_ext=${dollar}{r1_file##*.}
        compression=""
        case "${dollar}{r1_ext}" in
            gz)
            compression="${dollar}{compression}--gunzip"
            ;;
            bz | bz2)
            compression="${dollar}{compression}--bunzip2"
            ;;
        esac
        gsnap --nthreads 4 --dir="$(dirname "${index[0]}")" \
        --db=ref \
        --novelsplicing=1 \
        ${dollar}{compression} \
        ${"-s " + sites} \
        --localsplicedist=2000 \
        --format=sam --npaths=20 \
        ${R1} ${R2} | samtools sort -@ 4 - > "gsnap.bam"
    >>>

    output {
        File bam = "gsnap.bam"
    }
}

task Hisat {
    Array[File] index
    File? sites
    File R1
    File? R2
    String strand = "fr-firststrand"
    String dollar = "$"

    command <<<
        strandness=""
        case "${strand}" in
            fr-firststrand)
            strandness="--rna-strandness=RF"
            ;;
            fr-secondstrand)
            strandness = "--rna-strandness=FR"
            ;;
            f)
            strandness = "--rna-strandness=F"
            ;;
            r)
            strandness = "--rna-strandness=R"
            ;;
        esac
    hisat2 -p 4 -x ${sub(index[0], "\\.\\d\\.ht2l?", "")} \
    ${dollar}{strandness} \
    --min-intronlen=20 \
    --max-intronlen=2000 \
    ${"--known-splicesite-infile " + sites} \
    -1 ${R1} ${"-2 " + R2} | samtools sort -@ 4 - > "hisat.bam"
    >>>

    output {
        File bam = "hisat.bam"
    }
}

task Star {
    Array[File] index
    File? annotation
    File R1
    File? R2
    String dollar = "$"
    
    command <<<
        r1_file=${R1}
        r1_ext=${dollar}{r1_file##*.}
        compression=""
        case "${dollar}{r1_ext}" in
            gz)
            compression="${dollar}{compression}--readFilesCommand \"gzip -dc\""
            ;;
            bz | bz2)
            compression="${dollar}{compression}--readFilesCommand \"bzip2 -dc\""
            ;;
        esac

            STAR --genomeDir "$(dirname ${index[0]})" \
    --runThreadN 4 \
    ${dollar}{compression} \
    --runMode alignReads \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM XS NM MD --outSAMstrandField intronMotif \
    --alignIntronMin 20 \
    --alignIntronMax 2000 \
    --alignMatesGapMax 2000 \
    ${"--sjdbGTFfile " + annotation} \
    --readFilesIn ${R1} ${R2} && ln -s Aligned.out.bam star.bam
    >>>

    output {
        File test_ok = "Aligned.out.bam"
        File bam = "star.bam"
    }
}

task Tophat {
    Array[File] index
    File R1
    File? R2
    File? annotation
    String strand
    command {
        tophat2 \
        --num-threads 4 \
        --library-type=${strand} \
        --min-intron-length=20 \
        --max-intron-length=2000 \
        ${"--GTF " + annotation} \
        ${sub(index[0], "\\.\\d\\.bt2l?", "")} \
        ${R1} ${R2} && ln -s tophat_out/accepted_hits.bam tophat2_accepted.bam
    }

    output {
        File bam = "tophat2_accepted.bam"
    }
}
