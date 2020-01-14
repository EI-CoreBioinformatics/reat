workflow wf_align_long {
    File? LR
    File? annotation
    Array[File] star_index
    File reference

    if(defined(LR)) {
        call GMapIndex {
            input:
            reference = reference
        }

        if (defined(annotation)) {
            call GMapExonsIIT {
                input:
                annotation = annotation
            }
        }

        call GMapLong {
            input:
            reference = reference,
            LR = LR,
            gmap_index = GMapIndex.gmap_index,
            iit = GMapExonsIIT.iit
        }
        # call StarLong {
        #     input:
        #     LR = LR,
        #     annotation = annotation,
        #     index = star_index
        # }

        call Minimap2Long {
            input:
            LR = LR,
            reference = reference
        }

        call Minimap2Long2Gff {
            input:
            bam = Minimap2Long.bam
        }
    }

    output {
        Array[File?] bams = [Minimap2Long.bam]
        Array[File?] gff = [GMapLong.gff, Minimap2Long2Gff.gff]
    }

}

task GMapIndex {
    File reference

    command {
        gmap_build --dir=gmapIndex --db=test_genome ${reference}
    }

    output {
        Array[File] gmap_index = glob("gmapIndex/test_genome/*")
    }
}

task GMapExonsIIT {
    File? annotation

    command {
        gtf_genes ${annotation} | iit_store -o gmap_exons.iit
    }

    output {
        File? iit = "gmap_exons.iit"
    }
}

task GMapLong {
    File reference
    Array[File] gmap_index
    File? LR
    File? iit
    Int? min_trimmed_coverage
    Int? min_identity
    String? strand

    command {
        gzcat ${LR} | $(determine_gmap.py ${reference}) --dir="$(dirname ${gmap_index[0]})" --db=test_genome \
        --min-intronlength=20 --intronlength=2000 \
        ${"-m " + iit} \
        ${"--min-trimmed-coverage=" + min_trimmed_coverage} \
        ${"--min-identity" + min_identity} \
        ${"-z " + strand} \
        --format=gff3_match_cdna \
        --nthreads=4 > gmap.out.gff
    }

    output {
        File gff = "gmap.out.gff"
    }
}

task StarLong {
    File? LR
    Array[File] index
    File? annotation
    String dollar = "$"
    
    command <<<
        lr_file=${LR}
        lr_ext=${dollar}{lr_file##*.}
        compression=""
        case "${dollar}{lr_ext}" in
            gz)
            compression="${dollar}{compression}--readFilesCommand \"gzip -dc\""
            ;;
            bz | bz2)
            compression="${dollar}{compression}--readFilesCommand \"bzip2 -dc\""
            ;;
        esac

            STARlong --genomeDir "$(dirname ${index[0]})" \
    --runThreadN 4 \
    ${dollar}{compression} \
    --runMode alignReads \
    --readNameSeparator space \
    --outFilterMultimapScoreRange 1 \
    --outFilterMismatchNmax 2000 \
    --scoreGapNoncan -20 \
    --scoreGapGCAG -4 \
    --scoreGapATAC -8 \
    --scoreDelOpen -1 \
    --scoreDelBase -1 \
    --scoreInsOpen -1 \
    --scoreInsBase -1 \
    --alignEndsType Local \
    --seedSearchStartLmax 50 \
    --seedPerReadNmax 100000 \
    --seedPerWindowNmax 1000 \
    --alignTranscriptsPerReadNmax 100000 \
    --alignTranscriptsPerWindowNmax 10000 \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI NM MD \
    --outSAMstrandField intronMotif \
    --alignIntronMin 20 \
    --alignIntronMax 2000 \
    --alignMatesGapMax 2000 \
    ${"--sjdbGTFfile " + annotation} \
    --outFileNamePrefix STARlong.out \
    --readFilesIn ${LR}
    >>>

    output {
        File bam = "STARlong.out.bam"
    }
}


task Minimap2Long {
    File? LR
    File reference

    command {
        minimap2 \
        -x splice \
        --cs=long \
        -G 2000 \
        -u b \
        -t 4 \
        -a -L --MD \
        --eqx -2 \
        --secondary=no \
        ${reference} ${LR} | samtools view -bS - | samtools sort -@ 4 --reference ${reference} -T minimap2.sort -o minimap2.out.bam -
    }

    output {
        File bam = "minimap2.out.bam"
    }
}

task Minimap2Long2Gff {
    File bam

    command {
        paftools.js splice2bed -m <(samtools view -h ${bam}) | correct_bed12_mappings.py > output.bed12
    }

    output {
        File gff = "output.bed12"
    }
}