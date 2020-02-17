version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_align_short {
    input {
        Array[PRSample] samples
        File? annotation
        Array[File] gsnap_index
        Array[File] hisat_index
        Array[File] star_index
    }

    if (defined(annotation)) {
        call hisat2SpliceSites {
            input: annotation = annotation
        }
    }

    call Hisat {
        input:
        sites = hisat2SpliceSites.sites,
        sample = samples[0],
        index = hisat_index
    }

    call Star {
        input:
        annotation = annotation,
        sample = samples[0],
        index = star_index
    }

    Array[AlignedSample] aligned_samples = [Hisat.aligned_sample, Star.aligned_sample]

    scatter (aligned_sample in aligned_samples) {
        call Sort {
            input:
            sample = aligned_sample
        }
    }

    scatter (aligned_sample in aligned_samples) {
        call Stats {
            input:
            sample = aligned_sample
        }
    }

    output {
        Array[IndexedAlignedSample] indexed_aligned_samples = Sort.indexed_aligned_sample
        Array[Pair[File,File]] indexed_bams = Sort.indexed_bam
        Array[File] stats = Stats.stats
        Array[Array[File]] plots = Stats.plots
    }
}

task Sort{
    input {
        AlignedSample sample
        String name = basename(sample.bam, ".bam")
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        Pair[File,File] indexed_bam = (name + ".sorted.bam", name + ".sorted.bam.bai")
        IndexedAlignedSample indexed_aligned_sample = { "name": sample.name, 
                                                        "strand": sample.strand, 
                                                        "aligner": sample.aligner, 
                                                        "bam": sample.name+"."+sample.aligner+".sorted.bam", 
                                                        "index": sample.name+"."+sample.aligner+".sorted.bam.bai"
                                                    }
    }

    command <<<
        set -euxo pipefail
        samtools sort ~{sample.bam} > ~{sample.name + "." + sample.aligner + ".sorted.bam"}
        samtools index ~{sample.name + "." + sample.aligner + ".sorted.bam"}
    >>>
}

task Stats {
    input {
        AlignedSample sample
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File stats = sample.name + "." + sample.aligner + ".stats"
        Array[File] plots = glob("plot/*.png")
    }

    command <<<
        set -euxo pipefail
        samtools stats ~{sample.bam} > ~{sample.name + "." + sample.aligner + ".stats"} && \
        plot-bamstats -p "plot/~{sample.name + "." + sample.aligner}" ~{sample.name + "." + sample.aligner + ".stats"}
    >>>
}

task GSnapSpliceSites {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File sites = "gsnapSplicesites.iit"
    }

    command <<<
        set -euxo pipefail
        gtf_splicesites ~{annotation} | iit_store -o gsnapSplicesites.iit
    >>>
}

task hisat2SpliceSites {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File sites = "hisat2Splicesites.txt"
    }

    command <<<
        set -euxo pipefail
        hisat2_extract_splice_sites.py ~{annotation} > "hisat2Splicesites.txt"
    >>>
}

task GSnap {
    input {
        Array[File] index
        File? sites
        PRSample sample
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        AlignedSample aligned_sample = {"name": sample.name, "strand": sample.strand, "aligner": "hisat", "bam": sample.name+".gsnap.bam"}
    }

    command <<<
        set -euxo pipefail
        r1_file=~{sample.R1}
        r1_ext=${r1_file##*.}
        compression=""
        case "${r1_ext}" in
            gz)
            compression="${compression}--gunzip"
            ;;
            bz | bz2)
            compression="${compression}--bunzip2"
            ;;
        esac
        gsnap --nthreads ~{runtime_attr.cpu_cores} --dir="$(dirname "~{index[0]}")" \
        --db=ref \
        --novelsplicing=1 \
        ${compression} \
        ~{"-s " + sites} \
        --localsplicedist=2000 \
        --format=sam --npaths=20 \
        ~{sample.R1} ~{sample.R2} | samtools sort -@ 4 - > "~{sample.name}.gsnap.bam"
    >>>

}

task Hisat {
    input {
        Array[File] index
        File? sites
        PRSample sample
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        AlignedSample aligned_sample = {"name": sample.name, "strand": sample.strand, "aligner": "hisat", "bam": sample.name+".hisat.bam"}
    }

    command <<<
        set -euxo pipefail
        case "~{sample.strand}" in
            fr-firststrand)
            strandness="--rna-strandness=RF"
            ;;
            fr-secondstrand)
            strandness="--rna-strandness=FR"
            ;;
            f)
            strandness="--rna-strandness=F"
            ;;
            r)
            strandness="--rna-strandness=R"
            ;;
        esac
    hisat2 -p 4 -x ~{sub(index[0], "\\.\\d\\.ht2l?", "")} \
    ${strandness} \
    --min-intronlen=20 \
    --max-intronlen=2000 \
    ~{"--known-splicesite-infile " + sites} \
    -1 ~{sample.R1} ~{"-2 " + sample.R2} | samtools sort -@ 4 - > "~{sample.name}.hisat.bam"
    >>>
}

task Star {
    input {
        Array[File] index
        File? annotation
        PRSample sample
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File bam = "Aligned.out.bam"
        AlignedSample aligned_sample = {"name": sample.name, "strand": sample.strand, "aligner": "star", "bam": sample.name+".star.bam"}
    }

    command <<<
        set -euxo pipefail
        r1_file=~{sample.R1}
        r1_ext=${r1_file##*.}
        compression=""
        case "${r1_ext}" in
            gz)
            compression="--readFilesCommand \"gzip -dc\""
            ;;
            bz | bz2)
            compression="--readFilesCommand \"bzip2 -dc\""
            ;;
        esac

            STAR --genomeDir "$(dirname ~{index[0]})" \
    --runThreadN 4 \
    "${compression}" \
    --runMode alignReads \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM XS NM MD --outSAMstrandField intronMotif \
    --alignIntronMin 20 \
    --alignIntronMax 2000 \
    --alignMatesGapMax 2000 \
    ~{"--sjdbGTFfile " + annotation} \
    --readFilesIn ~{sample.R1} ~{sample.R2} && ln -s Aligned.out.bam "~{sample.name}.star.bam"
    >>>
}

task Tophat {
    input {
        Array[File] index
        PRSample sample
        File? annotation
        String strand
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File bam = "tophat2_accepted.bam"
        AlignedSample aligned_sample = {"name": sample.name, "strand": sample.strand, "aligner": "tophat2", "bam": sample.name+".tophat2.bam"}

    }

    command <<<
        set -euxo pipefail
        tophat2 \
        --num-threads 4 \
        --library-type=~{strand} \
        --min-intron-length=20 \
        --max-intron-length=2000 \
        ~{"--GTF " + annotation} \
        ~{sub(index[0], "\\.\\d\\.bt2l?", "")} \
        ~{sample.R1} ~{sample.R2} && ln -s tophat_out/accepted_hits.bam "~{sample.name}.tophat2.bam"
    >>>
}
