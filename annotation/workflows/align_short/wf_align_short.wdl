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
        String aligner = "hisat"
    }
    
    if (defined(annotation)) {
        call hisat2SpliceSites {
            input: annotation = annotation
        }
    }

    if (aligner == "hisat") {
        scatter (sample in samples) {
            scatter(PR in sample.read_pair) {
                call Hisat {
                    input:
                    sites = hisat2SpliceSites.sites,
                    strand = sample.strand,
                    name = sample.name,
                    sample = PR,
                    index = hisat_index
                }
            }
            AlignedSample hisat_aligned_sample = object {"name": sample.name, "strand": sample.strand, "aligner": "hisat", "bam": Hisat.aligned_pair}
        }
    }

    if (aligner == "star") {
        scatter (sample in samples) {
            scatter(PR in sample.read_pair) {
                call Star {
                    input:
                    annotation = annotation,
                    strand = sample.strand,
                    name = sample.name,
                    sample = PR,
                    index = star_index
                }
            }
            AlignedSample star_aligned_sample = object {"name": sample.name, "strand": sample.strand, "aligner": "star", "bam": Star.aligned_pair}
        }
    }

    Array[AlignedSample] def_aligned_samples = select_first([hisat_aligned_sample, star_aligned_sample])

    scatter (aligned_sample in def_aligned_samples) {
        scatter (bam in aligned_sample.bam) {
            call Sort {
                input:
                bam = bam
            }
        }
        # IndexedAlignedSample indexed_aligned_sample = object{"name": aligned_sample.name, "strand": aligned_sample.strand, "aligner": aligned_sample.aligner, "index_bam": Sort.indexed_bam}
        AlignedSample sorted_aligned_sample = object{"name": aligned_sample.name, "strand": aligned_sample.strand, "aligner": aligned_sample.aligner, "bam": Sort.sorted_bam}
    }

    scatter (aligned_sample in def_aligned_samples) {
        scatter (bam in aligned_sample.bam) {
            call Stats {
                input:
                bam = bam
            }
        }
    }

    output {
        # Array[IndexedAlignedSample] indexed_aligned_samples = indexed_aligned_sample
        Array[AlignedSample] aligned_samples = sorted_aligned_sample
        Array[Array[File]] stats = Stats.stats
        Array[Array[Array[File]]] plots = Stats.plots
    }
}

task Sort{
    input {
        File bam
        String name = basename(bam, ".bam")
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        IndexedBam indexed_bam = { "bam": name + ".sorted.bam", "index": name + ".sorted.bam.bai" }
        File sorted_bam = name + ".sorted.bam"
    }

    command <<<
        set -euxo pipefail
        samtools sort -@~{cpus} ~{bam} > ~{name + ".sorted.bam"}
        samtools index ~{name + ".sorted.bam"}
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task Stats {
    input {
        File bam
        String name = basename(bam, ".bam")
        RuntimeAttr? runtime_attr_override
    }

    output {
        File stats = name + ".stats"
        Array[File] plots = glob("plot/*.png")
    }

    command <<<
        set -euxo pipefail
        samtools stats ~{bam} > ~{name + ".stats"} && \
        plot-bamstats -p "plot/~{name}" ~{name + ".stats"}
    >>>
    
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
}

task GSnapSpliceSites {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File sites = "gsnapSplicesites.iit"
    }

    command <<<
        set -euxo pipefail
        gtf_splicesites ~{annotation} | iit_store -o gsnapSplicesites.iit
    >>>

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
}

task hisat2SpliceSites {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File sites = "hisat2Splicesites.txt"
    }

    command <<<
        set -euxo pipefail
        hisat2_extract_splice_sites.py ~{annotation} > "hisat2Splicesites.txt"
    >>>

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
}

task Hisat {
    input {
        Array[File] index
        File? sites
        ReadPair sample
        String strand
        String name
        String rp_name = name+"."+basename(sample.R1)
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        File aligned_pair = rp_name + ".hisat.bam"
    }

    command <<<
        set -euxo pipefail
        case "~{strand}" in
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
    hisat2 -p ~{cpus} -x ~{sub(index[0], "\\.\\d\\.ht2l?", "")} \
    ${strandness} \
    --min-intronlen=20 \
    --max-intronlen=2000 \
    ~{"--known-splicesite-infile " + sites} \
    -1 ~{sample.R1} ~{"-2 " + sample.R2} | samtools sort -@ 4 - > "~{rp_name}.hisat.bam"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task Star {
    input {
        Array[File] index
        File? annotation
        ReadPair sample
        String strand
        String name
        String rp_name = name+"."+basename(sample.R1)
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        File aligned_pair = rp_name + ".star.bam"
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
    --runThreadN ~{cpus} \
    "${compression}" \
    --runMode alignReads \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM XS NM MD --outSAMstrandField intronMotif \
    --alignIntronMin 20 \
    --alignIntronMax 2000 \
    --alignMatesGapMax 2000 \
    ~{"--sjdbGTFfile " + annotation} \
    --readFilesIn ~{sample.R1} ~{sample.R2} && ln -s Aligned.out.bam "~{rp_name}.star.bam"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
