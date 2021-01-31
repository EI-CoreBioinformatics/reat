version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"
workflow wf_align_long {
    input {
        File reference
        Array[LRSample] long_samples
        Boolean is_hq
        File? bed_junctions
        String aligner = "minimap2"
        String? aligner_extra_parameters
        Float min_identity = 0.9
        Int? min_intron_len = 20
        Int? max_intron_len = 200000
        Int? max_intron_len_ends = 100000
        RuntimeAttr? alignment_resources
        RuntimeAttr? indexing_resources
    }

    parameter_meta {
        reference: {description:"Genome target to align the reads", category: "required"}
        long_samples: {description:"Long read samples, each item is defined by a name, it's strand and one or more long read files.", category: "required"}
        is_hq: {description:"Selects high quality parameters for the alignment program.", category: "required"}
        bed_junctions: {description:"Where possible uses a user provided set of junctions to guide the alignments.", category: "required"}
        aligner: {description:"Selects the aligner program, the options are: minimap2 and gmap.", category: "required"}
        alignment_resources: {description:"Computational resources to override the defaults for running the alignments.", category: "required"}
        indexing_resources: {description:"Computational resources to generate the genome target index previous to alignment, overrides defaults.", category: "required"}
    }
    
    # Add aligner option
    if (aligner == "minimap2") {
        scatter (sample in long_samples) {
            scatter (LR in sample.LR) {
                call Minimap2Long {
                    input:
                    LR = LR,
                    is_hq = is_hq,
                    strand = sample.strand,
                    name = sample.name,
                    reference = reference,
                    bed_junctions = bed_junctions,
                    max_intron_len = select_first([max_intron_len, 200000]),
                    extra_parameters = aligner_extra_parameters,
                    runtime_attr_override = alignment_resources
                }
            }
            AlignedSample mm2_aligned_sample = object {name: sample.name, strand:sample.strand, aligner:"minimap2", bam: Minimap2Long.bam, merge:false,
                                               is_ref: select_first([sample.is_ref, false]), exclude_redundant: select_first([sample.exclude_redundant, false]),
                                               score: select_first([sample.score, 0])}
        }
    }

    if (aligner == "gmap") {
        call GMapIndex {
            input:
            reference = reference,
            runtime_attr_override = indexing_resources
        }
        scatter (sample in long_samples) {
            scatter (LR in sample.LR) {
                call GMapLong {
                    input:
                    LR = LR,
                    gmap_index = GMapIndex.gmap_index,
                    strand = sample.strand,
                    name = sample.name,
                    reference = reference,
                    min_identity = min_identity,
                    min_intron_len = select_first([min_intron_len,20]),
                    max_intron_len = select_first([max_intron_len,2000]),
                    max_intron_len_ends = select_first([max_intron_len_ends, 4000]),
                    extra_parameters = aligner_extra_parameters,
                    runtime_attr_override = alignment_resources
                }
            }
            AlignedSample gmap_aligned_sample = object {name: sample.name, strand:sample.strand, aligner:"gmap", bam: GMapLong.bam, merge:false,
                                                score: select_first([sample.score, 1]),
                                                is_ref: select_first([sample.is_ref, false]),
                                                exclude_redundant: select_first([sample.exclude_redundant, false])}
        }
    }

    Array[AlignedSample] def_alignments = select_first([mm2_aligned_sample, gmap_aligned_sample])

    scatter (aligned_sample in def_alignments) {
        scatter (bam in aligned_sample.bam) {
            call LongAlignmentStats {
                input:
                bam = bam
            }
        }

        call SummaryAlignmentStats {
            input:
            stats = LongAlignmentStats.stats,
            output_prefix = aligned_sample.name
        }
        File summary_alignment_stats = SummaryAlignmentStats.summary_stats
    }

    call CollectAlignmentStats {
        input:
        summary_stats = summary_alignment_stats,
        prefix = if(is_hq) then "HQ" else "LQ"
    }
    output {
        Array[AlignedSample] bams = def_alignments
        Array[File] summary_stats = summary_alignment_stats
        File summary_stats_table = CollectAlignmentStats.summary_stats_table
    }

}

task CollectAlignmentStats {
    input {
        Array[File] summary_stats
        String prefix
    }

    output {
        File summary_stats_table = "align_long." + prefix +"_read_samples.summary.stats.tsv"
    }

    command <<<
    short_read_summary_stats_table ~{sep=" " summary_stats} > align_long.~{prefix}_read_samples.summary.stats.tsv
    >>>
}

task SummaryAlignmentStats {
    input {
        Array[File] stats
        String output_prefix
    }

    output {
        File summary_stats = "alignments/" + output_prefix + ".summary.stats.tsv"
    }

    command <<<
    mkdir alignments
    cd alignments
    alignment_summary_stats ~{sep=" " stats} > ~{output_prefix}.summary.stats.tsv
    >>>
}

task LongAlignmentStats {
    input {
        File bam
        String name = basename(bam, ".bam")
        RuntimeAttr? runtime_attr_override
    }

    output {
        File stats = "align_long_stats/" + name + ".stats"
        Array[File] plots = glob("align_long_stats/plot/*")
    }

    command <<<
        set -euxo pipefail
        mkdir align_long_stats
        cd align_long_stats
        samtools stats ~{bam} > ~{name + ".stats"} && \
        plot-bamstats -p "plot/~{name}" ~{name + ".stats"}
    >>>
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}


task GMapIndex {
    input {
        File reference
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

    output {
        Array[File] gmap_index = glob("gmapIndex/test_genome/*")
    }


    command <<<
        set -euxo pipefail
        gmap_build --dir=gmapIndex --db=test_genome ~{reference}
    >>>
}

task GMapExonsIIT {
    input {
        File? annotation
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

    output {
        File iit = "gmap_exons.iit"
    }

    command <<<
        set -euxo pipefail
        gtf_genes ~{annotation} | iit_store -o gmap_exons.iit
    >>>
}

task GMapLong {
    input {
        File reference
        Array[File] gmap_index
        File LR
        String LR_basename = sub(basename(LR), "\.[^/.]+$", "")
        String name
        String strand
        String? extra_parameters
        File? iit
        Int? min_trimmed_coverage
        Float min_identity
        Int max_intron_len_ends
        Int max_intron_len
        Int min_intron_len
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 16
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File bam = "alignments/gmap."+name+"."+LR_basename+".bam"
    }

    command <<<
        set -euxo pipefail
        filename=$(basename -- "~{LR}")
        extension="${filename##*.}"

        in_pipe="cat ~{LR}"
        if [ "$extension" == "bam" ]
        then
            in_pipe="samtools fastq ~{LR}"
        elif [ "$extension" == "gz" ]
        then
            in_pipe="gunzip -c ~{LR}"
        fi

        strand_opt=""
        if [ "~{strand}" == "fr-firststrand" ]
        then
            strand_opt="-z antisense_force"
        fi

        if [ "~{strand}" == "fr-secondstrand" ]
        then
            strand_opt="-z sense_force"
        fi

        mkdir alignments
        cd alignments

        $in_pipe | $(determine_gmap.py ~{reference}) --dir="$(dirname ~{gmap_index[0]})" --db=test_genome \
        ~{"--min-intronlength=" + min_intron_len} ~{"--max-intronlength-middle=" + max_intron_len} \
        ~{"--max-intronlength-ends=" + max_intron_len_ends} --npaths=1 \
        ~{"-m " + iit} ${strand_opt} \
        ~{"--min-trimmed-coverage=" + min_trimmed_coverage} \
        ~{"--min-identity=" + min_identity} \
        --format=samse ~{extra_parameters} \
        --nthreads="~{task_cpus}" | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ 4 --reference ~{reference} -T gmap.sort -o gmap.~{name}.~{LR_basename}.bam -
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }

}

task Minimap2Long {
    input {
        File LR
        String LR_basename = sub(basename(LR), "\.[^/.]+$", "")
        Boolean is_hq
        String name
        String strand
        File reference
        Int max_intron_len
        String? extra_parameters
        File? bed_junctions
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 16
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File bam = "alignments/minimap2."+name+"."+LR_basename+".bam"
    }

    command <<<
        set -euxo pipefail
        # Replace long_sample.LR with samtools fastq if suffix is bam
        filename=$(basename -- "~{LR}")
        extension="${filename##*.}"

        in_pipe="cat ~{LR}"
        if [ "$extension" == "bam" ]
        then
            in_pipe="samtools fastq ~{LR}"
        elif [ "$extension" == "gz" ]
        then
            in_pipe="gunzip -c ~{LR}"
        fi

        strand_opt="-ub"
        if [ "~{strand}" == "fr-secondstrand" ]
        then
            strand_opt="-uf"
        fi
        
        mkdir alignments
        cd alignments
        $in_pipe | \
        minimap2 ~{extra_parameters} \
        -ax ~{if (is_hq) then "splice:hq" else "splice"} \
        ~{"--junc-bed " + bed_junctions} \
        --cs=long \
        -G ~{max_intron_len} \
        ${strand_opt} \
        -t ~{task_cpus} \
        -L --MD \
        --eqx -2 \
        --secondary=no \
        ~{reference} - | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ ~{task_cpus/2} --reference ~{reference} -T minimap2.sort -o minimap2.~{name}.~{LR_basename}.bam -
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}