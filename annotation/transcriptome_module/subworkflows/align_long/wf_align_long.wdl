version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"
import "wf_minimap.wdl" as mm2
import "wf_gmap.wdl" as gmap

workflow wf_align_long {
    input {
        IndexedReference indexed_reference
        File? annotation
        Array[LRSample] long_samples
        Boolean is_hq
	    Boolean skip_2pass_alignment = false
        File? extra_junctions
        File? portcullis_junctions
        String aligner = "minimap2"
        String? aligner_extra_parameters
        Float min_identity = 0.9
        Int? min_intron_len = 20
        Int? max_intron_len = 200000
        Int? max_intron_len_ends = 100000
        RuntimeAttr? twopass_resources
        RuntimeAttr? twopass_merge_resources
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

    if (defined(annotation)) {
        call gff2bed {
            input:
            annotation = select_first([annotation])
        }
    }

    if (defined(annotation) || defined(portcullis_junctions) || defined(extra_junctions)) {
        call JunctionUnion {
            input:
            junction_files = select_all([gff2bed.bed, portcullis_junctions, extra_junctions])
        }
    }

    if (aligner == "minimap2" || aligner == "2pass" || aligner == "2pass_merged") {
        call Minimap2Index {
            input:
            is_hq = is_hq,
            reference = indexed_reference.fasta,
            indexing_resources = indexing_resources
        }
        File mmi_2pass = select_first([Minimap2Index.index])
        scatter (sample in long_samples) {
            call mm2.wf_mm2 as minimap2 {
                input:
                reference = mmi_2pass,
                bed_junctions = JunctionUnion.combined_junctions,
                LRS = sample.LR,
                is_hq = is_hq,
                name = "base_alignment",
                strand = sample.strand,
                score = sample.score,
                is_ref = sample.is_ref,
                exclude_redundant = sample.exclude_redundant,
                max_intron_len = max_intron_len,
                aligner_extra_parameters = aligner_extra_parameters,
                alignment_resources = alignment_resources
            }
            if (aligner == "2pass" || aligner == "2pass_merged") {
                call twopass_score {
                    input:
                    strand = sample.strand,
                    name = sample.name,
                    alignment = minimap2.aligned_sample.bam,
                    reference = indexed_reference.fasta,
                    reference_fai = indexed_reference.fai,
                    runtime_attr_override = twopass_resources
                }
            }
            if (aligner == "2pass" && !skip_2pass_alignment) {
                File twopass_junctions = select_first([twopass_score.filtered_junctions])
                call mm2.wf_mm2 as minimap2_2pass {
                    input:
                    reference = mmi_2pass,
                    bed_junctions = twopass_junctions,
                    LRS = sample.LR,
                    is_hq = is_hq,
                    name = "2pass",
                    strand = sample.strand,
                    score = sample.score,
                    is_ref = sample.is_ref,
                    exclude_redundant = sample.exclude_redundant,
                    max_intron_len = max_intron_len,
                    aligner_extra_parameters = aligner_extra_parameters,
                    alignment_resources = alignment_resources
                }
            }
            AlignedSample def_minimap_1pass = select_first([minimap2_2pass.aligned_sample, minimap2.aligned_sample])
        }

		if (aligner == "2pass_merged") {
			# Merge the scored junctions from all samples
			Array[File] scored_junctions = select_all(select_first([twopass_score.scored_junctions]))
			call twopass_merge {
				input:
				reference = indexed_reference.fasta,
				reference_fai = indexed_reference.fai,
				scored_beds = scored_junctions,
				runtime_attr_override = twopass_merge_resources
			}
			# And align all samples with the new merged junctions
			if (!skip_2pass_alignment){
				scatter (sample in long_samples) {
					call mm2.wf_mm2 as minimap2_2pass_merged {
						input:
						reference = mmi_2pass,
						bed_junctions = twopass_merge.merged_filtered_junctions,
						LRS = sample.LR,
						is_hq = is_hq,
						name = "2pass_merged",
						strand = sample.strand,
						score = sample.score,
						is_ref = sample.is_ref,
						exclude_redundant = sample.exclude_redundant,
						max_intron_len = max_intron_len,
						aligner_extra_parameters = aligner_extra_parameters,
						alignment_resources = alignment_resources
					}
				}
			}
		}
		Array[AlignedSample] def_mm2 = select_first([minimap2_2pass_merged.aligned_sample, def_minimap_1pass])
    }

    if (aligner == "gmap") {
        call GMapIndex {
            input:
            reference = indexed_reference.fasta,
            runtime_attr_override = indexing_resources
        }
        scatter (sample in long_samples) {
            call gmap.wf_gmap {
                input:
                LRS = sample.LR,
                gmap_index = GMapIndex.gmap_index,
                name = sample.name,
                strand = sample.strand,
                score = sample.score,
                is_ref = sample.is_ref,
                exclude_redundant = sample.exclude_redundant,
                reference = indexed_reference.fasta,
                min_identity = min_identity,
                min_intron_len = select_first([min_intron_len,20]),
                max_intron_len = select_first([max_intron_len,2000]),
                max_intron_len_ends = select_first([max_intron_len_ends, 4000]),
                aligner_extra_parameters = aligner_extra_parameters,
                alignment_resources = alignment_resources
            }
        }
    }

    Array[AlignedSample] def_alignments = select_first([def_mm2, wf_gmap.aligned_sample])

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

	if (aligner == "2pass" || aligner == "2pass_merged") {
		Array[File] filtered_junctions = select_all(select_first([twopass_score.filtered_junctions]))
	}

    if (defined(JunctionUnion.combined_junctions) ||
        defined(twopass_junctions) ||
        defined(twopass_merge.merged_filtered_junctions)) {
        call MergeAllJunctions as final_junctions {
            input:
				combined_junctions = JunctionUnion.combined_junctions,
				twopass_indv = filtered_junctions,
				twopass_merged = twopass_merge.merged_filtered_junctions
        }
    }
    output {
        Array[AlignedSample] bams = def_alignments
        Array[File] summary_stats = summary_alignment_stats
        File summary_stats_table = CollectAlignmentStats.summary_stats_table
        File? junctions = final_junctions.merged_junctions
    }

}

task MergeAllJunctions {
	input {
		File? combined_junctions
		Array[File]? twopass_indv
		File? twopass_merged
	}

	output {
		File merged_junctions = "final_merged_junctions.bed"
	}

	command <<<
		cat ~{combined_junctions} ~{sep=" " twopass_indv} ~{twopass_merged} > merged_junctions.bed
		junctools convert -if bed -of ebed -s -d merged_junctions.bed -o final_merged_junctions.bed
	>>>
}

task Minimap2Index {
	input {
		Boolean is_hq
		File reference
		RuntimeAttr? indexing_resources
	}

	output {
		File index = basename(reference) + ".mmi"
	}

    Int cpus = 16
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 8,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([indexing_resources, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

	command <<<
	minimap2 -ax ~{if (is_hq) then "splice:hq" else "splice"} \
		-t ~{task_cpus} \
		-d ~{basename(reference)}.mmi ~{reference}

	>>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
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
        Array[File] gmap_index = glob("gmapIndex/reference/*")
    }


    command <<<
        set -euxo pipefail
        gmap_build -D gmapIndex -d reference ~{reference}
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

task JunctionUnion {
    input {
        Array[File] junction_files
    }

    output {
        File combined_junctions = "combined_junctions.bed"
    }

    command <<<
    cat ~{sep=" " junction_files} > "all.bed"
    junctools convert -if bed -of ebed -s -d all.bed -o combined_junctions.bed
    >>>
}

task gff2bed {
	input {
		File annotation
	}

	output {
		File bed = "annotation_junctions.bed"
	}

	command <<<
	gffread --bed ~{annotation} > annotation_junctions.bed
	>>>
}

task twopass_score {
	input {
		Array[File] alignment
		File reference
		File reference_fai
		String name
		String strand
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 16,
		mem_gb: 8,
		max_retries: 1,
		queue: ""
	}

	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	Int task_cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
	Int task_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + task_cpus/2)


	output {
		File filtered_junctions = name+".filt_junc.bed"
		File scored_junctions = name+".score_junc.bed"
	}

	command <<<
		strand_opt="--unstranded"
		if [ "~{strand}" == "fr-secondstrand" ]
		then
			strand_opt="--stranded"
		fi

		ln -s ~{reference} reference.fasta
		ln -s ~{reference_fai} reference.fasta.fai

		samtools merge --write-index aln.bam ~{sep=" " alignment}
		2passtools score --processes ~{task_cpus} ${strand_opt} -o ~{name + ".score_junc.bed"} -f reference.fasta aln.bam
		2passtools filter --exprs 'decision_tree_2_pred' -o tmp.filt_junc.bed ~{name + ".score_junc.bed"}
		junctools convert -if bed -of ebed <(junctools convert -if bed -of ibed tmp.filt_junc.bed) > ~{name + ".filt_junc.bed"}
	>>>

	runtime {
		cpu: task_cpus
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}
}

task twopass_merge {
	input {
		Array[File] scored_beds
		File reference
		File reference_fai
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 16,
		mem_gb: 8,
		max_retries: 1,
		queue: ""
	}

	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	Int task_cpus = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
	Int task_mem = round(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + task_cpus/2)


	output {
		File merged_filtered_junctions = "merged.filt_junc.bed"
		File merged_scored_junctions = "merged.score_junc.bed"
	}

	command <<<
		ln -s ~{reference} reference.fasta
		ln -s ~{reference_fai} reference.fasta.fai

		2passtools merge --processes ~{task_cpus} -o merged.score_junc.bed -f reference.fasta ~{sep=" " scored_beds}
		2passtools filter --exprs 'decision_tree_2_pred' -o tmp.filt_junc.bed merged.score_junc.bed
		junctools convert -if bed -of ebed <(junctools convert -if bed -of ibed tmp.filt_junc.bed) > merged.filt_junc.bed
	>>>

	runtime {
		cpu: task_cpus
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}
}
