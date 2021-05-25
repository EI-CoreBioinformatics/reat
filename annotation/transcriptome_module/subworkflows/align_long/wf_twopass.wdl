version 1.0

import "../common/rt_struct.wdl"
import "../common/structs.wdl"

workflow wf_twopass {
	input {
		File reference
		File reference_fai
		File reference_index
		Array[File]+ LRS
		File? bed_junctions
		Int? max_intron_len = 200000
		Boolean is_hq
		Boolean merge_juncs = false
		String strand
		String name
		Int? score
		Boolean? is_ref
		Boolean? exclude_redundant
		String? aligner_extra_parameters
        RuntimeAttr? twopass_resources
        RuntimeAttr? twopass_merge_resources
		RuntimeAttr? alignment_resources
	}

	scatter (LR in LRS) {
		call Minimap2Pass as base_alignment {
			input:
			LR = LR,
			is_hq = is_hq,
			strand = strand,
			name = name,
			reference = reference_index,
			bed_junctions = bed_junctions,
			max_intron_len = select_first([max_intron_len, 200000]),
			extra_parameters = aligner_extra_parameters,
			runtime_attr_override = alignment_resources
		}

		call twopass {
			input:
			strand = strand,
			alignment = base_alignment.bam,
			index = base_alignment.csi,
			reference = reference,
			reference_fai = reference_fai,
			runtime_attr_override = twopass_resources
		}

		if (!merge_juncs) {
			call Minimap2Pass as second_pass {
				input:
				LR = LR,
				is_hq = is_hq,
				strand = strand,
				name = name,
				reference = reference_index,
				bed_junctions = twopass.filtered_junctions,
				max_intron_len = select_first([max_intron_len, 200000]),
				extra_parameters = aligner_extra_parameters,
				runtime_attr_override = alignment_resources
			}
		}
	}

	if (merge_juncs) {
		call twopass_merge {
			input:
			reference = reference,
			reference_fai = reference_fai,
			strand = strand,
			scored_beds = twopass.scored_junctions,
			runtime_attr_override = twopass_merge_resources
		}

		scatter (LR in LRS) {
			call Minimap2Pass as merged_second_pass {
				input:
				LR = LR,
				is_hq = is_hq,
				strand = strand,
				name = name,
				reference = reference_index,
				bed_junctions = twopass_merge.merged_filtered_junctions,
				max_intron_len = select_first([max_intron_len, 200000]),
				extra_parameters = aligner_extra_parameters,
				runtime_attr_override = alignment_resources
			}
		}
	}

	if (!merge_juncs) {
		Array[File] non_merged_junc_bams = select_all(second_pass.bam)
	}
	if (merge_juncs) {
		Array[File] merged_junc_bams = select_first([merged_second_pass.bam])
	}

	if (!defined(twopass_merge.merged_filtered_junctions)) {
		call CombineJunctions {
			input:
			junctions = twopass.filtered_junctions
		}
	}

	String aligner_name = if merge_juncs then "minimap2-2pass_merged" else "minimap2-2pass"

	call ReformatJunctions {
		input:
		junctions = select_first([twopass_merge.merged_filtered_junctions, CombineJunctions.combined_junctions]),
		reference = reference
	}

	output {
		AlignedSample aligned_sample = object {name: name, strand:strand, aligner:aligner_name, bam: select_first([non_merged_junc_bams, merged_junc_bams]), merge:false,
									   is_ref: select_first([is_ref, false]), exclude_redundant: select_first([exclude_redundant, false]),
									   score: select_first([score, 0])}
		File filtered_unique_junctions = ReformatJunctions.final_junctions
	}
}

task Minimap2Pass {
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
		File bam = "alignments/minimap2."+name+"."+LR_basename+".bam"
		File csi = "alignments/minimap2."+name+"."+LR_basename+".bam.csi"
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
		-L \
		--eqx -2 \
		--secondary=no \
		~{reference} - | samtools view -F 4 -F 0x900 -bS - | samtools sort --write-index -@ ~{task_cpus/2} -m 1G --reference ~{reference} -T minimap2.sort -o minimap2.~{name}.~{LR_basename}.bam -
	>>>

	runtime {
		cpu: task_cpus
		memory: task_mem + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}
}

task twopass {
	input {
		File alignment
		File index
		File reference
		File reference_fai
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
		File filtered_junctions = basename(alignment)+".filt_junc.bed"
		File scored_junctions = basename(alignment)+".score_junc.bed"
	}

	command <<<
		strand_opt="--unstranded"
		if [ "~{strand}" == "fr-secondstrand" ]
		then
			strand_opt="--stranded"
		fi

		ln -s ~{reference} reference.fasta
		ln -s ~{reference_fai} reference.fasta.fai

		ln -s ~{alignment} aln.bam
		ln -s ~{index} aln.bam.csi
		2passtools score --processes ~{task_cpus} ${strand_opt} -o ~{basename(alignment) + ".score_junc.bed"} -f reference.fasta aln.bam
		2passtools filter --exprs 'decision_tree_2_pred' -o ~{basename(alignment) + ".filt_junc.bed"} ~{basename(alignment) + ".score_junc.bed"}
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
		File merged_filtered_junctions = "merged.filt_junc.bed"
		File merged_scored_junctions = "merged.score_junc.bed"
	}

	command <<<
		strand_opt="--unstranded"
		if [ "~{strand}" == "fr-secondstrand" ]
		then
			strand_opt="--stranded"
		fi

		ln -s ~{reference} reference.fasta
		ln -s ~{reference_fai} reference.fasta.fai

		2passtools merge --processes ~{task_cpus} -o merged.score_junc.bed -f reference.fasta ~{sep=" " scored_beds}
		2passtools filter --exprs 'decision_tree_2_pred' -o merged.filt_junc.bed merged.score_junc.bed
	>>>

	runtime {
		cpu: task_cpus
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
		queue: select_first([runtime_attr.queue, default_attr.queue])
	}
}

task CombineJunctions {
	input {
		Array[File] junctions
	}

	output {
		File combined_junctions = "2pass_indv_combined.bed"
	}

	command <<<
		cat ~{sep=" " junctions} | sort -k1,1V -k2,3n | uniq > 2pass_indv_combined.bed
	>>>
}

task ReformatJunctions {
	input {
		File junctions
		File reference
	}

	output {
		File final_junctions = "final_2pass_junctions.bed"
	}

	command <<<
		cat ~{junctions} > final_2pass_junctions.bed
	>>>
}