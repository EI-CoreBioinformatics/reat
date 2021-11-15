version 1.0

import "../common/rt_struct.wdl"
import "../common/structs.wdl"

workflow wf_mm2 {
	input {
		File reference
		Array[File]+ LRS
		File? bed_junctions
        Int? max_intron_len = 200000
		Boolean is_hq
		String strand
		String name
		Int? score
		Boolean? is_ref
		Boolean? exclude_redundant
        String? aligner_extra_parameters
        RuntimeAttr? alignment_resources
	}

	scatter (LR in LRS) {
		call Minimap2Long {
			input:
			LR = LR,
			is_hq = is_hq,
			strand = strand,
			name = name,
			reference = reference,
			bed_junctions = bed_junctions,
			max_intron_len = select_first([max_intron_len, 200000]),
			extra_parameters = aligner_extra_parameters,
			runtime_attr_override = alignment_resources
		}
	}

	output {
		AlignedSample aligned_sample = object {name: name, strand:strand, aligner:"minimap2", bam: Minimap2Long.bam, merge:false,
									   is_ref: select_first([is_ref, false]), exclude_redundant: select_first([exclude_redundant, false]),
									   score: select_first([score, 0])}
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
        ~{reference} - | samtools view -F 4 -F 0x900 -bS - | samtools sort -@ ~{task_cpus/2} -m 1G -T minimap2.sort -o minimap2.~{name}.~{LR_basename}.bam -
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
		File reference
	}

	output {
		File junctions = basename(alignment)+".junc.bed"
	}

	command <<<
		2passtools score -o ~{basename(alignment) + ".score.bed"} -f ~{reference} ~{alignment}
		2passtools filter --exprs 'decision_tree_2_pred' -o ~{basename(alignment) + ".junc.bed"} ~{basename(alignment) + ".score.bed"}


	>>>
}