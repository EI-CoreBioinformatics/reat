version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_merge_long {
	input {
		File? reference_annotation
        AlignedSample aligned_sample
		String? extra_parameters
        String output_directory
        RuntimeAttr? runtime_attr_override
	}

	# Generate a list of chromosomes and regions to split the main job for parallel merging
	# bedtools merge
	# chunk regions
	call SplitInput {
		input:
		aligned_sample = aligned_sample
	}

	Array[File] regions = SplitInput.split_bams

	# scatter over the regions
	  # gffread
	scatter (region in regions) {
		call GffreadMerge {
			input:
			bam = region,
			extra_parameters = extra_parameters,
			runtime_attr_override = runtime_attr_override,
			output_directory = output_directory
		}
	}
	# Join the resuls either using 'cat' or 'gffread' again ensuring the output is sorted

	call JoinResults {
		input:
		gff_files = GffreadMerge.gff,
		output_file = aligned_sample.name+"."+aligned_sample.aligner+".gffread_merge.gtf"
	}

	output {
		File gff = JoinResults.gff
	}
}

task JoinResults {
	input {
		Array[File] gff_files
		String output_file
	}

	output {
		File gff = output_file
	}

	command <<<
		cat ~{sep=" " gff_files} > ~{output_file}
	>>>
}

task SplitInput {
	input {
		AlignedSample aligned_sample
		Int num_reads = 100000
	}

	output {
		Array[File] split_bams = glob('split_bundles/*.bam')
	}

	command <<<
		samtools merge - ~{sep=' ' aligned_sample.bam} | seqkit bam -N ~{num_reads} -o split_bundles
	>>>
}

task GffreadMerge {
    input {
        String output_directory
        File bam
        String? extra_parameters
        RuntimeAttr? runtime_attr_override
    }

    output {
        File gff = output_directory + "/" + basename(bam) + ".gtf"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        sam2gff -s ~{bam} | gffread -T -M -K ~{extra_parameters} -o ~{basename(bam)}.gtf
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
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