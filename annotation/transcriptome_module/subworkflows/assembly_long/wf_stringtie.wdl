version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_stringtie_long {
	input {
		File? reference_annotation
        Boolean collapse = false
        String collapse_string = if (collapse==true) then "-R " else "-L "
        AlignedSample aligned_sample
        String output_directory
        RuntimeAttr? runtime_attr_override
	}

	scatter (bam in aligned_sample.bam) {
		call StringtieLong {
			input:
			name = aligned_sample.name,
			aligner = aligned_sample.aligner,
			bam = bam,
			reference_annotation = reference_annotation,
			output_directory = "assembly_long",
			collapse = collapse,
			runtime_attr_override = runtime_attr_override
		}
	}

	if (length(StringtieLong.gff) > 1) {
		call StringtieMerge {
			input:
			label = aligned_sample.name + "_" + aligned_sample.aligner,
			assemblies = StringtieLong.gff
		}
	}

	output {
		File gff = select_first([StringtieMerge.gff, StringtieLong.gff[0]])
	}
}


task StringtieMerge {
	input {
		String label
		Array[File] assemblies
	}

	output {
		File gff = "merged.gff"
	}

	command <<<
		stringtie --merge -l ~{label}  ~{sep=" " assemblies} -o merged.gtf
	>>>
}

task StringtieLong {
    input {
		File? reference_annotation
        Boolean collapse = false
		String name
		String aligner
        File bam
        String output_directory
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File gff = output_directory + "/" + name+"."+aligner+".stringtie.gtf"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        stringtie ~{bam} -p "~{cpus}"~{" -G " + reference_annotation} ~{if(collapse) then "-R " else "-L "} -o "~{name}.~{aligner}.stringtie.gtf"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
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
