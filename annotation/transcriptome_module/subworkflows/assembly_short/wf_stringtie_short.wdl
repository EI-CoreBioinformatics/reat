version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_stringtie_short {
	input {
		File? reference_annotation
        Boolean collapse = false
        String collapse_string = if (collapse==true) then "-R " else "-L "
        AlignedSample aligned_sample
		String? extra_parameters
        String output_directory
        RuntimeAttr? runtime_attr_override
	}

	scatter (bam in aligned_sample.bam) {
		call Stringtie {
			input:
			aligned_sample = bam,
			strand = aligned_sample.strand,
			extra_parameters = extra_parameters,
			output_directory = output_directory,
			runtime_attr_override = runtime_attr_override
		}
	}

	Array[File] stringtie_assemblies = select_first([Stringtie.assembled])
	# If this contains more than one ReadPair, and has been marked as "merge",
	# then it will contain a single bam file as the result of the
	# alignment part of the workflow, so it will not need merging assemblies.
	# Otherwise, this conditional merges the assemblies generated for the multiple bam files
	if (length(stringtie_assemblies) > 1) {
		call Merge {
			input:
			name = aligned_sample.name,
			aligner_name = aligned_sample.aligner,
			assemblies = stringtie_assemblies,
			output_directory = output_directory,
			runtime_attr_override = runtime_attr_override
		}
	}

	File def_stringtie_assembly = select_first([Merge.assembly, stringtie_assemblies[0]]) # Indexing the array is OK because we expect a single item in it

	output {
		File gff = def_stringtie_assembly
	}
}

task Merge {
    input {
        String name
        String aligner_name
        String output_directory
        Array[File] assemblies
        File? annotation
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 1

    output {
        File assembly = output_directory + "/" + name+"."+aligner_name+".stringtie.gtf"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        stringtie --merge \
        ~{"-G " + annotation} \
        -l ~{name}"_STRG" \
        -o "~{name+"."+aligner_name}.stringtie.gtf" \
        ~{sep=" " assemblies}
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

task Stringtie {
    input {
        File aligned_sample
        String output_directory
        String prefix = basename(aligned_sample, ".bam")
        String strand
        String? extra_parameters
        File? reference_annotation
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
        mem_gb: 16,
        max_retries: 1,
        queue: ""
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int task_cpus = select_first([runtime_attr.cpu_cores, cpus])

    output {
        File assembled = output_directory + "/"+prefix+".stringtie.gtf"
    }

    command <<<
        set -euxo pipefail
        strandness=""
        case ~{strand} in
            fr-firststrand)
            strandness="--rf"
            ;;
            fr-secondstrand)
            strandness="--fr"
            ;;
        esac

        mkdir ~{output_directory}
        cd ~{output_directory}
        stringtie ~{aligned_sample} \
        -p "~{task_cpus}" \
        ${strandness} \
        ~{"-G " + reference_annotation} \
        ~{extra_parameters} -o "~{prefix}.stringtie.gtf"
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
    }
}
