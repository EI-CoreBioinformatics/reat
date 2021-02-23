version 1.0

import "../common/rt_struct.wdl"
import "../common/structs.wdl"

workflow wf_gmap {
	input {
		File reference
		Array[File] gmap_index
		Array[File]+ LRS
        Float min_identity = 0.9
        Int? min_intron_len = 20
        Int? max_intron_len = 200000
        Int? max_intron_len_ends = 100000
		String strand
		String name
		Int? score
		Boolean? is_ref
		Boolean? exclude_redundant
        String? aligner_extra_parameters
        RuntimeAttr? alignment_resources
	}

	scatter (LR in LRS) {
		call GMapLong {
			input:
			LR = LR,
			gmap_index = gmap_index,
			strand = strand,
			name = name,
			reference = reference,
			min_identity = min_identity,
			min_intron_len = select_first([min_intron_len,20]),
			max_intron_len = select_first([max_intron_len,2000]),
			max_intron_len_ends = select_first([max_intron_len_ends, 4000]),
			extra_parameters = aligner_extra_parameters,
			runtime_attr_override = alignment_resources
		}
	}

	output {
		AlignedSample aligned_sample = object {name: name, strand:strand, aligner:"gmap", bam: GMapLong.bam, merge:false,
										score: select_first([score, 1]),
										is_ref: select_first([is_ref, false]),
										exclude_redundant: select_first([exclude_redundant, false])}
	}
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
