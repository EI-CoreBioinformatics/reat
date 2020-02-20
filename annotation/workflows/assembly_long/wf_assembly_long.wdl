version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_long {
    input {
        File? reference_annotation
        Array[AlignedSample] aligned_samples
        String assembler = "None"
    }

    scatter (sample in aligned_samples) {
        if (assembler == "None") {
            call sam2gff {
                input:
                aligned_sample = sample
            }
        }
        if (assembler == "merge") {
            call gffread_merge {
                input:
                aligned_sample = sample
            }
        }

        if (assembler == "stringtie") {
            call stringtie_long {
                input:
                reference_annotation = reference_annotation,
                aligned_sample = sample
            }
        }
        File def_gff = select_first([sam2gff.gff, gffread_merge.gff, stringtie_long.gff])
    }

    output {
        Array[File] gff = def_gff
    }

}

task stringtie_long {
    input {
        File? reference_annotation
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File gff = "result.gff"
    }

    command <<<
        set -euxo pipefail
    stringtie -p "~{cpus}" ~{"-G " + reference_annotation} -L ~{sep=" " aligned_sample.bam} -o "result.gff"
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

task sam2gff {
    input {
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File gff = "result.gff"
    }

    command <<<
        set -euxo pipefail
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff -s ~{aligned_sample.name} > result.gff
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

task gffread_merge {
    input {
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }

    output {
        File gff = "result.gff"
    }

    command <<<
        set -euxo pipefail
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff -s ~{aligned_sample.name} | gffread -T -M -K -o result.gff
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}