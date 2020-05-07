version 1.0

import "../common/tasks.wdl"
import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_long {
    input {
        File? reference_annotation
        Array[AlignedSample] aligned_samples
        String assembler = "filter"
        RuntimeAttr? assembly_resources
    }

    scatter (sample in aligned_samples) {
        if (assembler == "filter") {
            call Sam2gff {
                input:
                aligned_sample = sample,
                runtime_attr_override = assembly_resources
            }
            call FilterGFF {
                input:
                gff = Sam2gff.gff
            }
        }
        if (assembler == "merge") {
            call GffreadMerge {
                input:
                aligned_sample = sample,
                runtime_attr_override = assembly_resources
            }
        }
        if (assembler == "stringtie") {
            call StringtieLong as stringtie_assemble {
                input:
                reference_annotation = reference_annotation,
                aligned_sample = sample,
                runtime_attr_override = assembly_resources
            }
        }
        if (assembler == "stringtie_collapse") {
            call StringtieLong as stringtie_collapse {
                input:
                reference_annotation = reference_annotation,
                collapse = true,
                aligned_sample = sample,
                runtime_attr_override = assembly_resources
            }
        }

        File def_gff = select_first([FilterGFF.filtered_gff, GffreadMerge.gff, stringtie_assemble.gff, stringtie_collapse.gff])
        AssembledSample assembled_long = object { name: sample.name+"."+sample.aligner+assembler, strand: sample.strand, assembly: def_gff}
        if (assembler != "filter") {
            call tasks.Stats {
                input:
                gff = def_gff
            }
        }
    }
    
    Array[File] assembled_samples_stats = select_all(Stats.stats)

    if (length(assembled_samples_stats)>1) {
        call tasks.SummaryStats {
            input:
            stats = assembled_samples_stats
        }
    }

    output {
        Array[AssembledSample] gff = assembled_long
        Array[File?] stats = Stats.stats
        File? summary_stats = SummaryStats.summary
    }
}

task FilterGFF {
    input {
        File gff
        String min_coverage = "80"
        String min_identity = "95"
    }

    output {
        File filtered_gff = basename(gff)+"."+min_identity+"id"+min_coverage+"cov.gff"
    }

    command <<<
    filter_gmap_hardFilter_v0.1.pl --gff ~{gff} --identity ~{min_identity} --coverage ~{min_coverage} > ~{basename(gff)}.~{min_identity}id~{min_coverage}cov.gff
    >>>
}

task StringtieLong {
    input {
        File? reference_annotation
        Boolean collapse = false
        String collapse_string = if (collapse==true) then "-R " else "-L "
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File gff = aligned_sample.name+"."+aligned_sample.aligner+".stringtie.gtf"
    }

    command <<<
        set -euxo pipefail
        stringtie -p "~{cpus}" ~{"-G " + reference_annotation} ~{collapse_string} ~{sep=" " aligned_sample.bam} -o "~{aligned_sample.name}.~{aligned_sample.aligner}.stringtie.gtf"
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

task Sam2gff {
    input {
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File gff = aligned_sample.name+"."+aligned_sample.aligner+".sam2gff.gff"
    }

    command <<<
        set -euxo pipefail
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff -s ~{aligned_sample.name} > ~{aligned_sample.name}.~{aligned_sample.aligner}.sam2gff.gff
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

task GffreadMerge {
    input {
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }

    output {
        File gff = aligned_sample.name+"."+aligned_sample.aligner+".gffread_merge.gtf"
    }

    command <<<
        set -euxo pipefail
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff -s ~{aligned_sample.name} | gffread -T -M -K -o ~{aligned_sample.name}.~{aligned_sample.aligner}.gffread_merge.gtf
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