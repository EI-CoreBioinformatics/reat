version 1.0

import "../common/tasks.wdl"
import "../common/structs.wdl"
import "../common/rt_struct.wdl"
import "wf_stringtie.wdl" as wstl
import "wf_merge.wdl" as wmrg

workflow wf_assembly_long {
    input {
        File? reference_annotation
        Array[AlignedSample] aligned_samples
        String assembler = "filter"
        String stats_output_prefix
        String? assembler_extra_parameters
        Int? min_coverage
        Int? min_identity
        RuntimeAttr? assembly_resources
    }
    String output_directory = "assembly_long"

    scatter (sample in aligned_samples) {
        if (assembler == "filter") {
            call Sam2gff {
                input:
                aligned_sample = sample,
                min_coverage = min_coverage,
                min_identity = min_identity,
                runtime_attr_override = assembly_resources,
                output_directory = output_directory
            }
        }
        if (assembler == "merge") {
                call wmrg.wf_merge_long as GffreadMerge {
                    input:
                    aligned_sample = sample,
                    extra_parameters = assembler_extra_parameters,
                    runtime_attr_override = assembly_resources,
                    output_directory = output_directory
                }
        }
        if (assembler == "stringtie") {
            call wstl.wf_stringtie_long as stringtie_assemble {
                input:
                reference_annotation = reference_annotation,
                aligned_sample = sample,
                extra_parameters = assembler_extra_parameters,
                runtime_attr_override = assembly_resources,
                output_directory = output_directory
            }
        }
        if (assembler == "stringtie_collapse") {
            call wstl.wf_stringtie_long as stringtie_collapse {
                input:
                reference_annotation = reference_annotation,
                collapse = true,
                aligned_sample = sample,
                extra_parameters = assembler_extra_parameters,
                runtime_attr_override = assembly_resources,
                output_directory = output_directory
            }
        }

        File def_gff = select_first([Sam2gff.gff, GffreadMerge.gff, stringtie_assemble.gff, stringtie_collapse.gff])
        AssembledSample assembled_long = object { name: sample.name+"."+sample.aligner+"."+assembler, strand: sample.strand, assembly: def_gff,
                                         score: sample.score,
                                         is_ref: sample.is_ref,
                                         exclude_redundant: sample.exclude_redundant}
        if (assembler != "filter") {
            call tasks.TranscriptAssemblyStats {
                input:
                gff = def_gff,
                output_directory = output_directory
            }
        }
    }
    
    Array[File] assembled_samples_stats = select_all(TranscriptAssemblyStats.stats)

    # If at least one sample was assembled then generate the summary stats.
    # Otherwise, since all samples were just filtered, there will be no stats present in the array
    if (length(assembled_samples_stats)>=1) {
        call tasks.TranscriptAssemblySummaryStats {
            input:
            stats = assembled_samples_stats,
            output_prefix = "assembly_long." + stats_output_prefix
        }
    }

    output {
        Array[AssembledSample] gff = assembled_long
        Array[File] stats = assembled_samples_stats
        File? summary_stats = TranscriptAssemblySummaryStats.summary
    }
}

task Sam2gff {
    input {
        String output_directory
        Int? min_coverage = 80
        Int? min_identity = 95
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    output {
        File gff = output_directory + "/" + aligned_sample.name+"."+aligned_sample.aligner+".sam2gff.gtf"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff --gtf -s ~{aligned_sample.name} \
        -u ~{aligned_sample.name}.~{aligned_sample.aligner}.unfiltered.gtf ~{'--min_coverage=' + min_coverage} ~{'--min_identity=' + min_identity} > ~{aligned_sample.name}.~{aligned_sample.aligner}.sam2gff.gtf
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

task GffreadMerge {
    input {
        String output_directory
        AlignedSample aligned_sample
        String? extra_parameters
        RuntimeAttr? runtime_attr_override
    }

    output {
        File gff = output_directory + "/" + aligned_sample.name+"."+aligned_sample.aligner+".gffread_merge.gtf"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        for bam in ~{sep=" " aligned_sample.bam}; do
        samtools view -F 4 -F 0x900 $bam; done | sam2gff -s ~{aligned_sample.name} | gffread -T -M -K ~{extra_parameters} -o ~{aligned_sample.name}.~{aligned_sample.aligner}.gffread_merge.gtf
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