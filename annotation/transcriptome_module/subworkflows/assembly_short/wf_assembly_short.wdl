version 1.0

import "../common/tasks.wdl"
import "../common/structs.wdl"
import "../common/rt_struct.wdl"
import "wf_stringtie_short.wdl" as wf_sgt

workflow wf_assembly_short {
    input {
        Array[AlignedSample] aligned_samples
        File? reference_annotation
        Boolean? skip_scallop
        String? stringtie_extra_parameters
        String? scallop_extra_parameters
        RuntimeAttr? stringtie_assembly_resources
        RuntimeAttr? scallop_assembly_resources
    }

    Boolean def_skip_scallop = select_first([skip_scallop, false])

    String output_directory = "assembly_short"

    scatter (aligned_sample in aligned_samples) {
        call wf_sgt.wf_stringtie_short as stringtie {
                input:
                aligned_sample = aligned_sample,
                reference_annotation = reference_annotation,
                extra_parameters = stringtie_extra_parameters,
                output_directory = output_directory,
                runtime_attr_override = stringtie_assembly_resources
        }

        AssembledSample stringtie_assembly = object {
            name: aligned_sample.name+"."+aligned_sample.aligner+".stringtie", 
            strand: aligned_sample.strand, 
            score: aligned_sample.score,
            is_ref: aligned_sample.is_ref,
            exclude_redundant: aligned_sample.exclude_redundant,
            assembly: stringtie.gff
            }
    }

    if (!def_skip_scallop) {
        scatter (aligned_sample in aligned_samples) {
            call Scallop {
                input:
                aligned_sample = aligned_sample,
                output_directory = output_directory,
                extra_parameters = scallop_extra_parameters,
                runtime_attr_override = scallop_assembly_resources
            }

            AssembledSample scallop_assembly = object {
                name: aligned_sample.name+"."+aligned_sample.aligner+".scallop",
                strand: aligned_sample.strand,
                score: aligned_sample.score,
                is_ref: aligned_sample.is_ref,
                exclude_redundant: aligned_sample.exclude_redundant,
                assembly: Scallop.assembled
            }
        }
    }
    if (defined(scallop_assembly)) {
        Array[AssembledSample] all_assemblies = flatten([stringtie_assembly, select_first([scallop_assembly])])
    }

    if (! defined(scallop_assembly)) {
        Array[AssembledSample] no_scallop_assemblies = flatten([stringtie_assembly])
    }

    Array[AssembledSample] def_all_assemblies = select_first([all_assemblies, no_scallop_assemblies])

    scatter (assembly in select_all(stringtie_assembly)) {
        call tasks.TranscriptAssemblyStats as Stringtie_Stats{
            input:
            gff = assembly.assembly,
            output_directory = output_directory
        }
    }
    call tasks.TranscriptAssemblySummaryStats as Stringtie_Summary_Stats{
        input:
        stats = Stringtie_Stats.stats,
        output_prefix = "assembly_short.stringtie"
    }

    if (defined(scallop_assembly)) {
        scatter (assembly in select_first([scallop_assembly])) {
            call tasks.TranscriptAssemblyStats as Scallop_Stats {
                input:
                gff = assembly.assembly,
                output_directory = output_directory
            }
        }

        call tasks.TranscriptAssemblySummaryStats as Scallop_Summary_Stats{
            input:
            stats = Scallop_Stats.stats,
            output_prefix = "assembly_short.scallop"
        }
    }

    if (defined(scallop_assembly)) {
        Array[File] def_scallop_stats = select_first([Scallop_Stats.stats])
        Array[File] all_stats = flatten([Stringtie_Stats.stats, def_scallop_stats])
    }

    if (! defined(scallop_assembly)) {
        Array[File] no_scallop_stats = Stringtie_Stats.stats
    }


    output {
        Array[AssembledSample] assemblies = def_all_assemblies
        Array[File] stats = select_first([all_stats, no_scallop_stats])
        File stringtie_summary_stats = Stringtie_Summary_Stats.summary
        File? scallop_summary_stats = Scallop_Summary_Stats.summary
    }
}

task Scallop {
    input {
        String output_directory
        String? extra_parameters
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    Int cpus = 1

    output {
        File assembled = output_directory+"/"+aligned_sample.name+"."+aligned_sample.aligner+".scallop.gtf"
    }

    command <<<
        set -euxo pipefail
            case "~{aligned_sample.strand}" in
            fr-firststrand)
            strandness="--library_type first"
            ;;
            fr-secondstrand)
            strandness="--library_type second"
            ;;
            f)
            strandness="--library_type second"
            ;;
            r)
            strandness="--library_type first"
            ;;
            fr-unstranded)
            strandness="--library_type unstranded"
            ;;
        esac

        mkdir ~{output_directory}
        cd ~{output_directory}

        scallop --verbose 0 -i ~{sep=" " aligned_sample.bam} ~{extra_parameters} -o "~{aligned_sample.name+"."+aligned_sample.aligner}.scallop.gtf" "${strandness}"
        sed -i -e "s/\"gene./\"~{aligned_sample.name}_SCLP./g" "~{aligned_sample.name+"."+aligned_sample.aligner}.scallop.gtf"
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