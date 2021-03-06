version 1.0

import "../common/tasks.wdl"
import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_short {
    input {
        Array[AlignedSample] aligned_samples
        File? reference_annotation
        Boolean? skip_scallop
        RuntimeAttr? stringtie_assembly_resources
        RuntimeAttr? scallop_assembly_resources
    }

    Boolean def_skip_scallop = select_first([skip_scallop, false])

    String output_directory = "assembly_short"

    scatter (aligned_sample in aligned_samples) {
        # There's some issue with scatter of scatter and optional parameters, 
        # the following if construct is to eliminate the parameter from the WF if it is not present
        # working around this issue
        # TODO: Check if this is still the case
        if (defined(reference_annotation)) {
            scatter (bam in aligned_sample.bam) {
                call Stringtie as annot{
                    input:
                    aligned_sample = bam,
                    strand = aligned_sample.strand,
                    reference_annotation = reference_annotation,
                    output_directory = output_directory,
                    runtime_attr_override = stringtie_assembly_resources
                }
            }
        }

        if (!defined(reference_annotation)) {
            scatter (bam in aligned_sample.bam) {
                call Stringtie as no_annot {
                    input:
                    aligned_sample = bam,
                    strand = aligned_sample.strand,
                    output_directory = output_directory,
                    runtime_attr_override = stringtie_assembly_resources
                }
            }
        }

        Array[File] stringtie_assemblies = select_first([no_annot.assembled, annot.assembled])
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
                runtime_attr_override = stringtie_assembly_resources
            }
        }

        File def_stringtie_assembly = select_first([Merge.assembly, stringtie_assemblies[0]]) # Indexing the array is OK because we expect a single item in it
        AssembledSample stringtie_assembly = object {
            name: aligned_sample.name+"."+aligned_sample.aligner+".stringtie", 
            strand: aligned_sample.strand, 
            score: aligned_sample.score,
            is_ref: aligned_sample.is_ref,
            exclude_redundant: aligned_sample.exclude_redundant,
            assembly: def_stringtie_assembly
            }
    }

    if (!def_skip_scallop) {
        scatter (aligned_sample in aligned_samples) {
            call Scallop {
                input:
                aligned_sample = aligned_sample,
                output_directory = output_directory,
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

task Stringtie {
    input {
        File aligned_sample
        String output_directory
        String prefix = basename(aligned_sample, ".bam")
        String strand
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
        -o "~{prefix}.stringtie.gtf"
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        queue: select_first([runtime_attr.queue, default_attr.queue])
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

task Scallop {
    input {
        String output_directory
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

        scallop --verbose 0 -i ~{sep=" " aligned_sample.bam} -o "~{aligned_sample.name+"."+aligned_sample.aligner}.scallop.gtf" "${strandness}"
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