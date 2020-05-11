version 1.0

import "../common/tasks.wdl"
import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_short {
    input {
        Array[AlignedSample] aligned_samples
        File? reference_annotation
        RuntimeAttr? assembly_resources
    }

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
                    runtime_attr_override = assembly_resources
                }
            }
        }

        if (!defined(reference_annotation)) {
            scatter (bam in aligned_sample.bam) {
                call Stringtie as no_annot {
                    input:
                    aligned_sample = bam,
                    strand = aligned_sample.strand,
                    runtime_attr_override = assembly_resources
                }
            }
        }

        Array[File] stringtie_assemblies = select_first([no_annot.assembled, annot.assembled])
        if (length(stringtie_assemblies) > 1) {
            call Merge {
                input:
                name = aligned_sample.name,
                aligner_name = aligned_sample.aligner,
                assemblies = stringtie_assemblies,
                runtime_attr_override = assembly_resources
            }
        }

        File def_stringtie_assembly = select_first([Merge.assembly, stringtie_assemblies[0]]) # Indexing the array is OK because we expect a single item in it
        AssembledSample stringtie_assembly = object { name: aligned_sample.name+"."+aligned_sample.aligner+".stringtie", strand: aligned_sample.strand, assembly: def_stringtie_assembly}
    }

    scatter (aligned_sample in aligned_samples) {
        call Scallop {
            input:
            aligned_sample = aligned_sample,
            runtime_attr_override = assembly_resources
        }
    }
    Array[AssembledSample] all_assemblies = flatten([stringtie_assembly, Scallop.assembly])

    scatter (assembly in stringtie_assembly) {
        call tasks.TranscriptAssemblyStats as Stringtie_Stats{
            input:
            gff = assembly.assembly
        }
    }

    scatter (assembly in Scallop.assembly) {
        call tasks.TranscriptAssemblyStats as Scallop_Stats {
            input:
            gff = assembly.assembly
        }
    }

    call tasks.TranscriptAssemblySummaryStats as Stringtie_Summary_Stats{
        input:
        stats = Stringtie_Stats.stats,
        output_prefix = "stringtie"
    }

    call tasks.TranscriptAssemblySummaryStats as Scallop_Summary_Stats{
        input:
        stats = Scallop_Stats.stats,
        output_prefix = "scallop"
    }

    output {
        Array[AssembledSample] assemblies = all_assemblies
        Array[File] stats = flatten([Stringtie_Stats.stats, Scallop_Stats.stats])
        File stringtie_summary_stats = Stringtie_Summary_Stats.summary
        File scallop_summary_stats = Scallop_Summary_Stats.summary
    }
}

task Stringtie {
    input {
        File aligned_sample
        String prefix = basename(aligned_sample, ".bam")
        String strand
        File? reference_annotation
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File assembled = prefix+".stringtie.gtf"
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

        stringtie ~{aligned_sample} \
        -p "~{cpus}" \
        "${strandness}" \
        ~{"-G " + reference_annotation} \
        -o "~{prefix}.stringtie.gtf"
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

task Merge {
    input {
        String name
        String aligner_name
        Array[File] assemblies
        File? annotation
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File assembly = name+"."+aligner_name+".stringtie.gtf"
    }

    command <<<
        set -euxo pipefail
        stringtie --merge \
        ~{"-G " + annotation} \
        -o "~{name+"."+aligner_name}.stringtie.gtf" \
        ~{sep=" " assemblies}
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

task Scallop {
    input {
        AlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    Int cpus = 2

    output {
        File assembled = aligned_sample.name+"."+aligned_sample.aligner+".scallop.gtf"
        AssembledSample assembly = {"name": aligned_sample.name+"."+aligned_sample.aligner+".scallop", "strand": aligned_sample.strand, "assembly": aligned_sample.name+"."+aligned_sample.aligner+".scallop.gtf"}
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

        scallop --verbose 0 -i ~{sep=" " aligned_sample.bam} -o "~{aligned_sample.name+"."+aligned_sample.aligner}.scallop.gtf" "${strandness}"
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