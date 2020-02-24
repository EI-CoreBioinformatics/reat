version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_short {
    input {
        Array[AlignedSample] aligned_samples
#        File? annotation
    }

    scatter (aligned_sample in aligned_samples) {
        scatter (bam in aligned_sample.bam) {
            call Stringtie {
                input:
                aligned_sample = bam,
                strand = aligned_sample.strand
#                annotation = annotation
            }
        }
        call Merge {
            input:
            name = aligned_sample.name,
            aligner_name = aligned_sample.aligner,
            assemblies = Stringtie.assembled
        }
        AssembledSample stringtie_assembly = object { name: aligned_sample.name+"."+aligned_sample.aligner+".stringtie", strand: aligned_sample.strand, assembly: Merge.assembly}
    }

    scatter (aligned_sample in aligned_samples) {
        call Scallop {
            input:
            aligned_sample = aligned_sample
        }
    }

    output {
        Array[AssembledSample] assemblies = flatten([stringtie_assembly, Scallop.assembly])
        Array[File] fassemblies = flatten([Merge.assembly, Scallop.assembled])
    }
}

task Stringtie {
    input {
        File aligned_sample
        String prefix = basename(aligned_sample, ".bam")
        String strand
        File? annotation
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
        ~{"-G " + annotation} \
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