version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow wf_assembly_short {
    input {
        Array[IndexedAlignedSample] aligned_samples
        File? annotation
    }

    scatter (aligned_sample in aligned_samples) {
        call Stringtie{
            input:
            aligned_sample = aligned_sample,
            annotation = annotation
        }

        call Scallop {
            input:
            aligned_sample = aligned_sample
        }
    }

    output {
        Array[AssembledSample] assemblies = flatten([Stringtie.assembly, Scallop.assembly])
        Array[File] fassemblies = flatten([Stringtie.assembled, Scallop.assembled])
    }
}

task Stringtie {
    input {
        IndexedAlignedSample aligned_sample
        File? annotation
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8

    output {
        File assembled = aligned_sample.name+"."+aligned_sample.aligner+".stringtie.gtf"
        AssembledSample assembly = {"name": aligned_sample.name+"."+aligned_sample.aligner+".stringtie", "strand": aligned_sample.strand, "assembly": aligned_sample.name+"."+aligned_sample.aligner+".stringtie.gtf"}
    }

    command <<<
        case ~{aligned_sample.strand} in
            fr-firststrand)
            strandness="--rf"
            ;;
            fr-secondstrand)
            strandness="--fr"
            ;;
        esac

        stringtie ~{aligned_sample.bam} \
        -p "~{cpus}" \
        "${strandness}" \
        ~{"-G " + annotation} \
        -o "~{aligned_sample.name+"."+aligned_sample.aligner}.stringtie.gtf"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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


# Needs to have the tool available... Not built yet for OSX
task Scallop {
    input {
        IndexedAlignedSample aligned_sample
        RuntimeAttr? runtime_attr_override
    }
    
    Int cpus = 2

    output {
        File assembled = aligned_sample.name+"."+aligned_sample.aligner+".scallop.gtf"
        AssembledSample assembly = {"name": aligned_sample.name+"."+aligned_sample.aligner+".scallop", "strand": aligned_sample.strand, "assembly": aligned_sample.name+"."+aligned_sample.aligner+".scallop.gtf"}
    }

    command <<<
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

        scallop --verbose 0 -i ~{aligned_sample.bam} -o "~{aligned_sample.name+"."+aligned_sample.aligner}.scallop.gtf" "${strandness}"
    >>>

    RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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