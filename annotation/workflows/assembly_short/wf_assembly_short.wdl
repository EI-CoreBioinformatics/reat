version 1.0

import "../common/structs.wdl"

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
    }

    output {
        File assembled = aligned_sample.name+"."+aligned_sample.aligner+".stringtie.gtf"
        AssembledSample assembly = {"name": aligned_sample.name+"."+aligned_sample.aligner+".stringtie", "strand": aligned_sample.strand, "assembly": aligned_sample.name+"."+aligned_sample.aligner+".stringtie.gtf"}
    }

    command <<<
        case "~{aligned_sample.strand}" in
            fr-firststrand)
            strandness="--rf"
            ;;
            fr-secondstrand)
            strandness="--fr"
            ;;
        esac

        stringtie ~{aligned_sample.bam} \
        -p 4 \
        "${strandness}" \
        ~{"-G " + annotation} \
        -o "~{aligned_sample.name+"."+aligned_sample.aligner}.stringtie.gtf"
    >>>
}


# Needs to have the tool available... Not built yet for OSX
task Scallop {
    input {
        IndexedAlignedSample aligned_sample
    }

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
}