version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow portcullis {
    input {
        File reference
        Map[String,Array[String]]? group_to_samples # Consider applying "localization_optional" to the preparation task
        Array[AlignedSample] aligned_samples
        File? annotation
    }

    parameter_meta {
        reference: "Genome reference file."
        annotation: "Reference annotation of junctions in BED format."
        sample_groups: "Pre-defined groupings for the aligned_samples, these need to match the sample names. Each portcullis run groups one or more samples."
        aligned_samples: "List of aligned samples."
    }

    if (defined(annotation)) {
        File def_annotation = select_first([annotation])
        call PrepareRef {
            input:
            annotation = def_annotation
        }
    }

    if (defined(group_to_samples)) {
        call PrepareGroups {
            input:
            group_to_samples = select_first([group_to_samples]),
            aligned_samples = aligned_samples
        }
        scatter (sample in PrepareGroups.portcullis_grouped_sample) {
            call Full as groupedFull{
                input:
                reference = reference,
                sample = sample,
                reference_bed = PrepareRef.refbed
            }
        }
        Array[File] def_grouped = groupedFull.pass
    }

    if (!defined(group_to_samples)) {
        scatter (aligned_sample in aligned_samples) {
            call Full {
                input:
                reference = reference,
                sample = aligned_sample,
                reference_bed = PrepareRef.refbed
            }
        }
        Array[File] def_ungrouped = Full.pass
    }

    Array[File] to_merge = select_first([def_grouped, def_ungrouped])

    call Merge {
        input:
        tabs = to_merge
    }

    output {
        File tab = Merge.tab
        File bed = Merge.bed
        File gff3 = Merge.gff3
    }
}

task PrepareGroups {
    input {
        Map[String,Array[String]] group_to_samples
        Array[AlignedSample] aligned_samples
    }

    output {
        Array[AlignedSample] portcullis_grouped_sample = read_json(stdout())
    }

    command <<<
        portcullis_sample_grouping ~{write_objects(aligned_samples)} ~{write_map(group_to_samples)}
    >>>
}

task Full {
    input {
        File reference
        File? reference_bed
        AlignedSample sample
        Int min_coverage = 2
        RuntimeAttr? runtime_attr_override
    }

    output {
        File pass = "portcullis_out/3-filt/portcullis_filtered.pass.junctions.tab"
    }

    command <<<
    set -euxo pipefail
    strandedness="UNKNOWN"
    case "~{sample.strand}" in
        fr-firststrand)
        strandedness="firststrand"
        ;;
        fr-secondstrand)
        strandedness="secondstrand"
        ;;
        fr-unstranded)
        strandedness="unstranded"
        ;;
    esac
    portcullis full --save_bad --exon_gff --intron_gff --strandedness=$strandedness \
    ~{"-r " + reference_bed} --canonical=C,S --min_cov=~{min_coverage} \
    ~{reference} ~{sep=" " sample.bam}
    >>>
}

task PrepareRef {
    input {
        File annotation
        RuntimeAttr? runtime_attr_override
    }
    
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

    output {
        File refbed = "reference.refbed"
    }

    command <<<
        set -euxo pipefail
        junctools convert -if gtf -of ebed -o "reference.refbed" ~{annotation}
    >>>

}

task Prepare {
    input {
        File? reference
        AlignedSample sample
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        Array[File] prep_dir = glob("portcullis_prep/*")
    }

    command <<<
        set -euxo pipefail
        portcullis prep -c -o portcullis_prep -t ~{cpus} ~{reference} ~{sep=" "sample.bam}
    >>>

   RuntimeAttr default_attr = object {
        cpu_cores: "~{cpus}",
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

task Junction {
    input {
        Array[File] prep_dir
        String strand = "firststrand"
        RuntimeAttr? runtime_attr_override
    }
    
    Int cpus = 8

    output {
        Array[File] junc_dir = glob("portcullis_junc/*")
        File tab = "portcullis_junc/portcullis.junctions.tab"
    }

    command <<<
        set -euxo pipefail
        strandness="UNKNOWN"
        case "~{strand}" in
            fr-firststrand)
            strandness="firststrand"
            ;;
            fr-secondstrand)
            strandness="secondstrand"
            ;;
            fr-unstranded)
            strandness="unstranded"
            ;;
        esac
        prep_dir_path="$(dirname ~{prep_dir[0]})"
        portcullis junc -c --strandedness="${strandness}" -t ~{cpus}  "${prep_dir_path}"
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

task Filter {
    input {
        Array[File] prep_dir
        Array[File] junc_dir
        File? reference_bed
        File tab
        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    
    output {
        File pass = "portcullis_filter.pass.junctions.tab"
    }

    command <<<
        set -euxo pipefail
        prep_dir_path="$(dirname ~{prep_dir[0]})"

        portcullis filter -o portcullis_filter --canonical=OFF \
        --max_length=2000 ~{"--reference " + reference_bed } \
        --threads=~{cpus} "${prep_dir_path}" ~{tab}
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
        Array[File] tabs
        RuntimeAttr? runtime_attr_override
    }
    
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

    output {
        File tab = "portcullis.merged.tab"
        File bed = "portcullis.merged.bed"
        File gff3 = "portcullis.merged.gff3"
    }

    command <<<
        set -euxo pipefail
        (junctools set --prefix=portcullis_merged --output=portcullis.merged.tab --operator=mean union ~{sep=" " tabs} || touch portcullis.merged.tab)
        junctools convert -if portcullis -of ebed --output=portcullis.merged.bed portcullis.merged.tab
        junctools convert -if portcullis -of igff --output=portcullis.merged.gff3 portcullis.merged.tab
    >>>
}