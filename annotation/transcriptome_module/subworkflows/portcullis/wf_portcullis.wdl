version 1.0

import "../common/structs.wdl"
import "../common/rt_struct.wdl"

workflow portcullis {
    input {
        File reference
        String merge_operator = "max"
        Object? group_to_samples # Consider applying "localization_optional" to the preparation task
        Array[AlignedSample] aligned_samples
        File? annotation
        RuntimeAttr? portcullis_resources
    }

    parameter_meta {
        reference: {description: "Genome reference file.", category: "required"}
        annotation: {description:"Reference annotation of junctions in BED format.", category:"required"}
        sample_groups: {description:"Pre-defined groupings for the aligned_samples, these need to match the sample names. Each portcullis run groups one or more samples.", category:"required"}
        aligned_samples: {description:"List of aligned samples.", category:"required"}
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
        scatter (sample_files in PrepareGroups.groupfile_to_bams) {
            call Full as groupedFull{
                input:
                reference = reference,
                sample = object { name: basename(sample_files,'.sample'), strand: sub(basename(sample_files,'.sample'),'.+\_',''), aligner: "*", bam: read_lines(sample_files), merge:false},
                reference_bed = PrepareRef.refbed,
                runtime_attr_override = portcullis_resources
            }
        }
        Array[File] def_grouped_pass = groupedFull.pass_tab
        Array[File] def_grouped_fail = groupedFull.pass_tab
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
        Array[File] def_ungrouped_pass = Full.pass_tab
        Array[File] def_ungrouped_fail = Full.pass_tab
    }

    Array[File] to_merge_pass = select_first([def_grouped_pass, def_ungrouped_pass])
    Array[File] to_merge_fail = select_first([def_grouped_fail, def_ungrouped_fail])

    call Merge as PassMerge {
        input:
        merge_operator = merge_operator,
        tabs = to_merge_pass,
        output_directory = "portcullis",
        is_pass = true
    }

    call Merge as FailMerge {
        input:
        merge_operator = merge_operator,
        tabs = to_merge_fail,
        output_directory = "portcullis",
        is_pass = false
    }

    output {
        File pass_tab = PassMerge.tab
        File pass_bed = PassMerge.bed
        File pass_gff3 = PassMerge.gff3

        File fail_tab = FailMerge.tab
        File fail_bed = FailMerge.bed
        File fail_gff3 = FailMerge.gff3
    }
}

task PrepareGroups {
    input {
        Object group_to_samples
        Array[AlignedSample] aligned_samples
    }

    output {
        Array[File] groupfile_to_bams = glob('*.sample')
    }

    command <<<
        portcullis_sample_grouping ~{write_objects(aligned_samples)} ~{write_json(group_to_samples)}
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
        File pass_tab = "portcullis_out/3-filt/portcullis_filtered.pass.junctions.tab"

        File fail_tab = "portcullis_out/3-filt/portcullis_filtered.fail.junctions.tab"
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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 32,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

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
        String merge_operator
        Boolean is_pass
        String output_directory
        RuntimeAttr? runtime_attr_override
    }
    
    String type_text = if(is_pass) then "pass" else "fail"
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
        File tab = output_directory + "/portcullis." + type_text + ".merged.tab"
        File bed = output_directory + "/portcullis." + type_text + ".merged.bed"
        File gff3 = output_directory + "/portcullis." + type_text + ".merged.gff3"
    }

    command <<<
        set -euxo pipefail
        mkdir ~{output_directory}
        cd ~{output_directory}
        (junctools set --prefix=portcullis_merged --output=portcullis.~{type_text}.merged.tab --operator=~{merge_operator} union ~{sep=" " tabs} || touch portcullis.~{type_text}.merged.tab)
        junctools convert -if portcullis -of ebed --output=portcullis.~{type_text}.merged.bed portcullis.~{type_text}.merged.tab
        junctools convert -if portcullis -of igff --output=portcullis.~{type_text}.merged.gff3 portcullis.~{type_text}.merged.tab
    >>>
}