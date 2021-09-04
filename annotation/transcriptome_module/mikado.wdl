version 1.0

import "subworkflows/common/structs.wdl"
import "subworkflows/common/tasks.wdl"
import "subworkflows/mikado/wf_mikado.wdl" as mikado

workflow wf_main_mikado {
    input {
        IndexedReference reference_genome
        File? annotation

        Int annotation_score = 1
        String mode
        String genetic_code
        Int mikado_genetic_code
        Boolean check_reference

        File all_scoring_file
        File long_scoring_file
        File? long_lq_scoring_file
        Array[AssembledSample]? LQ_assemblies
        Array[AssembledSample]? HQ_assemblies
        Array[AssembledSample]? SR_assemblies
        File? junctions_bed
        File? HQ_junctions_bed
        File? LQ_junctions_bed
        File? homology_proteins
        File? orf_calling_proteins
        Boolean separate_LQ = false
        Boolean exclude_LQ_junctions = false
        Boolean skip_mikado_long = false
        Boolean filter_HQ_assemblies = false
        Boolean filter_LQ_assemblies = false
        String? orf_calling_program

        File? all_prepare_cfg
        File? all_serialise_cfg
        File? all_pick_cfg
        File? long_prepare_cfg
        File? long_serialise_cfg
        File? long_pick_cfg
        File? long_lq_prepare_cfg
        File? long_lq_serialise_cfg
        File? long_lq_pick_cfg

        RuntimeAttr? orf_calling_resources
        RuntimeAttr? protein_index_resources
        RuntimeAttr? protein_alignment_resources
        RuntimeAttr? homology_index_resources
        RuntimeAttr? homology_alignment_resources

        RuntimeAttr? mikado_all_pick_resources
        RuntimeAttr? mikado_all_serialise_resources
        RuntimeAttr? mikado_all_prepare_resources

        RuntimeAttr? mikado_long_pick_resources
        RuntimeAttr? mikado_long_serialise_resources
        RuntimeAttr? mikado_long_prepare_resources

        RuntimeAttr? mikado_long_lq_pick_resources
        RuntimeAttr? mikado_long_lq_serialise_resources
        RuntimeAttr? mikado_long_lq_prepare_resources
    }

    parameter_meta {
        reference_genome: "Reference genome to align against"
        paired_samples: "Paired short read samples, each item is defined by a biological replicate name with one or more technical replicates. Technical replicates are defined by a name, R1, R2 and strand."
        LQ_long_read_samples: "Low quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        HQ_long_read_samples: "High quality long read samples, each item is defined by a name, it's strand and one or more long read files."
        mikado_scoring_file: "Mikado scoring file"
    }

    RuntimeAttr default_runtime_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1,
        queue: ""
    }

    Boolean run_mikado_homology = defined(homology_proteins)

    # combine junctions
    if (!exclude_LQ_junctions) {
        call CombineAllJunctions as all {
            input:
            portcullis_junctions = junctions_bed,
            HQ_junctions = HQ_junctions_bed,
            LQ_junctions = LQ_junctions_bed
        }
    }

    if (exclude_LQ_junctions) {
        call CombineAllJunctions as noLQ {
            input:
            portcullis_junctions = junctions_bed,
            HQ_junctions = HQ_junctions_bed
        }
    }

    File def_junctions = select_first([all.merged_junctions, noLQ.merged_junctions])

    if (filter_HQ_assemblies && defined(HQ_assemblies)) {
        scatter (sample in select_first([HQ_assemblies])){
            call FilterWithJunctions as HQ_filtered {
                input:
                bed_junctions = def_junctions,
                gtf = sample.assembly,
                name = sample.name
            }
            AssembledSample hq_filtered = object {
                                                 name: sample.name+".filtered",
                                                 strand: sample.strand,
                                                 assembly: HQ_filtered.filtered_gtf,
                                                 score: sample.score,
                                                 is_ref: sample.is_ref,
                                                 exclude_redundant: sample.exclude_redundant
                                             }
        }
    }

    if (filter_LQ_assemblies && defined(LQ_assemblies)) {
        scatter (sample in select_first([LQ_assemblies])){
            call FilterWithJunctions as LQ_filtered {
                input:
                bed_junctions = def_junctions,
                gtf = sample.assembly,
                name = sample.name
            }
            AssembledSample lq_filtered = object {
                                                 name: sample.name+".filtered",
                                                 strand: sample.strand,
                                                 assembly: LQ_filtered.filtered_gtf,
                                                 score: sample.score,
                                                 is_ref: sample.is_ref,
                                                 exclude_redundant: sample.exclude_redundant
                                             }
        }
    }

    if (defined(HQ_assemblies)) {
        Array[AssembledSample] def_hq_assemblies = select_first([hq_filtered, HQ_assemblies])
    }

    if (defined(LQ_assemblies)) {
        Array[AssembledSample] def_lq_assemblies = select_first([lq_filtered, LQ_assemblies])
    }

    # The user can choose to run the LQ-LR datasets separately
    if (separate_LQ)
    {
        call mikado.wf_mikado as Mikado_short_and_long_noLQ {
            input:
            annotation = annotation,
            annotation_score = annotation_score,
            mode = mode,
            check_reference = check_reference,
            indexed_reference =  reference_genome,
            SR_assemblies = SR_assemblies,
            HQ_assemblies = def_hq_assemblies,
            scoring_file = all_scoring_file,
            orf_calling_proteins = orf_calling_proteins,
            orf_caller = orf_calling_program,
            genetic_code = genetic_code,
            mikado_genetic_code = mikado_genetic_code,
            mikado_do_homology_assessment = run_mikado_homology,
            homology_proteins = homology_proteins,
            junctions = def_junctions,
            output_prefix = "mikado_all_noLQ",
            prepare_extra_config = all_prepare_cfg,
            serialise_extra_config = all_serialise_cfg,
            pick_extra_config = all_pick_cfg,
            mikado_prepare_resources = mikado_all_prepare_resources,
            mikado_serialise_resources = mikado_all_serialise_resources,
            mikado_pick_resources = mikado_all_pick_resources,
            orf_calling_resources = select_first([orf_calling_resources,default_runtime_attr]),
            orf_protein_index_resources = select_first([protein_index_resources,default_runtime_attr]),
            orf_protein_alignment_resources = select_first([protein_alignment_resources,default_runtime_attr]),
            homology_index_resources = select_first([homology_index_resources,default_runtime_attr]),
            homology_alignment_resources = select_first([homology_alignment_resources,default_runtime_attr])
        }

        if (defined(HQ_assemblies)) {
            call mikado.wf_mikado as Mikado_longHQ {
                input:
                annotation = annotation,
                annotation_score = annotation_score,
                mode = mode,
                check_reference = check_reference,
                indexed_reference =  reference_genome,
                HQ_assemblies = def_hq_assemblies,
                scoring_file = long_scoring_file,
                orf_calling_proteins = orf_calling_proteins,
                orf_caller = orf_calling_program,
                genetic_code = genetic_code,
                mikado_genetic_code = mikado_genetic_code,
                mikado_do_homology_assessment = run_mikado_homology,
                homology_proteins = homology_proteins,
                junctions = def_junctions,
                output_prefix = "mikado_longHQ",
                prepare_extra_config = long_prepare_cfg,
                serialise_extra_config = long_serialise_cfg,
                pick_extra_config = long_pick_cfg,
                mikado_prepare_resources = mikado_long_prepare_resources,
                mikado_serialise_resources = mikado_long_serialise_resources,
                mikado_pick_resources = mikado_long_pick_resources,
                orf_calling_resources = select_first([orf_calling_resources,default_runtime_attr]),
                orf_protein_index_resources = select_first([protein_index_resources,default_runtime_attr]),
                orf_protein_alignment_resources = select_first([protein_alignment_resources,default_runtime_attr]),
                homology_index_resources = select_first([homology_index_resources,default_runtime_attr]),
                homology_alignment_resources = select_first([homology_alignment_resources,default_runtime_attr])
            }
        }

        if (defined(LQ_assemblies)) {
            call mikado.wf_mikado as Mikado_longLQ {
                input:
                annotation = annotation,
                annotation_score = annotation_score,
                mode = mode,
                check_reference = check_reference,
                scoring_file = select_first([long_lq_scoring_file]),
                indexed_reference =  reference_genome,
                LQ_assemblies = def_lq_assemblies,
                junctions = def_junctions,
                orf_calling_proteins = orf_calling_proteins,
                orf_caller = orf_calling_program,
                genetic_code = genetic_code,
                mikado_genetic_code = mikado_genetic_code,
                mikado_do_homology_assessment = run_mikado_homology,
                homology_proteins = homology_proteins,
                output_prefix = "mikado_longLQ",
                prepare_extra_config = long_lq_prepare_cfg,
                serialise_extra_config = long_lq_serialise_cfg,
                pick_extra_config = long_lq_pick_cfg,
                mikado_prepare_resources = mikado_long_lq_prepare_resources,
                mikado_serialise_resources = mikado_long_lq_serialise_resources,
                mikado_pick_resources = mikado_long_lq_pick_resources,
                orf_calling_resources = select_first([orf_calling_resources,default_runtime_attr]),
                orf_protein_index_resources = select_first([protein_index_resources,default_runtime_attr]),
                orf_protein_alignment_resources = select_first([protein_alignment_resources,default_runtime_attr]),
                homology_index_resources = select_first([homology_index_resources,default_runtime_attr]),
                homology_alignment_resources = select_first([homology_alignment_resources,default_runtime_attr])
            }
        }
    }
    if (!separate_LQ)
    {
        call mikado.wf_mikado as Mikado_short_and_long {
            input:
            annotation = annotation,
            annotation_score = annotation_score,
            mode = mode,
            check_reference = check_reference,
            scoring_file = all_scoring_file,
            indexed_reference =  reference_genome,
            SR_assemblies = SR_assemblies,
            LQ_assemblies = def_lq_assemblies,
            HQ_assemblies = def_hq_assemblies,
            junctions = def_junctions,
            orf_calling_proteins = orf_calling_proteins,
            orf_caller = orf_calling_program,
            genetic_code = genetic_code,
            mikado_genetic_code = mikado_genetic_code,
            mikado_do_homology_assessment = run_mikado_homology,
            homology_proteins = homology_proteins,
            output_prefix = "mikado_all",
            prepare_extra_config = all_prepare_cfg,
            serialise_extra_config = all_serialise_cfg,
            pick_extra_config = all_pick_cfg,
            mikado_prepare_resources = mikado_all_prepare_resources,
            mikado_serialise_resources = mikado_all_serialise_resources,
            mikado_pick_resources = mikado_all_pick_resources,
            orf_calling_resources = select_first([orf_calling_resources, default_runtime_attr]),
            orf_protein_index_resources = select_first([protein_index_resources, default_runtime_attr]),
            orf_protein_alignment_resources = select_first([protein_alignment_resources, default_runtime_attr]),
            homology_index_resources = select_first([homology_index_resources, default_runtime_attr]),
            homology_alignment_resources = select_first([homology_alignment_resources, default_runtime_attr])
        }

        if (!skip_mikado_long && (defined(LQ_assemblies) || defined(HQ_assemblies))) {
            call mikado.wf_mikado as Mikado_long {
                input:
                annotation = annotation,
                annotation_score = annotation_score,
                mode = mode,
                check_reference = check_reference,
                scoring_file = long_scoring_file,
                indexed_reference =  reference_genome,
                LQ_assemblies = def_lq_assemblies,
                HQ_assemblies = def_hq_assemblies,
                junctions = def_junctions,
                orf_calling_proteins = orf_calling_proteins,
                orf_caller = orf_calling_program,
                genetic_code = genetic_code,
                mikado_genetic_code = mikado_genetic_code,
                mikado_do_homology_assessment = run_mikado_homology,
                homology_proteins = homology_proteins,
                output_prefix = "mikado_long",
                prepare_extra_config = long_prepare_cfg,
                serialise_extra_config = long_serialise_cfg,
                pick_extra_config = long_pick_cfg,
                mikado_prepare_resources = mikado_long_prepare_resources,
                mikado_serialise_resources = mikado_long_serialise_resources,
                mikado_pick_resources = mikado_long_pick_resources,
                orf_calling_resources = select_first([orf_calling_resources, default_runtime_attr]),
                orf_protein_index_resources = select_first([protein_index_resources, default_runtime_attr]),
                orf_protein_alignment_resources = select_first([protein_alignment_resources, default_runtime_attr]),
                homology_index_resources = select_first([homology_index_resources, default_runtime_attr]),
                homology_alignment_resources = select_first([homology_alignment_resources, default_runtime_attr])
            }
        }
    }

    Array[File] mikado_stats_files = select_all([
            Mikado_long.stats, 
        Mikado_short_and_long.stats, 
        Mikado_longHQ.stats, 
        Mikado_longLQ.stats, 
        Mikado_short_and_long_noLQ.stats
        ])

    call tasks.TranscriptAssemblySummaryStats {
        input:
        stats = mikado_stats_files,
        output_prefix = "mikado"
    }

    output {
        File? long_orfs = Mikado_long.orfs
        File? long_loci = Mikado_long.loci
        File? long_scores = Mikado_long.scores
        File? long_metrics = Mikado_long.metrics
        File? long_stats = Mikado_long.stats

        File? short_orfs = Mikado_short_and_long.orfs
        File? short_loci = Mikado_short_and_long.loci
        File? short_scores = Mikado_short_and_long.scores
        File? short_metrics = Mikado_short_and_long.metrics
        File? short_stats = Mikado_short_and_long.stats

        File? short_and_long_noLQ_orfs = Mikado_short_and_long_noLQ.orfs
        File? short_and_long_noLQ_loci = Mikado_short_and_long_noLQ.loci
        File? short_and_long_noLQ_scores = Mikado_short_and_long_noLQ.scores
        File? short_and_long_noLQ_metrics = Mikado_short_and_long_noLQ.metrics
        File? short_and_long_noLQ_stats = Mikado_short_and_long_noLQ.stats

        File? longHQ_orfs = Mikado_longHQ.orfs
        File? longHQ_loci = Mikado_longHQ.loci
        File? longHQ_scores = Mikado_longHQ.scores
        File? longHQ_metrics = Mikado_longHQ.metrics
        File? longHQ_stats = Mikado_longHQ.stats

        File? longLQ_orfs = Mikado_longLQ.orfs
        File? longLQ_loci = Mikado_longLQ.loci
        File? longLQ_scores = Mikado_longLQ.scores
        File? longLQ_metrics = Mikado_longLQ.metrics
        File? longLQ_stats = Mikado_longLQ.stats

        File mikado_stats_summary = TranscriptAssemblySummaryStats.summary
    }
}

task CombineAllJunctions {
	input {
		File? portcullis_junctions
		File? HQ_junctions
		File? LQ_junctions
	}

	output {
		File merged_junctions = "all_junctions.bed"
	}

	command <<<
		cat ~{portcullis_junctions} ~{HQ_junctions} ~{LQ_junctions} > merged_junctions.bed
		junctools convert -if bed -of ebed -s -d merged_junctions.bed -o all_junctions.bed
	>>>
}

task FilterWithJunctions {
    input {
        String name
        File gtf
        File bed_junctions
    }

    output {
        File filtered_gtf = name + ".filtered.gtf"
    }

    command <<<
        junctools gtf filter -o ~{name}.filtered.gtf -j ~{bed_junctions} ~{gtf}
    >>>
}