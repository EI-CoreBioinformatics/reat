version 1.0

import "workflows/common/structs.wdl"
import "workflows/mikado/wf_mikado.wdl" as mikado

workflow wf_main_mikado {
    input {
        IndexedReference reference_genome
        File mikado_scoring_file
        Array[AssembledSample]? LQ_assemblies
        Array[AssembledSample]? HQ_assemblies
        Array[AssembledSample]? SR_assemblies
        File? annotation_bed
        File homology_proteins
        File orf_calling_proteins
        Boolean separate_LQ = false
        RuntimeAttr orf_calling_resources
        RuntimeAttr protein_index_resources
        RuntimeAttr protein_alignment_resources
        RuntimeAttr homology_index_resources
        RuntimeAttr homology_alignment_resources
    }

    parameter_meta {
        reference_genome: ""
        LQ_assemblies: ""
        HQ_assemblies: ""
        SR_assemblies: ""
        annotation_bed: ""
        mikado_scoring_file: ""
    }

    # The user can choose to run the LQ-LR datasets separately
    if (separate_LQ)
    {
        call mikado.wf_mikado as Mikado_short_and_long_noLQ {
            input:
            indexed_reference =  reference_genome,
            SR_assemblies = SR_assemblies,
            HQ_assemblies = HQ_assemblies,
            scoring_file = mikado_scoring_file,
            orf_calling_proteins = orf_calling_proteins,
            homology_proteins = homology_proteins,
            junctions = annotation_bed,
            orf_calling_resources = orf_calling_resources,
            orf_protein_index_resources = protein_index_resources,
            orf_protein_alignment_resources = protein_alignment_resources,
            homology_index_resources = homology_index_resources,
            homology_alignment_resources = homology_alignment_resources
        }

        call mikado.wf_mikado as Mikado_longHQ {
            input:
            indexed_reference =  reference_genome,
            HQ_assemblies = HQ_assemblies,
            scoring_file = mikado_scoring_file,
            orf_calling_proteins = orf_calling_proteins,
            homology_proteins = homology_proteins,
            junctions = annotation_bed,
            orf_calling_resources = orf_calling_resources,
            orf_protein_index_resources = protein_index_resources,
            orf_protein_alignment_resources = protein_alignment_resources,
            homology_index_resources = homology_index_resources,
            homology_alignment_resources = homology_alignment_resources
        }

        call mikado.wf_mikado as Mikado_longLQ {
            input:
            scoring_file = mikado_scoring_file,
            indexed_reference =  reference_genome,
            LQ_assemblies = LQ_assemblies,
            junctions = annotation_bed,
            orf_calling_proteins = orf_calling_proteins,
            homology_proteins = homology_proteins,
            orf_calling_resources = orf_calling_resources,
            orf_protein_index_resources = protein_index_resources,
            orf_protein_alignment_resources = protein_alignment_resources,
            homology_index_resources = homology_index_resources,
            homology_alignment_resources = homology_alignment_resources
        }
    }
    if (!separate_LQ)
    {
        call mikado.wf_mikado as Mikado_short_and_long {
            input:
            scoring_file = mikado_scoring_file,
            indexed_reference =  reference_genome,
            SR_assemblies = SR_assemblies,
            LQ_assemblies = LQ_assemblies,
            HQ_assemblies = HQ_assemblies,
            junctions = annotation_bed,
            orf_calling_proteins = orf_calling_proteins,
            homology_proteins = homology_proteins,
            orf_calling_resources = orf_calling_resources,
            orf_protein_index_resources = protein_index_resources,
            orf_protein_alignment_resources = protein_alignment_resources,
            homology_index_resources = homology_index_resources,
            homology_alignment_resources = homology_alignment_resources
        }

        call mikado.wf_mikado as Mikado_long {
            input:
            scoring_file = mikado_scoring_file,
            indexed_reference =  reference_genome,
            LQ_assemblies = LQ_assemblies,
            HQ_assemblies = HQ_assemblies,
            junctions = annotation_bed,
            orf_calling_proteins = orf_calling_proteins,
            homology_proteins = homology_proteins,
            orf_calling_resources = orf_calling_resources,
            orf_protein_index_resources = protein_index_resources,
            orf_protein_alignment_resources = protein_alignment_resources,
            homology_index_resources = homology_index_resources,
            homology_alignment_resources = homology_alignment_resources
        }
    }

    output {
        File? long_config = Mikado_long.mikado_config
        File? long_orfs = Mikado_long.orfs

        File? short_config = Mikado_short_and_long.mikado_config
        File? short_orfs = Mikado_short_and_long.orfs

        File? short_and_long_noLQ_config = Mikado_short_and_long_noLQ.mikado_config
        File? short_and_long_noLQ_orfs = Mikado_short_and_long_noLQ.orfs

        File? longHQ_config = Mikado_longHQ.mikado_config
        File? longHQ_orfs = Mikado_longHQ.orfs

        File? longLQ_config = Mikado_longLQ.mikado_config
        File? longLQ_orfs = Mikado_longLQ.orfs

    }
}