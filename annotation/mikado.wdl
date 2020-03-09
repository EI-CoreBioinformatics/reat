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
    }

    parameter_meta {
        reference_genome: ""
        LQ_assemblies: ""
        HQ_assemblies: ""
        SR_assemblies: ""
        annotation_bed: ""
        mikado_scoring_file: ""
    }

    call mikado.wf_mikado as Mikado_long {
        input:
        scoring_file = mikado_scoring_file,
        indexed_reference =  reference_genome,
        SR_assemblies = SR_assemblies,
        LQ_assemblies = LQ_assemblies,
        HQ_assemblies = HQ_assemblies,
        junctions = annotation_bed,
        orf_calling_proteins = orf_calling_proteins,
        homology_proteins = homology_proteins
    }

    call mikado.wf_mikado as Mikado_short {
        input:
        scoring_file = mikado_scoring_file,
        indexed_reference =  reference_genome,
        SR_assemblies = SR_assemblies,
        junctions = annotation_bed,
        orf_calling_proteins = orf_calling_proteins,
        homology_proteins = homology_proteins
    }

    output {
        File mikado_long_config = Mikado_long.mikado_config
        File? mikado_long_orfs = Mikado_long.orfs

        File mikado_short_config = Mikado_short.mikado_config
        File? mikado_short_orfs = Mikado_short.orfs
    }
}