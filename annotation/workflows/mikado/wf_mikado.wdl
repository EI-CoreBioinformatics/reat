version 1.0

import "../structs/structs.wdl"

workflow wf_mikado {
    input {
        Array[AssembledSample] assemblies
        Array[AssembledSample]? long_assemblies
        File reference_fasta
        String gencode
        String orf_caller
        File? scoring_file
    }

    scatter (sr_assembly in assemblies) {
        call GenerateModelsList as sr_models {
            input:
            assembly = sr_assembly
        }
    }

    Array[AssembledSample] def_long_assemblies = select_first([long_assemblies])

    scatter (lr_assembly in def_long_assemblies) {
        call GenerateModelsList as lr_models {
            input:
            assembly = lr_assembly
        }
    }

    File result = write_lines(flatten([sr_models.models, lr_models.models]))

    call MikadoConfigure {
        input:
        reference_fasta = reference_fasta,
        models = result,
        scoring_file = scoring_file
    }

    if (orf_caller == "Prodigal") {
        call Prodigal {
            input:
            gencode = gencode,
            reference = reference_fasta
        }
    }

    if (orf_caller == "GTCDS") {
        call GTCDS {
            input:
            reference = reference_fasta,
            gtf = MikadoConfigure.prepared_gtf
        }
    }

    File def_orfs = select_first([Prodigal.orfs, GTCDS.orfs])

    output {
        File mikado_config = MikadoConfigure.mikado_config
        File orfs = def_orfs
    }

}

task GTCDS {
    input {
        File reference
        File gtf
    }

    output {
        File orfs = "mikado_prepared.gt_cds.trans.bed12"
        File gff3 = "mikado_prepared.gt_cds.gff3"
    }

    command <<<
        awk '$3!~\"(CDS|UTR)\"' ~{gtf} \
        | mikado util convert -if gtf -of gff3 - \
        | gt gff3 -tidy -retainids -addids | gt cds -seqfile ~{reference} - matchdesc \
        | gff3_name_to_id.py - mikado_prepared.gt_cds.gff3 && \
        mikado util convert -t -of bed12 "mikado_prepared.gt_cds.gff3" "mikado_prepared.gt_cds.trans.bed12"
    >>>
}

task Prodigal {
    input {
        File reference
        String gencode
    }

    output {
        File orfs = "transcripts.fasta.prodigal.gff3"
    }

    command <<<
        code_id=$(python -c "import Bio.Data; print(CodonTable.generic_by_name[~{gencode}].id")
        prodigal -f gff -g "${code_id}" -i "~{reference}" -o "transcripts.fasta.prodigal.gff3"
    >>>
}

task GenerateModelsList {
    input {
        AssembledSample assembly
    }

    output {
        String models = read_string(stdout())
    }

    command <<<
        echo -e "~{assembly.assembly}\t~{assembly.name}\tTrue"
    >>>
}

task MikadoConfigure {
    input {
        File models
        File reference_fasta
        File? scoring_file
    }

    output {
        File mikado_config = "mikado.yaml"
        File prepared_fasta = "mikado_prepare/mikado_prepared.fasta"
        File prepared_gtf = "mikado_prepare/mikado_prepared.gtf"
    }

    command <<<
        mikado configure \
        ~{"--scoring=" + scoring_file} \
        --list=~{models} \
        ~{"--reference=" + reference_fasta} \
        mikado.yaml

        mikado prepare --procs=4 --json-conf=mikado.yaml -od mikado_prepare --strip_cds
    >>>
}