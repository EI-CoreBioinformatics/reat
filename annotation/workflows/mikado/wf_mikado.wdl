version 1.0

import "../structs/structs.wdl"
import "orf_caller/wf_transdecoder.wdl" as tdc
import "homology/wf_homology.wdl" as hml

workflow wf_mikado {
    input {
        Array[AssembledSample] assemblies
        Array[AssembledSample]? long_assemblies
        File reference_fasta
        File junctions
        String gencode = "Universal"
        String orf_caller = "None"
        Boolean mikado_do_homology_assessment = false
        File? scoring_file
        File? dbs
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

    call MikadoPrepare {
        input:
        reference_fasta = reference_fasta,
        models = result,
        scoring_file = scoring_file
    }

    # ORF Calling
    if (orf_caller != "None") {
        if (orf_caller == "Prodigal") {
            call Prodigal {
                input:
                gencode = gencode,
                reference = MikadoPrepare.prepared_fasta
            }
        }

        if (orf_caller == "GTCDS") {
            call GTCDS {
                input:
                reference = MikadoPrepare.prepared_fasta,
                gtf = MikadoPrepare.prepared_gtf
            }
        }

        if (orf_caller == "Transdecoder") {
            call tdc.wf_transdecoder as Transdecoder {
                input:
                prepared_fa = MikadoPrepare.prepared_fasta
            }
        }

        File maybe_orfs = select_first([Prodigal.orfs, GTCDS.orfs, Transdecoder.orfs])
    }

    # Mikado Homology
    if (defined(mikado_do_homology_assessment)) {
        call hml.wf_homology as Homology {
            input:
            program = "blastx",
            reference = MikadoPrepare.prepared_fasta,
            protein_db = dbs
        }
    }

    call MikadoSerialise {
        input:
        homology_alignments = Homology.homology,
        clean_seqs_db = Homology.homology_clean_db,
        junctions = junctions,
        orfs = maybe_orfs,
        transcripts = MikadoPrepare.prepared_fasta,
        indexed_reference = reference_fasta,
        config = MikadoPrepare.mikado_config
    }

    output {
        File mikado_config = MikadoPrepare.mikado_config
        File? orfs = maybe_orfs
        Array[File]? homologies = Homology.homology
        Array[File] serialise_out = MikadoSerialise.out
    }
}

task MikadoSerialise {
    input {
        Array[File]? homology_alignments
        File? clean_seqs_db
        File junctions
        File? orfs
        File transcripts
        File indexed_reference
        File config
    }

    output {
        Array[File] out = glob("mikado_serialise/*")
    }

    command <<<
    mikado serialise ~{sep=" " homology_alignments} ~{clean_seqs_db} ~{junctions} ~{orfs} 
    ~{"--transcripts=" + transcripts} ~{"--genome_fai="+indexed_reference}
    ~{"--json-conf=" + config} --force --start-method=spawn -od mikado_serialise --procs=4
    >>>
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
        strand="False"
        if [ ~{assembly.strand} != "fr-unstranded" ]; then
            strand="True"
        fi
        echo -e "~{assembly.assembly}\t~{assembly.name}\t${strand}"
    >>>
}

task MikadoPrepare {
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
