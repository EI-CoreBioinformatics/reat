version 1.0

import "../common/structs.wdl"
import "orf_caller/wf_transdecoder.wdl" as tdc
import "homology/wf_homology.wdl" as hml

workflow wf_mikado {
    input {
        Array[AssembledSample]? short_assemblies
        Array[AssembledSample]? long_assemblies
        IndexedReference indexed_reference
        File? junctions
        String gencode = "Universal"
        String orf_caller = "None"
        Boolean mikado_do_homology_assessment = false
        File? scoring_file
        File? dbs
    }

    if (defined(short_assemblies)) {
        Array[AssembledSample] def_short_assemblies = select_first([short_assemblies])
        scatter (sr_assembly in def_short_assemblies) {
            call GenerateModelsList as sr_models {
                input:
                assembly = sr_assembly
            }
        }
    }

    if (defined(long_assemblies)) {
        Array[AssembledSample] def_long_assemblies = select_first([long_assemblies])

        scatter (lr_assembly in def_long_assemblies) {
            call GenerateModelsList as lr_models {
                input:
                assembly = lr_assembly,
                long_score_bias = 1
            }
        }
    }

    File result = write_lines(flatten(select_all([sr_models.models, lr_models.models])))

    call MikadoPrepare {
        input:
        reference_fasta = indexed_reference.fasta,
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
        indexed_reference = indexed_reference,
        config = MikadoPrepare.mikado_config
    }

    call MikadoPick {
        input:
        config_file = MikadoPrepare.mikado_config,
        mikado_db = MikadoSerialise.mikado_db,
        transcripts = MikadoPrepare.prepared_gtf
    }

    output {
        File mikado_config = MikadoPrepare.mikado_config
        File? orfs = maybe_orfs
        Array[File]? homologies = Homology.homology
        Array[File] serialise_out = MikadoSerialise.out
        File pick_metrics = MikadoPick.metrics
        File pick_loci = MikadoPick.loci
        File pick_scores = MikadoPick.scores
        File pick_stats = MikadoPick.stats
    }
}

task MikadoPick {
#modes = ("permissive", "stringent", "nosplit", "split", "lenient")
    input {
        File config_file
        File transcripts
        File mikado_db
        String mode = "permissive"
        Int flank = 200
    }

    output {
        File metrics = "mikado_pick/mikado-" + mode + ".loci.metrics.tsv"
        File loci = "mikado_pick/mikado-" + mode + ".loci.gff3"
        File scores = "mikado_pick/mikado-" + mode + ".loci.scores.tsv"
        File loci_index = "mikado_pick/mikado-"+mode+".loci.gff3.midx"
        File index_log = "mikado_pick/index_loci.log"
        File stats = "mikado_pick/mikado-" + mode + ".loci.gff3.stats"
    }

    command <<<
    export TMPDIR=/tmp
    mkdir -p mikado_pick
    mikado pick ~{"--source Mikado_" + mode} ~{"--mode " + mode} --procs=4 \
    ~{"--flank " + flank} --start-method=spawn ~{"--json-conf=" + config_file} \
    -od mikado_pick --loci-out mikado-~{mode}.loci.gff3 -lv INFO ~{"-db " + mikado_db} \
    ~{transcripts}
    mikado compare -r mikado_pick/mikado-~{mode}.loci.gff3 -l mikado_pick/index_loci.log --index
    mikado util stats  mikado_pick/mikado-~{mode}.loci.gff3 mikado_pick/mikado-~{mode}.loci.gff3.stats
    >>>
}

task MikadoSerialise {
    input {
        IndexedReference indexed_reference
        File config
        File transcripts
        Array[File]? homology_alignments
        File? clean_seqs_db
        File? junctions
        File? orfs
    }

    output {
        Array[File] out = glob("mikado_serialise/*")
        File mikado_db = "mikado_serialise/mikado.db"
    }

    String xml_prefix = if defined(homology_alignments) then "--xml=" else ""

    command <<<
    fasta=~{indexed_reference.fasta}
    fai=~{indexed_reference.fai}
    
    ln -s ${fasta} .
    ln -s ${fai} .
    mikado serialise ~{xml_prefix}~{sep="," homology_alignments} ~{"--blast_targets="+clean_seqs_db} ~{"--junctions="+junctions} ~{"--orfs="+orfs} \
    ~{"--transcripts=" + transcripts} --genome_fai=${fai} \
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
        Int long_score_bias = 0
    }

    output {
        String models = read_string(stdout())
    }

    command <<<
        strand="False"
        if [ ~{assembly.strand} != "fr-unstranded" ]; then
            strand="True"
        fi
        echo -e "~{assembly.assembly}\t~{assembly.name}\t${strand}\t~{long_score_bias}"
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
