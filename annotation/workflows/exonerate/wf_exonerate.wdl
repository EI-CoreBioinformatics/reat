version 1.0

import "../common/structs.wdl"
import "../common/tasks.wdl" as tasks

workflow wf_exonerate {
    input {
        IndexedReference masked_reference_genome
        Array[LabeledFasta] related_species_protein
    }

    call ExonerateDatabase {
        input:
        target = masked_reference_genome
    }

    scatter (specie_proteins in related_species_protein) {
        call tasks.sanitizeFasta {
            input:
            reference = specie_proteins.fasta
        }

        call tasks.SplitSequences {
            input:
            prefix = specie_proteins.label,
            sequences_file = sanitizeFasta.sanitised_reference,
            num_out_files = 2
        }
        scatter (seq_file in SplitSequences.seq_files) {
            call Exonerate {
                input:
                query = seq_file,
                db_esi = ExonerateDatabase.esi,
                db_esd = ExonerateDatabase.esd
            }
        }

        Array[File] specie_exonerate_result = Exonerate.hits
    }

    Array[Array[File]] exonerate_hits = specie_exonerate_result

    output {
        Array[Array[File]] exonerate_results = exonerate_hits
    }
}


task ExonerateDatabase {
    input {
        IndexedReference target
    }

    output {
        File esi = basename(target.fasta)+".esi"
        File esd = basename(target.fasta)+".esd"
    }

    command <<<
    fasta2esd --softmask yes "~{target.fasta}" "~{basename(target.fasta)}.esd" && \
    esd2esi --translate yes "~{basename(target.fasta)}.esd" "~{basename(target.fasta)}.esi"
    >>>
}

task Exonerate {
    input {
        File db_esd
        File db_esi
        File query
    }

    output {
        File hits = "exonerate.hits"
    }

    command <<<
    ln ~{db_esd} .
    ln ~{db_esi} .
    exonerate_wrapper.py -M 2 -ir 20 2000 -t 4 --geneseed 250 --hspfilter 100 --score 50 --percent 30 \
    --serverlog server.log --log mapping.log ~{basename(db_esi)} ~{query} exonerate.hits
    >>>
}