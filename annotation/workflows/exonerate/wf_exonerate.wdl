version 1.0

import "../structs/structs.wdl"
import "../structs/tasks.wdl" as tasks

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
            sequences_file = sanitizeFasta.sanitised_reference
        }

        scatter (protein_chunk in SplitSequences.seq_files) {
            call Exonerate {
                input:
                query = protein_chunk,
                target_db = ExonerateDatabase.db
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
        File db = basename(target.fasta)+".esi"
    }

    command <<<
    fasta2esd --softmask yes ~{target.fasta} ~{basename(target.fasta)}".esd"
    esd2esi --translate yes ~{basename(target.fasta)}".esd" ~{basename(target.fasta)}".esi"
    >>>
}

task Exonerate {
    input {
        File target_db
        File query
    }

    output {
        File hits = basename(query)+".hits"
    }

    command <<<
    exonerate_wrapper.py -ir 20 2000 -t 4 --geneseed 250 --hspfilter 100 --score 50 --percent 30 \
    --serverlog server.log --log mapping.log ~{target_db} ~{basename(query)}.hits
    >>>
}