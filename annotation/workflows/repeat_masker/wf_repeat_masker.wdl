version 1.0

workflow wf_repeat_masker {
    input {
        File reference_fasta
        File repeat_db
        Boolean run_modeller
        Boolean retrieve_known
        Array[File]? safe_cds_sequences
        String? clade
        String? specie
    }

    if (run_modeller) {
        call BuildModellerDB {
            input:
            reference_fasta = reference_fasta
        }

        call RepeatModeller {
            input:
            db = BuildModellerDB.db
        }
    }

    if (retrieve_known) {
        String def_one = select_first([clade, specie])
        call RetrieveLibraries {
            input:
            clade = clade,
            specie = specie,
            any = def_one
        }
    }


}

task BuildModellerDB {
    input {
        File reference_fasta
    }

    output {
        Array[File] db = glob("genome.*")
    }

    command <<<
    BuildDatabase -name genome -engine ncbi ~{reference_fasta}
    >>>
}

task RepeatModeller {
    input {
        Array[File] db
    }

    output {
        File families = "families.txt"
    }

    command <<<
    ln -s ~{sep=" " db} -t .
    RepeatModeler -engine ncbi -pa 4 -database genome
    >>>
}

task RetrieveLibraries {
    input {
        String? clade
        String? specie
        String any
    }

    output {
        File retrieved = "retrieved.fa"
    }

    command <<<
    queryRepeatDatabase.pl ~{"-species " + specie} ~{"-clade " + clade} > "retrieved.fa"
    >>>
}