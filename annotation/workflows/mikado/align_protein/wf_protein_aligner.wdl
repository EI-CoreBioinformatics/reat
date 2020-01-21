version 1.0

task SanitiseProteinBlastDB {
    input {
        File db
    }

    output {
        File clean_db = "output.db"
    }

    command <<<
        sanitize_sequence_db.py -cstop ~{db} | gt seqtransform -addstopaminos -width "60" > "output.db"

    >>>
}

task BlastIndex {
    input {
        File target
    }

    output {
        File index = "blast_index.db"
    }

    command <<<
    blast index ~{target} > "blast_index.db"
    >>>
}

task BlastAlign {
    input {
        File index
        File query
    }

    output {
        File out = "output.xml"
    }

    command <<<
        blast ~{index} ~{query} > "output.txt"
    >>>
}

task DiamondIndex {
    input {
        File target
    }

    output {
        File index = "diamond_index.db"
    }

    command <<<
    diamond index ~{target} > "diamond_index.db"
    >>>
}

task DiamondAlign {
    input {
        File index
        File query
    }

    output {
        File out = "output.txt"
    }

    command <<<
        diamond blast ~{index} ~{query} > "output.xml"
    >>>
}
