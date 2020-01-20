version 1.0

import "../structs/structs.wdl"

workflow wf_transdecoder {
    input {
        File prepared_fa
        File prepared_gtf
        Array[File]? dbs
    }

    call TransdecoderLongOrf {
        input:
            reference = prepared_fa
    }

    if (defined(dbs)) {
        call SanitiseProteinBlastDB { # Check if this is an input for both homology + orf_calling or just the same task
            input:
                dbs = dbs
        }
    }


    output {
        File long_orfs = TransdecoderLongOrf.orfs
        File? clean_db = SanitiseProteinBlastDB.clean_db
    }
}


task SanitiseProteinBlastDB {
    input {
        Array[File]? dbs
    }

    output {
        File clean_db = "output.db"
    }

    command <<<
        sanitize_sequence_db.py -cstop ~{sep=" " dbs} | gt seqtransform -addstopaminos -width "60" > "output.db"

    >>>
}

task TransdecoderLongOrf {
    input {
        File reference
        Int minprot = 20 # Minimum protein length (TODO make parameter)
        String gencode = "Universal" # Genetic code (TODO Make parameter) 
        # More gencode info https://github.com/TransDecoder/TransDecoder/blob/master/PerlLib/Nuc_translator.pm
    }

    output {
        File orfs = "longest_orfs.gff3"
        File pep = "longest_orfs.pep"
    }

    command <<<
    Transdecoder.LongOrfs -m "~{minprot}" -t "~{reference}" -G "~{gencode}"
    >>>
}