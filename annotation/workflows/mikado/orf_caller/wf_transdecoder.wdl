version 1.0

import "../../structs/structs.wdl"
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

workflow wf_transdecoder {
    input {
        File prepared_fa
        String program = "blastx"
        File? db
    }

    call TransdecoderLongOrf {
        input:
            reference = prepared_fa
    }

    if (defined(db)) {
        File def_db = select_first([db])
        call prt_aln.SanitiseProteinBlastDB { # Check if this is an input for both homology + orf_calling or just the same task
            input:
                db = def_db
        }

        if (program == "blastx") {
            call prt_aln.BlastIndex as BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }

            call prt_aln.BlastAlign {
                input:
                query = TransdecoderLongOrf.pep,
                index = BlastIndex.index
            }
        }
        
        if (program == "diamond") {
            call prt_aln.DiamondIndex as DiamondIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }

            call prt_aln.DiamondAlign {
                input:
                query = TransdecoderLongOrf.pep,
                index = DiamondIndex.index
            }
        }

        File maybe_aligned_orfs = select_first([BlastAlign.out, DiamondAlign.out])
    }

    output {
        File long_orfs = TransdecoderLongOrf.orfs
        File? clean_db = SanitiseProteinBlastDB.clean_db
        File orfs = select_first([TransdecoderLongOrf.orfs, maybe_aligned_orfs])
    }
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