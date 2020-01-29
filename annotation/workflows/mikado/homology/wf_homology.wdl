version 1.0

import "../../common/structs.wdl"
import "../../common/tasks.wdl" as tasks
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

workflow wf_homology {
    input {
        File reference
        String program
        File? protein_db
    }

    # Split sequence file
    if (defined(protein_db)) {
        File def_db = select_first([protein_db])
        call tasks.SplitSequences {
            input:
            sequences_file = reference
        }

        call prt_aln.SanitiseProteinBlastDB {
            input:
            db = def_db
        }

        if (program == "blastx") {
            call prt_aln.BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }
            scatter (seq_file in SplitSequences.seq_files) {
                call prt_aln.BlastAlign {
                    input:
                    index = BlastIndex.index,
                    query = seq_file
                }
            }
        }

        if (program == "diamond") {
            call prt_aln.DiamondIndex {
                input:
                target = def_db
            }
            scatter (seq_file in SplitSequences.seq_files) {
                call prt_aln.DiamondAlign {
                    input:
                    index = DiamondIndex.index,
                    query = seq_file
                }
            }
        }
        Array[File] maybe_align = select_first([BlastAlign.out, DiamondAlign.out])
    }

    output {
        Array[File]? homology = maybe_align
        File? homology_clean_db = SanitiseProteinBlastDB.clean_db
    }
}