version 1.0

import "../../common/structs.wdl"
import "../../common/tasks.wdl" as tasks
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

workflow wf_homology {
    input {
        File reference
        String homology_alignment_program
        File? protein_db
        RuntimeAttr? index_resources
        RuntimeAttr? alignment_resources
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

        if (homology_alignment_program == "blast") {
            call prt_aln.BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db,
                runtime_attr_override = index_resources
            }
            scatter (seq_file in SplitSequences.seq_files) {
                call prt_aln.BlastAlign {
                    input:
                    index = BlastIndex.index,
                    blast_type = "blastx",
                    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop",
                    output_filename = "mikado_blast_homology.tsv",
                    query = seq_file,
                    runtime_attr_override = alignment_resources
                }
            }
        }

        if (homology_alignment_program == "diamond") {
            call prt_aln.DiamondIndex {
                input:
                target = def_db,
                runtime_attr_override = index_resources
            }
            scatter (seq_file in SplitSequences.seq_files) {
                call prt_aln.DiamondAlign {
                    input:
                    index = DiamondIndex.index,
                    blast_type = "blastx",
                    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop",
                    output_filename = "mikado_diamond_homology.tsv",
                    query = seq_file,
                    runtime_attr_override = alignment_resources
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