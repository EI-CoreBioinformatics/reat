version 1.0

import "../../common/structs.wdl"
import "../../common/tasks.wdl"
import "../align_protein/wf_protein_aligner.wdl" as prt_aln

struct TransdecoderChunk {
    File pep
    File cds
    File gff
}

workflow wf_transdecoder {
    input {
        File prepared_transcripts
        String program = "blast"
        File? orf_proteins
        Boolean refine_start_codons = true
    }

# Prepare the protein alignment database
    if (refine_start_codons && defined(orf_proteins)) {
        File def_db = select_first([orf_proteins])
        call prt_aln.SanitiseProteinBlastDB {
            input:
                db = def_db
        }
        if (program == "blast") {
            call prt_aln.BlastIndex as BlastIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }
        }
        if (program == "diamond") {
            call prt_aln.DiamondIndex as DiamondIndex {
                input:
                target = SanitiseProteinBlastDB.clean_db
            }
        }
    }


# Step 1 Subdivide the prepared_transcripts
    call tasks.SplitSequences {
        input:
        sequences_file = prepared_transcripts,
        prefix = "transdecoder_split"
    }

# Step 2 Generate the longest_orfs and nucleotide frequency table (just call LongOrf in a scatter)
    scatter (chunk in SplitSequences.seq_files) {
        call TransdecoderLongOrf {
            input:
            prepared_transcripts = chunk
        }
    }

# Step 3 Select top longest CDS entries from the scatter and merge all NT freqs table
# Step 3 Filter similar proteins (can use TDC exclude_similar)
# Step 3 Select top longest for the final set (can use TDC get_top_longest)
# Step 3 Train the markov model
    call Select_TrainMM_Score {
        input:
        cds_files = TransdecoderLongOrf.cds,
        base_freqs = TransdecoderLongOrf.base_freqs
    }

# Step 7 preparation
    if(refine_start_codons && defined(orf_proteins)) {
        call Train_PWM {
            input:
            prepared_transcripts = prepared_transcripts,
            top_cds = Select_TrainMM_Score.top_cds
        }
    }

# Step 6.1 Scatter again to score all cds entries using "hexamer_scores_file"

# Step 6.2 (optional) Run blastp on the chunks to filter
    scatter (chunk in TransdecoderLongOrf.chunks) {
        call Calculate_Scores {
            input:
            cds = chunk.cds,
            scores = Select_TrainMM_Score.hexamer_scores_file
        }

        if (refine_start_codons && defined(orf_proteins)) {
            if (program == "blast") {
                call prt_aln.BlastAlign {
                    input:
                    index = select_first([BlastIndex.index]),
                    query = chunk.pep,
                    extra = "-evalue 1e-5 -max_target_seqs 1",
                    outfmt = "6",
                    blast_type = "blastp",
                    output_filename = "mikado_blast_orfcalling.txt",
                }
            }
            if (program == "diamond") {
                call prt_aln.DiamondAlign {
                    input:
                    index = select_first([DiamondIndex.index]),
                    query = chunk.pep,
                    extra = "--evalue 1e-5 --max-target-seqs 1",
                    output_filename = "mikado_diamond_orfcalling.txt",
                    blast_type = "blastp"
                }
            }
            File alignments = select_first([BlastAlign.out, DiamondAlign.out])
        }

# Step 6.3 In the same scatter generate "best_candidates.gff3" per chunk
        call SelectBestOrfs {
            input:
            gff3 = chunk.gff,
            cds_scores = Calculate_Scores.scored_cds,
            RETAIN_LONG_ORFS_MIN_LENGTH = 10,
            blastp_hits_file = alignments,
            single_best_only = false
        }

# Step 7 (optional) Refine start codons
        if(refine_start_codons && defined(orf_proteins)) {
            call RefineStartCodons {
                input:
                transcripts = prepared_transcripts,
                refinement_model = select_first([Train_PWM.refinement_model]),
                gff = SelectBestOrfs.best_candidates,
            }
        }
        File chunk_final_models = select_first([RefineStartCodons.refined_models, SelectBestOrfs.best_candidates])
    }

# Step 8 concatenate all outputs into a final output file

    call tasks.MergeFiles {
        input:
            files_to_merge = chunk_final_models,
            output_filename = "best_candidates.gff3"
    }


    output {
        File final_models = MergeFiles.out
        File? clean_db = SanitiseProteinBlastDB.clean_db
        # File final_orfs = TransdecoderPredict.bed
    }
}

task TransdecoderLongOrf {
    input {
        File prepared_transcripts
        Int minprot = 100 # Minimum protein length
        String gencode = "universal" # Genetic code
# Options for gencode: 
    # Universal
    # Tetrahymena
    # Acetabularia
    # Ciliate
    # Dasycladacean
    # Hexamita
    # Candida
    # Euplotid
    # SR1_Gracilibacteria
    # Pachysolen_tannophilus
    # Mesodinium
    # Peritrich
    # Mitochondrial-Vertebrates
    # Mitochondrial-Yeast
    # Mitochondrial-Invertebrates
    # Mitochondrial-Protozoan
    # Mitochondrial-Echinoderm
    # Mitochondrial-Ascidian
    # Mitochondrial-Flatworm
    # Mitochondrial-Chlorophycean
    # Mitochondrial-Trematode
    # Mitochondrial-Scenedesmus_obliquus
    # Mitochondrial-Thraustochytrium
    # Mitochondrial-Pterobranchia
# For more gencode info https://github.com/TransDecoder/TransDecoder/blob/master/PerlLib/Nuc_translator.pm
        RuntimeAttr? runtime_attr_override
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        TransdecoderChunk chunks = object { cds: "transdecoder/longest_orfs.cds", pep: "transdecoder/longest_orfs.pep", gff: "transdecoder/longest_orfs.gff3" }
        File cds = "transdecoder/longest_orfs.cds"
        File gff3 = "transdecoder/longest_orfs.gff3"
        File pep = "transdecoder/longest_orfs.pep"
        File base_freqs = "transdecoder/base_freqs.dat"
    }

    command <<<
        set -euxo pipefail
        Transdecoder.LongOrfs -m "~{minprot}" -t "~{prepared_transcripts}" -G "~{gencode}" -O transdecoder
    >>>
}

task TransdecoderPredict {
    input {
        File prepared_transcripts
        Array[File] transdecoder_wd
        File? blast_aligned_proteins
        RuntimeAttr? runtime_attr_override
        String workdir = sub(transdecoder_wd[0], "\/(?!.*\/).*", "")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File bed = "mikado_prepared.fasta.transdecoder.bed"
    }

    command <<<
    set -euxo pipefail
    Transdecoder.Predict -t ~{prepared_transcripts} -O ~{workdir} ~{"--retain_blastp_hits " + blast_aligned_proteins}
    >>>

}

task Select_TrainMM_Score {
    input {
        Array[File] cds_files
        Array[File] base_freqs
        Int topORFs_train = 500
        Int red_num = topORFs_train * 10
        Int max_protein_size = 5000
    }

    command <<<
    set -euxo pipefail
    # TODO Write some code to get the longest from an array of files and merge base_freqs
    get_top_longest_from_multiple_files --cds=~{sep="," cds_files} -n ~{red_num} --max_size ~{max_protein_size} > "cds.longest_~{red_num}"
    collect_base_freqs --base_freqs=~{sep="," base_freqs} > "merged_base_freqs"
    exclude_similar_proteins.pl "cds.longest_~{red_num}" > "cds.longest_~{red_num}.nr"
    get_top_longest_fasta_entries.pl "cds.longest_~{red_num}.nr" ~{topORFs_train} > "top_~{topORFs_train}_longest"

    seq_n_baseprobs_to_loglikelihood_vals.pl "top_~{topORFs_train}_longest" "merged_base_freqs" > "hexamer.scores"
    >>>

    output {
        File top_cds = "top_"+topORFs_train+"_longest"
        File hexamer_scores_file = "hexamer.scores"
    }
}

task Train_PWM {
    input {
        File prepared_transcripts
        File top_cds
    }

    output {
        Array[File] refinement_model = glob("start_refinement*")
    }

    command <<<
    set -euxo pipefail
    train_start_PWM.pl --transcripts ~{prepared_transcripts} --selected_orfs ~{top_cds} --out_prefix start_refinement
    >>>
}

task Calculate_Scores {
    input {
        File cds
        File scores
    }

    output {
        File scored_cds = "longest_orfs.cds.scores"
    }

    command <<<
    set -euxo pipefail
    score_CDS_likelihood_all_6_frames.pl "~{cds}" "~{scores}" > "longest_orfs.cds.scores"
    >>>
}

task SelectBestOrfs {
    input {
        File gff3
        File cds_scores
        Int RETAIN_LONG_ORFS_MIN_LENGTH
        File? blastp_hits_file
        Boolean single_best_only
    }

    output {
        File best_candidates = "best_candidates.gff3"
    }

    command <<<
    select_best_ORFs_per_transcript.pl \
    --gff3_file ~{gff3} \
    --cds_scores ~{cds_scores} \
    --min_length_auto_accept ~{RETAIN_LONG_ORFS_MIN_LENGTH} \
    ~{"--blast_hits " + blastp_hits_file} \
    ~{if (single_best_only) then "--single_best_orf" else ""} > "best_candidates.gff3"
    >>>
}

task RefineStartCodons {
    input {
        File gff
        File transcripts
        Array[File] refinement_model
    }

    output {
        File refined_models = "best_candidates.gff3"
    }

    command <<<
    # TODO Recreate the directory required by TDC from refinement_model
    mkdir tmp_workdir
    ln -s ~{sep=" " refinement_model} tmp_workdir/
    start_codon_refinement.pl --transcripts ~{transcripts} --workdir tmp_workdir --gff3_file ~{gff} > best_candidates.gff3
    >>>
}
