version development

import "structs.wdl"

workflow wf_augustus {
	input {
		Array[File]? single_seqs
		Array[File]? many_seqs
		String species
		Boolean train_utr
		File extrinsic_config
		Directory augustus_config
		Int chunk_size = 3000000
		Int overlap_size = 100000
		String run_id
		File? intron_hints
		File? expressed_exon_hints
		Array[File]? homology_models
		File? repeats_gff
		File? gold_models
		File? silver_models
		File? bronze_models
		File? all_models
		File? hq_assembly_models
		File? lq_assembly_models
		File? hq_protein_alignment_models
		File? lq_protein_alignment_models
		File? hints_source_and_priority
	}

	if (defined(hints_source_and_priority)) {
		call ReadSourceAndPriority {
			input:
			hints_source_and_priority = select_first([hints_source_and_priority])
		}
	}

	String alignment_hints_source = select_first([ReadSourceAndPriority.alignment_hints_source, '#'])
	Int alignment_hints_priority = select_first([ReadSourceAndPriority.alignment_hints_priority, 0])
	String gold_intron_hints_source = select_first([ReadSourceAndPriority.gold_intron_hints_source, '#'])
	Int gold_intron_hints_priority = select_first([ReadSourceAndPriority.gold_intron_hints_priority, '0'])
	String silver_intron_hints_source = select_first([ReadSourceAndPriority.silver_intron_hints_source, '#'])
	Int silver_intron_hints_priority = select_first([ReadSourceAndPriority.silver_intron_hints_priority, '0'])
	String protein_hints_source = select_first([ReadSourceAndPriority.protein_hints_source, '#'])
	Int protein_hints_priority = select_first([ReadSourceAndPriority.protein_hints_priority, 0])
	String gold_model_hints_source = select_first([ReadSourceAndPriority.gold_model_hints_source, '#'])
	Int gold_model_hints_priority = select_first([ReadSourceAndPriority.gold_model_hints_priority, 0])
	String silver_model_hints_source = select_first([ReadSourceAndPriority.silver_model_hints_source, '#'])
	Int silver_model_hints_priority = select_first([ReadSourceAndPriority.silver_model_hints_priority, 0])
	String bronze_model_hints_source = select_first([ReadSourceAndPriority.bronze_model_hints_source, '#'])
	Int bronze_model_hints_priority = select_first([ReadSourceAndPriority.bronze_model_hints_priority, 0])
	String all_model_hints_source = select_first([ReadSourceAndPriority.all_model_hints_source, '#'])
	Int all_model_hints_priority = select_first([ReadSourceAndPriority.all_model_hints_priority, 0])

	String hq_assembly_hints_source = select_first([ReadSourceAndPriority.hq_assembly_hints_source, '#'])
	Int hq_assembly_hints_priority = select_first([ReadSourceAndPriority.hq_assembly_hints_priority, 0])
	String lq_assembly_hints_source = select_first([ReadSourceAndPriority.lq_assembly_hints_source, '#'])
	Int lq_assembly_hints_priority = select_first([ReadSourceAndPriority.lq_assembly_hints_priority, 0])
	String hq_protein_alignments_hints_source = select_first([ReadSourceAndPriority.hq_protein_alignment_hints_source, '#'])
	Int hq_protein_alignments_hints_priority = select_first([ReadSourceAndPriority.hq_protein_alignment_hints_priority, 0])
	String lq_protein_alignments_hints_source = select_first([ReadSourceAndPriority.lq_protein_alignment_hints_source, '#'])
	Int lq_protein_alignments_hints_priority = select_first([ReadSourceAndPriority.lq_protein_alignment_hints_priority, 0])

	if (defined(expressed_exon_hints) && alignment_hints_source != '#') {
		call UpdateExonPartSourceAndPriority {
			input:
			gff = select_first([expressed_exon_hints]),
			source = alignment_hints_source,
			priority = alignment_hints_priority
		}
	}

	if (defined(intron_hints)) {
		call PrepareIntronHints {
			input:
			intron_gff = select_first([intron_hints]),
			gold_priority = gold_intron_hints_priority,
			gold_source = gold_intron_hints_source,
			silver_priority = silver_intron_hints_priority,
			silver_source = silver_intron_hints_source
		}
	}

	if (defined(homology_models) && protein_hints_source != '#') {
		call PrepareProteinHints {
			input:
			aligned_proteins_gffs = select_first([homology_models]),
			priority = protein_hints_priority,
			source = protein_hints_source
		}
	}

	if (defined(gold_models) && gold_model_hints_source != '#') {
		call PrepareTranscriptHints as gold {
			input:
			transcripts = [select_first([gold_models])],
			category = "gold",
			source = gold_model_hints_source,
			priority = gold_model_hints_priority
		}
	}

	if (defined(silver_models) && silver_model_hints_source != '#') {
		call PrepareTranscriptHints as silver {
			input:
			transcripts = [select_first([silver_models])],
			category = "silver",
			source = silver_model_hints_source,
			priority = silver_model_hints_priority
		}
	}

	if (defined(bronze_models) && bronze_model_hints_source != '#') {
		call PrepareTranscriptHints as bronze {
			input:
			transcripts = [select_first([bronze_models])],
			category = "bronze",
			source = bronze_model_hints_source,
			priority = bronze_model_hints_priority
		}
	}

	if (defined(all_models) && all_model_hints_source != '#') {
		call PrepareTranscriptHints as all {
			input:
			transcripts = [select_first([all_models])],
			category = "all",
			source = all_model_hints_source,
			priority = all_model_hints_priority
		}
	}

	if (defined(hq_assembly_models) && hq_assembly_hints_source != '#') {
		call PrepareTranscriptHints as hq_assembly {
			input:
			transcripts = [select_first([hq_assembly_models])],
			category = "hq_assembly",
			source = hq_assembly_hints_source,
			priority = hq_assembly_hints_priority
		}
	}

	if (defined(lq_assembly_models) && lq_assembly_hints_source != '#') {
		call PrepareTranscriptHints as lq_assembly {
			input:
			transcripts = [select_first([lq_assembly_models])],
			category = "lq_assembly",
			source = lq_assembly_hints_source,
			priority = lq_assembly_hints_priority
		}
	}


	if (defined(hq_protein_alignment_models) && hq_protein_alignments_hints_source != '#') {
		call PrepareTranscriptHints as hq_protein_alignments {
			input:
			transcripts = [select_first([hq_protein_alignment_models])],
			category = "hq_protein_alignments",
			source = hq_protein_alignments_hints_source,
			priority = hq_protein_alignments_hints_priority
		}
	}

	if (defined(lq_protein_alignment_models) && lq_protein_alignments_hints_source != '#') {
		call PrepareTranscriptHints as lq_protein_alignments {
			input:
			transcripts = [select_first([lq_protein_alignment_models])],
			category = "lq_protein_alignments",
			source = lq_protein_alignments_hints_source,
			priority = lq_protein_alignments_hints_priority
		}
	}

	Array[File] hints_files = select_first([select_all([gold.result, silver.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   PrepareProteinHints.protein_hints,
							   UpdateExonPartSourceAndPriority.sp_gff,
							   hq_assembly.result, lq_assembly.result,
							   hq_protein_alignments.result, lq_protein_alignments.result,
							   repeats_gff]), []])
	Boolean with_hints = length(hints_files) > 0
	if (with_hints)
	{
		call cat as augustus_hints {
			input:
			files = select_all([gold.result, silver.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   PrepareProteinHints.protein_hints,
							   UpdateExonPartSourceAndPriority.sp_gff,
							   hq_assembly.result, lq_assembly.result,
							   hq_protein_alignments.result, lq_protein_alignments.result,
							   repeats_gff]),
			out_filename = "hints_" + run_id + ".gff"
		}
		if (defined(single_seqs)) {
			# Run chunk by chunk
			scatter (sequence in select_first([single_seqs])) {
				call SubdivideSequence {
					input:
					sequence = sequence,
					chunk_size = chunk_size,
					overlap = overlap_size
				}
				Array[String] chunks = read_lines(SubdivideSequence.chunks)
				scatter (chunk in chunks) {
					call AugustusByChunk {
					input:
					reference = select_first([sequence]),
					chunk = chunk,
					extrinsic_config = extrinsic_config,
					hints = augustus_hints.out,
					with_utr = train_utr,
					species = species,
					config_path = augustus_config,
					id = run_id
					}
				}
			}
		}

		if (defined(many_seqs)) {
			# Run file by file
			scatter (seqs in select_first([many_seqs])) {
				call Augustus {
					input:
					reference = seqs,
					extrinsic_config = extrinsic_config,
					hints = augustus_hints.out,
					with_utr = train_utr,
					species = species,
					config_path = augustus_config,
					id = run_id
				}
			}
		}
	}

	if (!with_hints) {
		if (defined(single_seqs)) {
			# Run chunk by chunk
			scatter (sequence in select_first([single_seqs])) {
				call SubdivideSequence as SubdivideSequence_noHints {
					input:
					sequence = sequence,
					chunk_size = chunk_size,
					overlap = overlap_size
				}
				Array[String] chunks_nohints = read_lines(SubdivideSequence_noHints.chunks)
				scatter (chunk_nohints in chunks_nohints) {
					call AugustusByChunk as AugustusByChunk_noHints {
					input:
					reference = select_first([sequence]),
					chunk = chunk_nohints,
					extrinsic_config = extrinsic_config,
					with_utr = train_utr,
					species = species,
					config_path = augustus_config,
					id = run_id
					}
				}
			}
		}

		if (defined(many_seqs)) {
			# Run file by file
			scatter (seqs in select_first([many_seqs])) {
				call Augustus as Augustus_noHints {
					input:
					reference = seqs,
					extrinsic_config = extrinsic_config,
					with_utr = train_utr,
					species = species,
					config_path = augustus_config,
					id = run_id
				}
			}
		}

	}

	if (defined(Augustus.predictions) || defined(Augustus_noHints.predictions)) {
		Array[File] single_aug = flatten(select_all([Augustus.predictions, Augustus_noHints.predictions]))
	}
	if (defined(AugustusByChunk_noHints.predictions) || defined(AugustusByChunk.predictions)) {
		Array[File] multi_aug = flatten(select_all(select_first([AugustusByChunk_noHints.predictions, AugustusByChunk.predictions])))
	}

	if (defined(single_aug) || defined(multi_aug)) {
		call cat as JoinAugustus {
			input:
			files = flatten(select_all([single_aug, multi_aug])),
			out_filename = "augustus.predictions.gff"
		}
	}

	output {
		File? hints = augustus_hints.out
		File? predictions = JoinAugustus.out
	}
}

task SubdivideSequence {
	input {
		File sequence
		Int chunk_size = 5000000
		Int overlap = 100000
	}

	output {
		File chunks = "chunks.txt"
	}

	command <<<
# Parse first line for the seq length
python3 -c "
chunks_file=open('chunks.txt', 'w')
sequence_file=open(~{sequence})
name, seq_length = sequence_file.readline().strip().split(maxsplit=2)
start=1
end=0
# Chunk this sequence and print the 'chunks' to a file
print("This sequence generates the following chunks:")
while end < seq_length:
	end = min(seq_length, start + ~{chunk_size})
	print(f"{start},{end}", file=chunks_file)
	print(f"{name}:{start}-{end}")
	start += end + 1 - ~{overlap}
chunks_file.close()
"
	>>>
}

task AugustusByChunk {
	input {
		File reference
		String chunk
		File extrinsic_config
		Directory config_path
		Boolean with_utr
		String species
		File? hints
		String id
	}

	output {
		File predictions = "augustus_"+id+".predictions.gff"
	}

	command <<<
		pred_start=$(echo "~{chunk}" | awk -F',' '{print $1}')
		pred_end=$(echo "~{chunk}" | awk -F',' '{print $2}')
		ln -s ~{config_path} config
		augustus --AUGUSTUS_CONFIG_PATH=~{config_path} --gff3=on \
		--predictionStart=$pred_start --predictionEnd=$pred_end \
		--UTR=~{if with_utr then "ON" else "OFF"} --stopCodonExcludedFromCDS=true --genemodel=partial \
		--alternatives-from-evidence=true ~{'--hintsfile=' + hints} --noInFrameStop=true \
		--allow_hinted_splicesites=atac --errfile=run~{id}.log --extrinsicCfgFile=~{extrinsic_config} \
		--species=~{species} ~{reference} | grep -v '^#' | awk -v 'OFS=\t' '$2="AUGUSTUS_RUN~{id}"' > augustus_~{id}.predictions.gff
	>>>
}

task cat {
	input {
		Array[File]+ files
		String out_filename
	}

	output {
		File out = out_filename
	}

	command <<<
		cat ~{sep=' ' files} > ~{out_filename}
	>>>
}

task UpdateExonPartSourceAndPriority {
	input {
		File gff
		String source
		Int priority
	}

	output {
		File sp_gff = basename(gff, '.gff') + 'S' + source + 'P' + priority + '.gff'
	}

	command <<<
		cat ~{gff} | sed 's/src=generic_source/src=~{source}/' | sed 's/pri=0/pri=~{priority}/' > ~{basename(gff, '.gff')}S~{source}P~{priority}.gff
	>>>
}

task PrepareIntronHints {
	input {
		File intron_gff
		Int gold_priority = 6
		String gold_source
		Int silver_priority = 4
		String silver_source
	}

	output {
		File? gold_intron_hints = 'gold_junctions_S'+gold_source+'P'+gold_priority+'.gff'
		File? silver_intron_hints = 'silver_junctions_S'+silver_source+'P'+silver_priority+'.gff'
	}

	command <<<
		if [[ '#' != '~{gold_source}' ]]; then
			cat ~{intron_gff} | awk -F "\t" '$6==1 {print $0";pri=~{gold_priority};"}' | \
			sed 's/\tportcullis\t/\tPortcullis_pass_gold_S~{gold_source}P~{gold_priority}\t/' > gold_junctions_S~{gold_source}P~{gold_priority}.gff
		fi

		if [[ '#' != '~{silver_source}' ]]; then
			cat ~{intron_gff} | awk -F "\t" '$6<1 {print $0";pri=~{silver_priority};"}' | \
			sed 's/\tportcullis\t/\tPortcullis_pass_silver_S~{silver_source}P~{silver_priority}\t/' > silver_junctions_S~{silver_source}P~{silver_priority}.gff
		fi
	>>>
}

task ReadSourceAndPriority {
	input {
		File hints_source_and_priority
	}
	# NOTE: Check the generate_augustus_hint_parameters script if the outputs or the workflow requires changes!
	output {
		String 	gold_model_hints_source = read_string("gold.s")
		Int 	gold_model_hints_priority = read_int("gold.p")
		String 	silver_model_hints_source = read_string("silver.s")
		Int 	silver_model_hints_priority = read_int("silver.p")
		String 	bronze_model_hints_source = read_string("bronze.s")
		Int 	bronze_model_hints_priority = read_int("bronze.p")
		String 	all_model_hints_source = read_string("all.s")
		Int 	all_model_hints_priority = read_int("all.p")
		String 	alignment_hints_source = read_string("alignment.s")
		Int 	alignment_hints_priority = read_int("alignment.p")
		String 	repeat_hints_source = read_string("repeat.s")
		Int 	repeat_hints_priority = read_int("repeat.p")
		String	protein_hints_source = read_string("protein.s")
		Int		protein_hints_priority = read_int("protein.p")
		String	gold_intron_hints_source = read_string("gold_intron.s")
		Int		gold_intron_hints_priority = read_int("gold_intron.p")
		String	silver_intron_hints_source = read_string("silver_intron.s")
		Int		silver_intron_hints_priority = read_int("silver_intron.p")
		String	hq_assembly_hints_source = read_string("hq_assembly.s")
		Int		hq_assembly_hints_priority = read_int("hq_assembly.p")
		String	lq_assembly_hints_source = read_string("lq_assembly.s")
		Int		lq_assembly_hints_priority = read_int("lq_assembly.p")
		String	hq_protein_alignment_hints_source = read_string("hq_protein_alignments.s")
		Int		hq_protein_alignment_hints_priority = read_int("hq_protein_alignments.p")
		String	lq_protein_alignment_hints_source = read_string("lq_protein_alignments.s")
		Int		lq_protein_alignment_hints_priority = read_int("lq_protein_alignments.p")
	}

	command <<<
		generate_augustus_hint_parameters ~{hints_source_and_priority}
	>>>
}

task PrepareProteinHints {
	input {
		Array[File] aligned_proteins_gffs
		String source = "P"
		Int priority = 4
	}

	output {
		File protein_hints = 'aligned_proteins.S'+source+'P'+priority+'.gff'
	}

	command <<<
		for i in ~{sep=' ' aligned_proteins_gffs}; do
		name=$(echo ${basename i} | sed 's/.gff//')
		gffread $i | gff_to_aug_hints -P ~{priority} -S ~{source} -s ${name}.protein -t CDS >> aligned_proteins.S~{source}P~{priority}.gff;
		done
	>>>

}

task PrepareTranscriptHints {
	input {
		Array[File] transcripts
		String category
		String source = "E"
		Int priority = 4
	}

	output {
		File result = category+'.transcripts.S'+source+'P'+priority+'.gff'
	}

	command <<<
		for i in ~{sep=' ' transcripts}; do
		name=$(echo $(basename $i) | sed 's/.gff//')
		cat $i | gff_to_aug_hints -P ~{priority} -S ~{source} -s ~{category}.transcripts -t exon >> ~{category}.transcripts.S~{source}P~{priority}.gff;
		done
	>>>
}

task Augustus {
	input {
		File reference
		File extrinsic_config
		Directory config_path
		Boolean with_utr
		String species
		File? hints
		String id
	}

	output {
		File predictions = "augustus_"+id+".predictions.gff"
	}

	command <<<
		ln -s ~{config_path} config
		augustus --AUGUSTUS_CONFIG_PATH=~{config_path} --gff3=on \
		--UTR=~{if with_utr then "ON" else "OFF"} --stopCodonExcludedFromCDS=true --genemodel=partial \
		--alternatives-from-evidence=true ~{'--hintsfile=' + hints} --noInFrameStop=true \
		--allow_hinted_splicesites=atac --errfile=run~{id}.log --extrinsicCfgFile=~{extrinsic_config} \
		--species=~{species} ~{reference} | grep -v '^#' | awk -v 'OFS=\t' '$2="AUGUSTUS_RUN~{id}"' > augustus_~{id}.predictions.gff
	>>>
}
