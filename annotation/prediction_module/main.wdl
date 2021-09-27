version 1.0

workflow ei_prediction {
	input {
		File reference_genome
		Array[File]? transcriptome_models
		Array[File]? homology_models

		File protein_validation_database
	}

	# Preprocess gene models
	if (defined(transcriptome_models)) {
		scatter (models in select_first([transcriptome_models])) {
			call PreprocessFiles as PreprocessTranscriptomic{
				input:
				models = models
			}
		}

		Array[File] processed_transcriptome = PreprocessTranscriptomic.out
	}

	if (defined(homology_models)) {
		scatter (models in select_first([homology_models])) {
			call PreprocessFiles as PreprocessHomology {
				input:
				models = models
			}
		}

		Array[File] processed_homology = PreprocessHomology.out
	}

	Array[File] all_models = flatten(select_all([processed_transcriptome, processed_homology]))

	# Generate protein files for the input models
	call GenerateModelProteins {
		input:
		genome = reference_genome,
		models = all_models
	}

	# TODO: Maybe split proteins for alignments?

	call IndexProteinsDatabase {
		input:
		db = protein_validation_database
	}

	call AlignProteins {
		input:
		proteins = GenerateModelProteins.proteins,
		db = IndexProteinsDatabase.diamond_index
	}

	# Check models for 'full-length'
	call LengthChecker {
		input:
		genome = reference_genome,
		models = all_models,
		protein_models = GenerateModelProteins.proteins,
		hits = AlignProteins.hits
	}

	call SelfBlastFilter {
		input:
		genome = reference_genome,
		clustered_models = LengthChecker.clustered_models,
		classification = LengthChecker.classification
	}

	# Generate CodingQuarry predictions

	call CodingQuarry {
		input:
		genome = reference_genome,
		transcripts = LengthChecker.clustered_models
	}

	# Generate SNAP predictions

	call SNAP {
		input:
		genome = reference_genome,
		transcripts = LengthChecker.clustered_models
	}

	# Generate GlimmerHMM predictions

	call GlimmerHMM {
		input:
		genome = reference_genome,
		transcripts = LengthChecker.clustered_models
	}

	# Feed all to Augustus in the various configurations

	output {

	}
}

task CodingQuarry {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "codingquarry.predictions.gff"
	}

	Int num_cpus = 8

	command <<<
		CodingQuarry -p ~{num_cpus} -f ~{genome} -t ~{transcripts}
		mv codingquarry.gff codingquarry.predictions.gff
	>>>
}

task PreprocessFiles {
	input {
		File models
	}

	output {
		File out = basename(models) + ".clustered.gtf"
	}

	command <<<
		gffread -F --cluster-only --keep-genes ~{models} > ~{basename(models)}.clustered.gtf
	>>>
}

task SNAP {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "snap.predictions.gff"
	}

	command <<<
		gff_to_zff -a ~{transcripts} > ~{transcripts}.ann
		fathom ~{transcripts}.ann ~{genome} -categorize 1000
		fathom uni.ann uni.dna -export 1000 -plus
		mkdir params
		cd params
		forge ../export.ann ../export.dna
		cd ..
		hmm-assembler.pl ~{genome} params > ~{genome}.hmm
		snap ~{genome}.hmm ~{genome}
		zff_to_gff -z snap.prediction.zff > snap.predictions.gff
	>>>
}

task GlimmerHMM {
	input {
		File genome
		File transcripts
	}

	output {
		File predictions = "glimmer.predictions.gff"
	}

	command <<<
		gff_to_glimmer -a ~{transcripts} > ~{transcripts}.cds
		trainGlimmerHMM ~{genome} ~{transcripts}.cds -d glimmer_training
		glimmhmm.pl $(which glimmerhmm) ~{genome} glimmer_training -g | gffread -g ~{genome} -vE --keep-genes -P > glimmer.predictions.gff
	>>>
}

task GenerateModelProteins {
	input {
		File genome
		Array[File] models
	}

	output {
		File proteins = "proteins.faa"
	}

	command <<<
		cat ~{sep=" " models} | gffread --stream -g ~{genome} -y proteins.faa
	>>>
}

task IndexProteinsDatabase {
	input {
		File db
	}

	output {
		File diamond_index = basename(db) + ".dmnd"
	}

	command <<<
		diamond makedb -d ~{basename(db)}.dmnd --in ~{db}
	>>>
}

task AlignProteins {
	input {
		File db
		File proteins
	}

	output {
		File hits = "diamond.hits.tsv"
	}

	command <<<
		diamond blastp -d ~{db} -q ~{proteins} -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > diamond.hits.tsv
	>>>
}

task LengthChecker {
	input{
		File genome
		Array[File] models
		File protein_models
		File hits
	}

	output {
		File classification = "transcripts.classification.tsv"
		File clustered_models = "all_models.clustered.gff"
		File gold = "gold.gff"
		File silver = "silver.gff"
		File bronze = "bronze.gff"
	}

	command <<<
		gffread -g ~{genome} -F --cluster-only --keep-genes -P ~{sep=" " models} > all_models.clustered.gff
		classify_transcripts -b ~{hits} -t all_models.clustered.gff
	>>>
}

task SelfBlastFilter{
	input {
		File genome
		File clustered_models
		File classification
	}

	output {
		File non_redundant_models = "nr_models_80cov_80id.gff"
	}

	command <<<
		gffread -y proteins.faa -g ~{genome} ~{clustered_models}
		diamond makedb --db self -p 8 --in proteins.faa
		diamond blastp -p 8 -d self -q proteins.faa -f6 > self.hits.tsv
		filter_self_hits -b self.hits.tsv -t ~{clustered_models} > nr_models_80cov_80id.gff
	>>>
}