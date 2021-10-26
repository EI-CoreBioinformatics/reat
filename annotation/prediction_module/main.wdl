version development

import "structs.wdl"
import "augustus.wdl"

workflow ei_prediction {
	input {
		IndexedReference reference_genome
		File extrinsic_config
		String species
		Directory augustus_config_path
		Array[File]? transcriptome_models  # Classify and divide into Gold, Silver and Bronze
		Array[File]? homology_models  # Take as is
		File? intron_hints  # Separate into gold == 1.0 score and silver (the rest)
		IndexedBAM? expressed_exon_hints  # Transform into gff passing by bigwig
		File? repeats_gff  # These are passed through to augustus TODO: Use them for the other predictors
		Int flank = 200
		Int kfold = 8
		Boolean optimise_augustus = false
		Boolean force_train = false
		Array[File]? augustus_runs # File with SOURCE PRIORITY pairs defining the augustus configurations
		File protein_validation_database
	}

	call SoftMaskGenome {
		input:
		genome = reference_genome
	}

	if (! defined(reference_genome.index)) {
		call IndexGenome {
			input:
			genome = SoftMaskGenome.fasta
		}

		IndexedReference new_reference_genome = object {fasta: SoftMaskGenome.soft_masked_genome, index: IndexGenome.index}
	}

	IndexedReference def_reference_genome = select_first(
											[
											new_reference_genome,
											object {fasta: SoftMaskGenome.soft_masked_genome, index: reference_genome.index}
											])

	# Preprocess gene models
	if (defined(transcriptome_models)) {
		scatter (models in select_first([transcriptome_models])) {
			call PreprocessFiles as PreprocessTranscriptomic {
				input:
				models = models
			}
		}

		Array[File] processed_transcriptome = PreprocessTranscriptomic.out
		# Generate protein files for the input models
		call GenerateModelProteins {
			input:
			genome = def_reference_genome,
			models = processed_transcriptome
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
			genome = def_reference_genome,
			models = processed_transcriptome,
			protein_models = GenerateModelProteins.proteins,
			hits = AlignProteins.hits
		}

		call SelfBlastFilter {
			input:
			genome = def_reference_genome,
			clustered_models = LengthChecker.clustered_models,
			classification = LengthChecker.nr_classification
		}

		# If we have enough models we train with UTR, otherwise we use the 'extended' training set and train without UTR
		if (SelfBlastFilter.num_utr_models >= 1200) {
			File maybe_utr = SelfBlastFilter.non_redundant_models_with_utr
			Int maybe_num_utr = SelfBlastFilter.num_utr_models
			Boolean maybe_train_utr = true
		}

		if (SelfBlastFilter.num_utr_models < 1200) {
			File maybe_no_utr = SelfBlastFilter.non_redundant_models_without_utr
			Int maybe_num_no_utr = SelfBlastFilter.num_noutr_models
			Boolean maybe_train_no_utr = false
		}

		File def_training_models = select_first([maybe_utr, maybe_no_utr])
		Int num_models = select_first([maybe_num_utr, maybe_num_no_utr])
		Boolean train_utr_ = select_first([maybe_train_utr, maybe_train_no_utr])
	}


	Boolean train_utr = select_first([train_utr_, false])

	if (defined(homology_models)) {
		scatter (models in select_first([homology_models])) {
			call PreprocessFiles as PreprocessHomology {
				input:
				models = models
			}
		}
		Array[File] processed_homology = PreprocessHomology.out
	}

	# Generate CodingQuarry predictions
	# Considers lowercase as masked by default
	if (num_models > 1500) {
		call CodingQuarry {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models])
		}
	}

	# Generate SNAP predictions
	# With the -lcmask option treats lowercase as masked nts
	if (num_models > 1300) {
		call SNAP {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models])
		}
	}

	# Generate GlimmerHMM predictions
	# Does not consider any masking
	if (num_models > 2000) {
		call GlimmerHMM {
			input:
			genome = def_reference_genome,
			transcripts = select_first([def_training_models])
		}
	}

	call SelectAugustusTestAndTrain {
		input:
		models = select_first([def_training_models]),
		reference = def_reference_genome,
		flank = flank
	}

	# Augustus
	# Using --softmasking=1 considers lowercase nts as masked
	if (num_models > 1000) {
		IndexedReference augustus_genome = object {fasta: SoftMaskGenome.unmasked_genome, index: def_reference_genome.index }
	# Feed all to Augustus in the various configurations
	# Generate training input for augustus (training + test sets)
	# Transform tranining set to GeneBank format

		call Species {
			input:
			base_config = augustus_config_path,
			species = species
		}

		if (!Species.existing_species || force_train) {
			call etraining as base_training {
				input:
				models = SelectAugustusTestAndTrain.train,
				species = species,
				config_path = Species.config_path,
				with_utr = train_utr
			}

			call augustus.Augustus as AugustusTest {
				input:
				reference = augustus_genome,
				extrinsic_config = extrinsic_config,
				with_utr = train_utr,
				species = species,
				config_path = base_training.improved_config_path,
				id = "test"
			}

			if (optimise_augustus) {
				call OptimiseAugustus {
					input:
					num_fold = kfold,
					with_utr = train_utr,
					species = species,
					config_path = base_training.improved_config_path,
					models = SelectAugustusTestAndTrain.train
				}

				call etraining as train_after_optimise {
					input:
					models = SelectAugustusTestAndTrain.train,
					species = species,
					config_path = OptimiseAugustus.optimised_config_path,
					with_utr = train_utr
				}
			}
		}

		Directory final_augustus_config = select_first([train_after_optimise.improved_config_path, base_training.improved_config_path, Species.config_path])

		if (defined(augustus_runs)) {
			Array[Int] counter = range(length(select_first([augustus_runs])))
			scatter (augustus_run_and_index in zip(select_first([augustus_runs]), counter)) {
				call augustus.wf_augustus {
					input:
					reference_genome = augustus_genome,
					species = species,
					augustus_config = final_augustus_config,
					extrinsic_config = extrinsic_config,
					intron_hints = intron_hints,
					expressed_exon_hints = expressed_exon_hints,
					homology_models = homology_models,
					repeats_gff = repeats_gff,
					train_utr = train_utr,
					gold_models = LengthChecker.gold,
					silver_models = LengthChecker.silver,
					bronze_models = LengthChecker.bronze,
					all_models = LengthChecker.clustered_models,
					hints_source_and_priority = augustus_run_and_index.left,
					run_id = augustus_run_and_index.right
				}
			}
		}
	}


	output {
		File? gold_models = LengthChecker.gold
		File? silver_models = LengthChecker.silver
		File? bronze_models = LengthChecker.bronze
		File? training_models = def_training_models
		File? coding_quarry_predictions = CodingQuarry.predictions
		File? snap_predictions = SNAP.predictions
		File? glimmerhmm_predictions = GlimmerHMM.predictions
		Array[File]? augustus_predictions = wf_augustus.predictions
		Directory? augustus_config = final_augustus_config
	}
}

task SoftMaskGenome {
	input {
		IndexedReference genome
		File? repeats_gff
	}

	output {
		File soft_masked_genome = sub(basename(genome.fasta), "\.(fasta|fa)", ".softmasked.fa")
		File hard_masked_genome = sub(basename(genome.fasta), "\.(fasta|fa)", ".hardmasked.fa")
		File unmasked_genome = sub(basename(genome.fasta), "\.(fasta|fa)", ".unmasked.fa")
	}

	command <<<
		cat ~{genome.fasta} | python3 -c "
		import sys;
		for line in sys.stdin:
			if line.startswith('>'):
				print(line, end='')
        		continue
			print(line.upper(), end='')
		" > ~{sub(basename(genome.fasta), "\.(fasta|fa)", ".unmasked.fa")}

		rep_file=~{repeats_gff}
		if [[ ~{repeats_gff} == "" ]];
		then
			touch repeats.gff
			rep_file="repeats.gff"
		fi

		bedtools maskfasta -soft -fi ~{genome.fasta} -bed <(gffread --bed $rep_file) -fo ~{sub(basename(genome.fasta), "\.(fasta|fa)", ".softmasked.fa")}
		bedtools maskfasta -mc 'N' -fi ~{genome.fasta} -bed <(gffread --bed $rep_file) -fo ~{sub(basename(genome.fasta), "\.(fasta|fa)", ".hardmasked.fa")}
	>>>
}

task IndexGenome {
	input {
		IndexedReference genome
	}

	output {
		File index = basename(genome.fasta) + ".fai"
	}

	command <<<
		ln -s ~{genome.fasta}
		samtools faidx ~{basename(genome.fasta)}
	>>>
}

task OptimiseAugustus {
	input {
		Int num_fold
		String species
		Directory config_path
		File models
		Boolean with_utr
	}

	output {
		Directory optimised_config_path = "config"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{config_path} config
		optimize_augustus.pl --AUGUSTUS_CONFIG_PATH=config --species=~{species} ~{models} --cpus=~{num_fold} --kfold=~{num_fold} ~{if with_utr then "--UTR=on" else ""}
	>>>
}

task Species {
	input {
		Directory base_config
		String species
	}

	output {
		Directory config_path = "config"
		Boolean existing_species = read_boolean("existed")
	}

	command <<<
		set -euxo pipefail
		if [ ! -d "~{base_config}/species/~{species}" ];
		then
			mkdir -p config/species
			cp -r ~{base_config}/species/generic config/species/generic
			new_species.pl --species=~{species} --AUGUSTUS_CONFIG_PATH=config
			sed '/^UTR/s/.*/UTR\ton/' config/species/~{species}/~{species}_parameters.cfg > config/species/~{species}/~{species}_parameters.cfg.utr
			mv config/species/~{species}/~{species}_parameters.cfg.utr config/species/~{species}/~{species}_parameters.cfg
			echo "false" > existed
		else
			cp -r ~{base_config}/species/~{species} config/species/~{species}
			echo "true" > existed
		fi
		cp -R ~{base_config}/{cgp,extrinsic,model,profile} config/
	>>>
}

task etraining {
	input {
		File models
		String species
		Boolean with_utr
		Directory config_path
	}

	output {
		Directory improved_config_path = "config"
	}

	command <<<
		ln -s ~{config_path} config
		etraining --AUGUSTUS_CONFIG_PATH=config --species=~{species} ~{models}
	>>>
}

task SelectAugustusTestAndTrain {
	input {
		File models
		IndexedReference reference
		Int flank = 200
	}

	output {
		File test = "test.gb"
		File train = "train.gb"
	}

	command <<<
		set -euxo pipefail
		generate_augustus_test_and_train ~{models}
		gff2gbSmallDNA.pl test.gff ~{reference.fasta} ~{flank} test.gb
		gff2gbSmallDNA.pl train.gff ~{reference.fasta} ~{flank} train.gb
	>>>
}

task CodingQuarry {
	input {
		IndexedReference genome
		File transcripts
	}

	output {
		File predictions = "codingquarry.predictions.gff"
	}

	Int num_cpus = 8

	command <<<
		set -euxo pipefail
		CodingQuarry -p ~{num_cpus} -f ~{genome.fasta} -t ~{transcripts}
		mv out/PredictedPass.gff3 codingquarry.predictions.gff
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
		set -euxo pipefail
		gffread -F --cluster-only --keep-genes ~{models} > ~{basename(models)}.clustered.gtf
	>>>
}

task SNAP {
	input {
		IndexedReference genome
		File transcripts
	}

	output {
		File predictions = "snap.predictions.gff"
	}

	command <<<
		set -euxo pipefail
		gff_to_zff -a ~{transcripts} > ~{transcripts}.ann
		fathom ~{transcripts}.ann ~{genome.fasta} -categorize 1000
		fathom uni.ann uni.dna -export 1000 -plus
		mkdir params
		cd params
		forge ../export.ann ../export.dna
		cd ..
		hmm-assembler.pl ~{genome.fasta} params > ~{genome.fasta}.hmm
		snap ~{genome.fasta}.hmm ~{genome.fasta} > snap.prediction.zff
		zff_to_gff -z snap.prediction.zff > snap.predictions.gff
	>>>
}

task GlimmerHMM {
	input {
		IndexedReference genome
		File transcripts
	}

	output {
		File predictions = "glimmer.predictions.gff"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		gff_to_glimmer -a ~{transcripts} > ~{transcripts}.cds
		trainGlimmerHMM ~{basename(genome.fasta)} ~{transcripts}.cds -d glimmer_training
		glimmhmm.pl $(which glimmerhmm) ~{basename(genome.fasta)} glimmer_training -g | gffread -g ~{basename(genome.fasta)} -vE --keep-genes -P > glimmer.predictions.gff
	>>>
}

task GenerateModelProteins {
	input {
		IndexedReference genome
		Array[File] models
	}

	output {
		File proteins = "proteins.faa"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		cat ~{sep=" " models} | gffread --stream -g ~{basename(genome.fasta)} -y proteins.faa
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
		set -euxo pipefail
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
		set -euxo pipefail
		diamond blastp -d ~{db} -q ~{proteins} -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > diamond.hits.tsv
	>>>
}

task LengthChecker {
	input{
		IndexedReference genome
		Array[File] models
		File protein_models
		File hits
	}

	output {
		File classification = "transcripts.classification.tsv"
		File nr_classification = "transcript.locus.classification.tsv"
		File clustered_models = "all_models.clustered.gff"
		File gold = "gold.gff"
		File silver = "silver.gff"
		File bronze = "bronze.gff"
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		gffread -g ~{basename(genome.fasta)} -F --cluster-only --keep-genes -P ~{sep=" " models} > all_models.clustered.gff
		classify_transcripts -b ~{hits} -t all_models.clustered.gff
	>>>
}

task SelfBlastFilter{
	input {
		IndexedReference genome
		File clustered_models
		File classification
	}

	output {
		File non_redundant_models_with_utr = "with_utr.gff"
		File non_redundant_models_without_utr = "without_utr.gff"
		Int num_utr_models = read_int("num_models_utr.int")
		Int num_noutr_models = read_int("num_models_noutr.int")
	}

	command <<<
		set -euxo pipefail
		ln -s ~{genome.fasta}
		ln -s ~{genome.index}
		gffread -y proteins.faa -g ~{basename(genome.fasta)} ~{clustered_models}

		diamond makedb --db self -p 8 --in proteins.faa
		diamond blastp -p 8 -d self -q proteins.faa -f6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop > self.hits.tsv
		filter_self_hits -b self.hits.tsv -t ~{clustered_models} -c ~{classification}
		awk '$3 == "gene"' with_utr.gff|wc -l > num_models_utr.int
		awk '$3 == "gene"' without_utr.gff|wc -l > num_models_noutr.int
	>>>
}
