version development

struct IndexedReference {
	File fasta
	File? index
}

struct IndexedBAM {
	File bam
	File? index
}

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
		File? repeats_gff  # These are passed through to augustus as is
		Int flank = 200
		Int kfold = 8
		Boolean optimise_augustus = false
		Boolean force_train = false

		File protein_validation_database
	}

	if (! defined(reference_genome.index)) {
		call IndexGenome {
			input:
			genome = reference_genome
		}

		IndexedReference new_reference_genome = object {fasta: reference_genome.fasta, index: IndexGenome.index}
	}

	IndexedReference def_reference_genome = select_first([new_reference_genome, reference_genome])

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
		genome = def_reference_genome,
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
		genome = def_reference_genome,
		models = all_models,
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
	Boolean train_utr = select_first([maybe_train_utr, maybe_train_no_utr])

	# Generate CodingQuarry predictions
	if (num_models > 1500) {
		call CodingQuarry {
			input:
			genome = def_reference_genome,
			transcripts = def_training_models
		}
	}

	# Generate SNAP predictions
	if (num_models > 1300) {
		call SNAP {
			input:
			genome = def_reference_genome,
			transcripts = def_training_models
		}
	}

	# Generate GlimmerHMM predictions
	if (num_models > 2000) {
		call GlimmerHMM {
			input:
			genome = def_reference_genome,
			transcripts = def_training_models
		}
	}

	call SelectAugustusTestAndTrain {
		input:
		models = def_training_models,
		reference = def_reference_genome,
		flank = flank
	}

	if (num_models > 1000) {

		if (defined(intron_hints)) {
			call PrepareIntronHints {
				input:
				intron_gff = select_first([intron_hints])
			}
		}

		if (defined(expressed_exon_hints)) {
			IndexedBAM def_exon = select_first([expressed_exon_hints])
			if (!defined(def_exon.index)) {
				call IndexBAM {
					input:
					bam = def_exon.bam
				}
				IndexedBAM def_calc_index = object {bam: def_exon.bam, index: IndexBAM.index}
			}

			IndexedBAM def_indexed_bam = select_first([def_calc_index, def_exon])

			call PrepareExonHints {
				input:
				expression_bam = def_indexed_bam,
				dUTP = true
			}
		}

		if (defined(homology_models)) {
			call PrepareProteinHints as protein_hints_P4{
				input:
				aligned_proteins_gffs = select_first([homology_models])
			}

			call PrepareProteinHints as protein_hints_P9{
				input:
				aligned_proteins_gffs = select_first([homology_models]),
				priority = 9
			}
		}

		call PrepareTranscriptHints as gold {
			input:
			transcripts = [LengthChecker.gold],
			category = "gold",
			source = 'M',
			priority = 10
		}

		call PrepareTranscriptHints as gold_e {
			input:
			transcripts = [LengthChecker.gold],
			category = "gold",
			source = 'E',
			priority = 10
		}

		call PrepareTranscriptHints as silver {
			input:
			transcripts = [LengthChecker.silver],
			category = "silver",
			source = 'F',
			priority = 9
		}

		call PrepareTranscriptHints as silver_e {
			input:
			transcripts = [LengthChecker.silver],
			category = "silver",
			source = 'E',
			priority = 9
		}

		call PrepareTranscriptHints as bronze {
			input:
			transcripts = [LengthChecker.bronze],
			category = "bronze",
			source = 'E',
			priority = 8
		}

		call PrepareTranscriptHints as all {
			input:
			transcripts = all_models,
			category = "all",
			source = 'E',
			priority = 7
		}

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

			call Augustus as AugustusTest {
				input:
				reference = def_reference_genome,
				extrinsic_config = extrinsic_config,
				with_utr = train_utr,
				species = species,
				config_path = base_training.improved_config_path
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

		# Prepare run1 evidence:
		# Mikado Gold Source M; Priority 10;
		# Mikado Silver Source M; Priority 9;
		# Mikado Bronze Source E; Priority 8;
		# PacBio Mikado Gold Source M; Priority 10;
		# PacBio Mikado Silver Source M; Priority 9;
		# PacBio Mikado Bronze Source F; Priority 8;
		# PacBio Canonical All Source E; Priority 7;
		# Portcullis Pass Gold Junctions Source E; Priority 6;
		# Portcullis Pass Silver Junctions Source E; Priority 4;
		# Proteins Source P; Priority 4;
		# Normalised Wig Hints Source W; Priority 3;
		# Repeats Source RM; Priority 1;

		call cat as run1_hints {
			input:
			files = select_all([gold.result, silver.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   protein_hints_P4.protein_hints,
							   PrepareExonHints.expression_gff,
							   repeats_gff
							   ]
					),
			out_filename = "run1_hints.gff"
		}
		call Augustus as run1 {
			input:
			reference = def_reference_genome,
			extrinsic_config = extrinsic_config,
			hints = run1_hints.out,
			with_utr = train_utr,
			species = species,
			config_path = final_augustus_config
		}

		# Prepare run2 evidence:
		# Mikado Gold Source M; Priority 10;
		# Mikado Silver Source M; Priority 9;
		# Mikado Bronze Source E; Priority 8;
		# PacBio Mikado Gold Source M; Priority 10;
		# PacBio Mikado Silver Source M; Priority 9;
		# PacBio Mikado Bronze Source F; Priority 8;
		# PacBio Canonical All Source E; Priority 7;
		# Portcullis Pass Gold Junctions Source E; Priority 6;
		# Portcullis Pass Silver Junctions Source E; Priority 4;
		# Proteins Source P; Priority 4;
		# Repeats Source RM; Priority 1;
		call cat as run2_hints {
			input:
			files = select_all([gold.result, silver.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   protein_hints_P4.protein_hints,
							   repeats_gff
							   ]
					),
			out_filename = "run2_hints.gff"
		}

		call Augustus as run2 {
			input:
			reference = def_reference_genome,
			extrinsic_config = extrinsic_config,
			hints = run2_hints.out,
			with_utr = train_utr,
			species = species,
			config_path = final_augustus_config
		}

		# Prepare run3 evidence:
		# Mikado Gold Source E; Priority 10;
		# Mikado Silver Source E; Priority 9;
		# Mikado Bronze Source E; Priority 8;
		# PacBio Mikado Gold Source E; Priority 10;
		# PacBio Mikado Silver Source E; Priority 9;
		# PacBio Mikado Bronze Source E; Priority 8;
		# PacBio Canonical All Source E; Priority 7;
		# Portcullis Pass Gold Junctions Source E; Priority 6;
		# Portcullis Pass Silver Junctions Source E; Priority 4;
		# Proteins Source P; Priority 4;
		# Repeats Source RM; Priority 1;
		call cat as run3_hints {
			input:
			files = select_all([gold_e.result, silver_e.result, bronze.result, all.result,
							   PrepareIntronHints.gold_intron_hints, PrepareIntronHints.silver_intron_hints,
							   protein_hints_P9.protein_hints,
							   repeats_gff
							   ]
					),
			out_filename = "run3_hints.gff"
		}

		call Augustus as run3 {
			input:
			reference = def_reference_genome,
			extrinsic_config = extrinsic_config,
			hints = run3_hints.out,
			with_utr = train_utr,
			species = species,
			config_path = final_augustus_config
		}
	}


	output {
		File gold_models = LengthChecker.gold
		File silver_models = LengthChecker.silver
		File bronze_models = LengthChecker.bronze
		File training_models = def_training_models
		File? coding_quarry_predictions = CodingQuarry.predictions
		File? snap_predictions = SNAP.predictions
		File? glimmerhmm_predictions = GlimmerHMM.predictions
		File? augustus_predictions = run1.predictions
		Directory? augustus_config = final_augustus_config
	}
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

task Augustus {
	input {
		IndexedReference reference
		File extrinsic_config
		Directory config_path
		Boolean with_utr
		String species
		File? hints
	}

	output {
		File predictions = "augustus.predictions.gff"
	}

	command <<<
		ln -s ~{config_path} config
		augustus --AUGUSTUS_CONFIG_PATH=~{config_path} \
		--UTR=~{if with_utr then "ON" else "OFF"} --stopCodonExcludedFromCDS=true --genemodel=partial \
		--alternatives-from-evidence=true ~{'--hintsfile=' + hints} --noInFrameStop=true \
		--allow_hinted_splicesites=atac --errfile=run1.log --extrinsicCfgFile=~{extrinsic_config} \
		--species=~{species} ~{reference.fasta} > augustus.predictions.gff
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

task PrepareIntronHints {
	input {
		File intron_gff
		Int gold_priority = 6
		Int silver_priority = 4
	}

	output {
		File gold_intron_hints = 'gold_junctions_SEP'+gold_priority+'.gff'
		File silver_intron_hints = 'silver_junctions_SEP'+silver_priority+'.gff'
	}

	command <<<
		cat ~{intron_gff} | awk -F "\t" '$6==1 {print $0";pri=~{gold_priority};"}' | \
		sed 's/\tportcullis\t/\tPortcullis_pass_gold_SEP~{gold_priority}\t/' > gold_junctions_SEP~{gold_priority}.gff

		cat ~{intron_gff} | awk -F "\t" '$6<1 {print $0";pri=~{silver_priority};"}' | \
		sed 's/\tportcullis\t/\tPortcullis_pass_silver_SEP~{silver_priority}\t/' > silver_junctions_SEP~{silver_priority}.gff
	>>>
}

task IndexBAM {
	input {
		File bam
	}

	output {
		File index = basename(bam)+'.csi'
	}

	command <<<
		ln -s ~{bam}
		samtools index -@ 4 -c ~{basename(bam)}
	>>>
}

task PrepareExonHints {
	input {
		IndexedBAM expression_bam
		Boolean dUTP
	}

	String bam_name = basename(expression_bam.bam, ".bam")
	String jid = basename(expression_bam.bam)

	output {
		File expression_gff = bam_name + '.exonhints.SWP3.augustus.gff'
	}

	# Uses 2 cpus

	command <<<
		ln -s ~{expression_bam.bam}
		ln -s ~{expression_bam.index}

		if [ ~{dUTP} ];
		then
			strand='1+-,1-+,2++,2--'
		else
			strand='1++,1--,2+-,2-+'
		fi
		touch ~{bam_name}_Forward.wig ~{bam_name}_Reverse.wig ~{bam_name}.wig
		samtools view -H ~{basename(expression_bam.bam)} | grep '^\@SQ'| cut -f2,3 | sed -e 's/SN://g' -e 's/LN://g' > lengths.txt
		bam2wig.py -i ~{basename(expression_bam.bam)} -s lengths.txt -o ~{bam_name} --strand=$strand"

		cat ~{bam_name}.Forward.wig | \
		wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=exonpart --radius=4.5 --pri=3 --strand='+' | \
		sed \"s/\\tw2h\\t/\\tw2h_~{jid}\\t/\" > ~{bam_name}.Forward.exonhints.SWP3.augustus.gff" &

		cat ~{bam_name}.Reverse.wig | \
		awk '/^variableStep/ {gsub(/-/,\"\");print;} !/^variableStep/ {print;}' \ |
		wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=exonpart --radius=4.5 --pri=3 --strand='-' | \
		sed \"s/\\tw2h\\t/\\tw2h_~{jid}\\t/\" > ~{bam_name}.Reverse.exonhints.SWP3.augustus.gff" &

		cat ~{bam_name}.Forward.exonhints.SWP3.augustus.gff ~{bam_name}.Reverse.exonhints.SWP3.augustus.gff > ~{bam_name}.exonhints.SWP3.augustus.gff
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
		name=$(echo ${i} | sed 's/.gff//')
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
		name=$(echo ${i} | sed 's/.gff//')
		cat $i | gff_to_aug_hints -P ~{priority} -S ~{source} -s ${name}.transcripts -t exon >> ~{category}.transcripts.S~{source}P~{priority}.gff;
		done
	>>>

}