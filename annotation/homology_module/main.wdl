version 1.0

struct GenomeAnnotation {
    File genome
    File annotation_gff
}

struct GenomeProteins {
    File genome
    File annotation_gff
    File protein_sequences
}

workflow ei_homology {
    input {
        Array[GenomeAnnotation] annotations
        File genome_to_annotate
    }

    call IndexGenome {
        genome = genome_to_annotate
    }

    scatter (annotation in annotations) {
        call PrepareAnnotations {
            input:
            annotation = annotation
        }

        GenomeProteins cleaned_up = object { genome: annotation.genome, annotation_gff: PrepareAnnotations.cleaned_up_gff, protein_sequences: PrepareAnnotations.proteins}

        call AlignProteins {
            input:
            genome_to_annotate = genome_to_annotate,
            genome_proteins = cleaned_up
        }

        call ScoreAlignments {
            input:
            alignments = AlignProteins.alignments,
            genome_to_annotate = genome_to_annotate,
            genome_proteins = cleaned_up
        }
    }
    output {
        Array[File] scored_alignments = ScoreAlignments.scored_alignments
    }
}

task IndexGenome {
    input {
        File genome
    }

    output {
        Array[File] genome_index = glob(sub(basename(genome),  "\.[^/.]+$", ""))
    }

    command <<<
        spaln -WP ~{genome_to_annotate}
    >>>
}

task PrepareAnnotations {
    input {
        GenomeAnnotation annotation
        String out_prefix = sub(basename(annotation.annotation_gff),  "\.[^/.]+$", "")
        Int min_cds_len = 20 # nts
    }

    output {
        File cleaned_up_gff = out_prefix + ".gff"
        File proteins = out_prefix + ".proteins.fa"
    }

    # TODO:
    # - Extract filtered proteins from annotations GTF format, the following only works on GFF
    command <<<
        xspecies_cleanup --annotation ~{annotation.annotation_gff} --genome ~{annotation.genome} --min_protein ~{min_cds_len} -y ~{out_prefix}.proteins.fa > ~{out_prefix}.gff
    >>>
}

task AlignProteins {
    input {
        File genome_to_annotate
        GenomeProteins genome_proteins
        String species = "Eudicoty" # TODO: Need a lookup table from the original file
        Int min_exon_len = 20
        Int min_coverage = 80
        Int min_identity = 50
    }

    Int num_cpus = 12
    String out_prefix = sub(basename(genome_proteins.genome), "\.[^/.]+$", "")

    output {
        File alignments = out_prefix+".alignment.gff"
    }

    command <<<
        spaln -t~{num_cpus} -KP -O0,12 -Q7 ~{"-T"+species} -d~{out_prefix} -o ~{out_prefix} -yL~{min_exon_len} ~{genome_proteins.protein_sequences}
        sortgrcd -O4 ~{out_prefix}.grd | tee ~{out_prefix}.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "spaln" > ~{out_prefix}.alignment.gff
    >>>


}

task ScoreAlignments {
    input {
        File alignments
        File genome_to_annotate
        GenomeProteins genome_proteins
        Int max_intron_len = 2000
    }

    output {
        File scored_alignments = "scored_alignments.gff"
    }

    # TODO:
    # - Annotate models with the count of introns exceeding max_intron_len per model
    # - Annotate exons with splice signal (to account for non-canonical splicing)
    # - Compare models to original structure (mgc approach preferably on-line)
    command <<<
        echo "xspecies_scoring --max_intron ~{max_intron_len} --genome ~{genome_to_annotate} --annotation ~{genome_proteins.annotation_gff} ~{alignments}" > scored_alignments.gff
    >>>
}