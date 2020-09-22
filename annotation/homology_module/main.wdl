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
    }

    call PrepareAnnotations {
        input:
        annotation = annotations
    }

    call AlignProteins {
        input:
        genome_proteins = PrepareAnnotations.genome_proteins
    }

    call ScoreAlignments {
        input:
        alignments = AlignProteins.alignments,
        genome_proteins = PrepareAnnotations.genome_proteins
    }

    output {
        File scored_alignments = ScoreAlignments.result
    }
}

task PrepareAnnotations {
    input {
        GenomeAnnotations annotation
        String out_prefix = sub(basename(annotation.genome),  "\.[^/.]+$", "")
        Int min_cds_len = 20 # nts
    }

    output {
        GenomeProteins genome_proteins = object {genome: annotation.genome, annotation_gff: annotation.annotation_gff, proteins: "proteins.fa"}
    }

    # TODO:
    # - Extract filtered proteins from annotations in GFF or GTF format
    command <<<
        extract_coding_sequences -genome ~{annotation.genome} -min_len ~{min_cds_len} -o ~{out_prefix} ~{annotation.annotation_gff}
    >>>
}

task AlignProteins {
    input {
        GenomeProteins genome_proteins
        String species = "Eudicoty" # TODO: Need a lookup table from the original file
        Int min_exon_len = 20
        Int min_coverage = 80
        Int min_identity = 50
    }

    output {
        File alignments = "alignment.gff"
    }

    Int num_cpus = 12

    command <<<
        spaln -WP ~{genome_proteins.genome}
        spaln -t~{num_cpus} -O0,12 -Q7 -T~{species} -d~{sub(basename(genome_proteins.genome), "\.[^/.]+$", "")} -o alignment -yL~{min_exon_len} ~{genome_proteins.protein_sequences}
        sortgrcd -O4 alignment.grd | tee alignment.s | spaln2gff --min_coverage ~{min_coverage} --min_identity ~{min_identity} -s "spaln" alignment.gff
    >>>
}

task ScoreAlignments {
    input {
        File alignments
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
        echo "TODO ~{genome_proteins.genome} ~{genome_proteins.annotation} ~{alignments} ~{max_intron_len}" > scored_alignments.gff
    >>>
}