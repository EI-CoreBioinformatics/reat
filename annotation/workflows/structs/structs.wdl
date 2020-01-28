version 1.0

struct LRSample {
    String name
    String strand
    File LR
}

struct PRSample {
    String name
    String strand
    File R1
    File R2
}

struct SESample {
    String name
    String strand
    File SR
}

struct AlignedSample {
    String name
    String strand
    String aligner
    File bam
}

struct IndexedAlignedSample {
    String name
    String strand
    String aligner
    File bam
    File index
}

struct AlignedSampleStats {
    String name
    File stats
    Array[File] plots
}

struct AssembledSample {
    String name
    String strand
    File assembly
}

struct IndexedReference {
    File fasta
    File fai
}

struct LabeledFasta {
    String label
    File fasta
}