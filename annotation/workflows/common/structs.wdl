version 1.0

struct LRSample {
    String name
    String strand
    Array[File]+ LR
    Int? score
    Boolean? is_ref
}

struct ReadPair {
    File R1
    File R2
}

struct PRSample {
    String name
    String strand
    Array[ReadPair]+ read_pair
    Int? score
    Boolean? is_ref
}

struct SESample {
    String name
    String strand
    Array[File]+ SR
    Int? score
    Boolean? is_ref
}

struct AlignedSample {
    String name
    String strand
    String aligner
    Array[File] bam
}

struct IndexedBam{
    File bam
    File index
}

struct IndexedAlignedSample {
    String name
    String strand
    String aligner
    Array[IndexedBam] index_bam
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