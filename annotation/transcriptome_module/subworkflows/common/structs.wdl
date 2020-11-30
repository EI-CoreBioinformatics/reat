version 1.0

struct LRSample {
    String name
    String strand
    Array[File]+ LR
    Int? score
    Boolean? is_ref
    Boolean? always_keep
}

struct ReadPair {
    File R1
    File R2
}

struct PRSample {
    String name
    String strand
    Array[ReadPair]+ read_pair
    Boolean merge
    Int? score
    Boolean? is_ref
    Boolean? always_keep
}

struct SESample {
    String name
    String strand
    Array[File]+ SR
    Boolean merge
    Int? score
    Boolean? is_ref
    Boolean? always_keep
}

struct AlignedSample {
    String name
    String strand
    String aligner
    Boolean merge
    Array[File] bam
}

struct IndexedBam{
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