version 1.0

struct LRSample {
    String name
    String strand
    Array[File]+ LR
    Int? score
    Boolean? is_ref
    Boolean? exclude_redundant
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
    Boolean? exclude_redundant
}

struct SESample {
    String name
    String strand
    Array[File]+ SR
    Boolean merge
    Int? score
    Boolean? is_ref
    Boolean? exclude_redundant
}

struct AlignedSample {
    String name
    String strand
    String aligner
    Boolean merge
    Array[File] bam
    Int score
    Boolean is_ref
    Boolean exclude_redundant
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
    Int score
    Boolean is_ref
    Boolean exclude_redundant
}

struct IndexedReference {
    File fasta
    File fai
}

struct LabeledFasta {
    String label
    File fasta
}