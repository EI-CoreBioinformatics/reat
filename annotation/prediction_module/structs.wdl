version development

struct RuntimeAttr {
    String? constraints
    Float? mem_gb
    Int cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? queue
}

struct IndexedReference {
	File fasta
	File? index
}

struct IndexedBAM {
	File bam
	File? index
}