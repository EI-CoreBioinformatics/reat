version development

import "structs.wdl"

workflow bam2hints {
	input {
		Array[File]? single_seqs
		Array[File]? many_seqs
		IndexedBAM? bam
		String? dUTP
		String output_prefix
	}

	if (defined(single_seqs)) {
		# Run chunk by chunk
		scatter (sequence_file in select_first([single_seqs])) {
			Array[String] sequences = read_lines(sequence_file)
			call SingleSeq {
				input:
				seq = sequences[0],
				bam = bam,
				dUTP = dUTP,
				output_prefix = "single"
			}
		}
	}

	if (defined(many_seqs)) {
		scatter (sequences in select_first([many_seqs])) {
			Array[String] lines = read_lines(sequences)
			call ManySeq {
				input:
				seqs = lines,
				bam = bam,
				dUTP = dUTP,
				output_prefix = "multi"
			}
		}
	}

	call JoinedHints {
		input:
		files = flatten(select_all([SingleSeq.hints, ManySeq.hints]))
	}

	output {
		File expression_gff = JoinedHints.hints
	}
}

task JoinedHints {
	input {
		Array[File] files
	}
	output {
		File hints = "exon.hints.gff"
	}
	command <<<
		cat ~{sep=" " files} > exon.hints.gff
	>>>
}

task SingleSeq {
	input {
		String seq
		IndexedBAM? bam
		String? dUTP
		String output_prefix
	}

	output {
		File hints = output_prefix + ".exonhints.gff"
	}

	command <<<
		touch unstranded.exonhints.augustus.gff \
		firststrand.Forward.exonhints.augustus.gff \
		firststrand.Reverse.exonhints.augustus.gff \
		secondstrand.Forward.exonhints.augustus.gff \
		secondstrand.Reverse.exonhints.augustus.gff

		if [ "~{dUTP}" == "secondstrand" ];
		then
				ln -s ~{select_first([bam]).bam}
				ln -s ~{select_first([bam]).index}
				samtools depth <(samtools view -b -f 128 -F 16 ~{basename(select_first([bam]).bam)} ~{seq}) <(samtools view -b -f 80 ~{basename(select_first([bam]).bam)} ~{seq}) | \
				awk -v OFS='\t' '
					BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
					{
					if (chrm != $1) {chrm=$1;}
					total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
					if (total_coverage < 10) next;
					if (start=="undef") start=$2;
					if ($2-start > 10 || length(values)>10) {
						tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
						print chrm, "w2h_secondstrand", "exonpart", start-1, prev_pos+1, ftot, "+", ".", "src=generic_source;pri=0;mult="int(ftot);
						start="undef";
						split("", values);
					} prev_pos = $2; values[$2] = total_coverage;}' > secondstrand.Forward.exonhints.augustus.gff &

				samtools depth <(samtools view -b -f 144 ~{basename(select_first([bam]).bam)} ~{seq}) <(samtools view -b -f 64 -F 16 ~{basename(select_first([bam]).bam)} ~{seq}) | \
				awk -v OFS='\t' '
					BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
					{
					if (chrm != $1) {chrm=$1;}
					total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
					if (total_coverage < 10) next;
					if (start=="undef") start=$2;
					if ($2-start > 10 || length(values)>10) {
						tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
						print chrm, "w2h_secondstrand", "exonpart", start-1, prev_pos+1, ftot, "-", ".", "src=generic_source;pri=0;mult="int(ftot);
						start="undef";
						split("", values);
					} prev_pos = $2; values[$2] = total_coverage;}' > secondstrand.Reverse.exonhints.augustus.gff &
		fi

		if [ "~{dUTP}" == "firststrand" ];
		then
				ln -s ~{select_first([bam]).bam}
				ln -s ~{select_first([bam]).index}
				samtools depth <(samtools view -b -f 128 -F 16 ~{basename(select_first([bam]).bam)} ~{seq}) <(samtools view -b -f 80 ~{basename(select_first([bam]).bam)} ~{seq}) | \
				awk -v OFS='\t' '
					BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
					{
					if (chrm != $1) {chrm=$1;}
					total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
					if (total_coverage < 10) next;
					if (start=="undef") start=$2;
					if ($2-start > 10 || length(values)>10) {
						tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
						print chrm, "w2h_firststrand", "exonpart", start-1, prev_pos+1, ftot, "+", ".", "src=generic_source;pri=0;mult="int(ftot);
						start="undef";
						split("", values);
					} prev_pos = $2; values[$2] = total_coverage;}' > firststrand.Reverse.exonhints.augustus.gff &

				samtools depth <(samtools view -b -f 144 ~{basename(select_first([bam]).bam)} ~{seq}) <(samtools view -b -f 64 -F 16 ~{basename(select_first([bam]).bam)} ~{seq}) | \
				awk -v OFS='\t' '
					BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
					{
					if (chrm != $1) {chrm=$1;}
					total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
					if (total_coverage < 10) next;
					if (start=="undef") start=$2;
					if ($2-start > 10 || length(values)>10) {
						tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
						print chrm, "w2h_firststrand", "exonpart", start-1, prev_pos+1, ftot, "-", ".", "src=generic_source;pri=0;mult="int(ftot);
						start="undef";
						split("", values);
					} prev_pos = $2; values[$2] = total_coverage;}' > firststrand.Forward.exonhints.augustus.gff &
		fi

		if [ "~{dUTP}" == "unstranded" ];
		then
			ln -s ~{select_first([bam]).bam}
			ln -s ~{select_first([bam]).index}
			samtools depth <(samtools view -b ~{basename(select_first([bam]).bam)} ~{seq}) | \
			awk -v OFS='\t' '
				BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
				{
				if (chrm != $1) {chrm=$1;}
				total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
				if (total_coverage < 10) next;
				if (start=="undef") start=$2;
				if ($2-start > 10 || length(values)>10) {
					tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
					print chrm, "w2h_unstranded", "exonpart", start-1, prev_pos+1, ftot, ".", ".", "src=generic_source;pri=0;mult="int(ftot);
					start="undef";
					split("", values);
				} prev_pos = $2; values[$2] = total_coverage;}' > unstranded.exonhints.augustus.gff &
		fi

		wait

		cat unstranded.exonhints.augustus.gff  \
		secondstrand.Forward.exonhints.augustus.gff firststrand.Forward.exonhints.augustus.gff \
		secondstrand.Reverse.exonhints.augustus.gff firststrand.Reverse.exonhints.augustus.gff > ~{output_prefix}.exonhints.gff
	>>>

}

task ManySeq {
	input {
		Array[String] seqs
		IndexedBAM? bam
		String? dUTP
		String output_prefix
	}

	output {
		File hints = output_prefix + ".exonhints.gff"
	}

	command <<<

		touch unstranded.exonhints.augustus.gff \
		firststrand.Forward.exonhints.augustus.gff \
		firststrand.Reverse.exonhints.augustus.gff \
		secondstrand.Forward.exonhints.augustus.gff \
		secondstrand.Reverse.exonhints.augustus.gff

		for name in ~{sep=" " seqs}; do
			if [ "~{dUTP}" == "secondstrand" ];
			then
					ln -s ~{select_first([bam]).bam}
					ln -s ~{select_first([bam]).index}
					samtools depth <(samtools view -b -f 128 -F 16 ~{basename(select_first([bam]).bam)} $name) <(samtools view -b -f 80 ~{basename(select_first([bam]).bam)} $name) | \
					awk -v OFS='\t' '
						BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
						{
						if (chrm != $1) {chrm=$1;}
						total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
						if (total_coverage < 10) next;
						if (start=="undef") start=$2;
						if ($2-start > 10 || length(values)>10) {
							tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
							print chrm, "w2h_secondstrand", "exonpart", start-1, prev_pos+1, ftot, "+", ".", "src=generic_source;pri=0;mult="int(ftot);
							start="undef";
							split("", values);
						} prev_pos = $2; values[$2] = total_coverage;}' > secondstrand.Forward.exonhints.augustus.gff &

					samtools depth <(samtools view -b -f 144 ~{basename(select_first([bam]).bam)} $name) <(samtools view -b -f 64 -F 16 ~{basename(select_first([bam]).bam)} $name) | \
					awk -v OFS='\t' '
						BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
						{
						if (chrm != $1) {chrm=$1;}
						total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
						if (total_coverage < 10) next;
						if (start=="undef") start=$2;
						if ($2-start > 10 || length(values)>10) {
							tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
							print chrm, "w2h_secondstrand", "exonpart", start-1, prev_pos+1, ftot, "-", ".", "src=generic_source;pri=0;mult="int(ftot);
							start="undef";
							split("", values);
						} prev_pos = $2; values[$2] = total_coverage;}' > secondstrand.Reverse.exonhints.augustus.gff &
			fi

			if [ "~{dUTP}" == "firststrand" ];
			then
					ln -s ~{select_first([bam]).bam}
					ln -s ~{select_first([bam]).index}
					samtools depth <(samtools view -b -f 128 -F 16 ~{basename(select_first([bam]).bam)} $name) <(samtools view -b -f 80 ~{basename(select_first([bam]).bam)} $name) | \
					awk -v OFS='\t' '
						BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
						{
						if (chrm != $1) {chrm=$1;}
						total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
						if (total_coverage < 10) next;
						if (start=="undef") start=$2;
						if ($2-start > 10 || length(values)>10) {
							tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
							print chrm, "w2h_firststrand", "exonpart", start-1, prev_pos+1, ftot, "+", ".", "src=generic_source;pri=0;mult="int(ftot);
							start="undef";
							split("", values);
						} prev_pos = $2; values[$2] = total_coverage;}' > firststrand.Reverse.exonhints.augustus.gff &

					samtools depth <(samtools view -b -f 144 ~{basename(select_first([bam]).bam)} $name) <(samtools view -b -f 64 -F 16 ~{basename(select_first([bam]).bam)} $name) | \
					awk -v OFS='\t' '
						BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
						{
						if (chrm != $1) {chrm=$1;}
						total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
						if (total_coverage < 10) next;
						if (start=="undef") start=$2;
						if ($2-start > 10 || length(values)>10) {
							tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
							print chrm, "w2h_firststrand", "exonpart", start-1, prev_pos+1, ftot, "-", ".", "src=generic_source;pri=0;mult="int(ftot);
							start="undef";
							split("", values);
						} prev_pos = $2; values[$2] = total_coverage;}' > firststrand.Forward.exonhints.augustus.gff &
			fi

			if [ "~{dUTP}" == "unstranded" ];
			then
				ln -s ~{select_first([bam]).bam}
				ln -s ~{select_first([bam]).index}
				samtools depth <(samtools view -b ~{basename(select_first([bam]).bam)} $name) | \
				awk -v OFS='\t' '
					BEGIN{chrm=""; split("", values); start="undef"; prev_pos=0}
					{
					if (chrm != $1) {chrm=$1;}
					total_coverage=0; for (i=3; i <= NF; i++) total_coverage += $i;
					if (total_coverage < 10) next;
					if (start=="undef") start=$2;
					if ($2-start > 10 || length(values)>10) {
						tot=0; for (v in values) tot += values[v]; ftot = tot/length(values);
						print chrm, "w2h_unstranded", "exonpart", start-1, prev_pos+1, ftot, ".", ".", "src=generic_source;pri=0;mult="int(ftot);
						start="undef";
						split("", values);
					} prev_pos = $2; values[$2] = total_coverage;}' > unstranded.exonhints.augustus.gff &
			fi

			wait

			cat unstranded.exonhints.augustus.gff  \
			secondstrand.Forward.exonhints.augustus.gff firststrand.Forward.exonhints.augustus.gff \
			secondstrand.Reverse.exonhints.augustus.gff firststrand.Reverse.exonhints.augustus.gff >> ~{output_prefix}.exonhints.gff
		done
	>>>

}