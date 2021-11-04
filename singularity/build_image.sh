set -euxo
version=0.4.6
rundir=$(dirname "$(realpath "$0")")
cd "$(mktemp -d)"
cp "${rundir}"/reat_singularity.def reat.def
sudo singularity build reat.img reat.def
mkdir -p /ei/software/testing/reat/${version}/x86_64/bin
cp reat.img /ei/software/testing/reat/${version}/x86_64/reat-${version}.img

cat << 'EOF' > /ei/software/testing/reat/${version}/x86_64/bin/singularity.exec
#!/bin/bash
DIR=$(dirname "$(readlink -f "$0")")
singularity exec "$DIR"/../reat-${version}.img $(basename "$0") $@
EOF

pushd /ei/software/testing/reat/${version}/x86_64/bin

ln -s singularity.exec hisat2
ln -s singularity.exec stringtie
ln -s singularity.exec samtools
ln -s singularity.exec gmap
ln -s singularity.exec gmapl
ln -s singularity.exec gsnap
ln -s singularity.exec gsnapl
ln -s singularity.exec gffread
ln -s singularity.exec minimap2
ln -s singularity.exec seqtk
ln -s singularity.exec spaln
ln -s singularity.exec sortgrcd
ln -s singularity.exec STAR
ln -s singularity.exec STARlong
ln -s singularity.exec TransDecoder.Predict
ln -s singularity.exec TransDecoder.LongOrfs
ln -s singularity.exec start_codon_refinement.py
ln -s singularity.exec seq_n_baseprobs_to_loglikelihood_vals.pl
ln -s singularity.exec seq_cache_populate.pl
ln -s singularity.exec select_best_ORFs_per_transcript.pl
ln -s singularity.exec score_CDS_likelihood_all_6_frames.pl
ln -s singularity.exec scallop
ln -s singularity.exec scallop-lr
ln -s singularity.exec train_start_PWM.pl
ln -s singularity.exec python
ln -s singularity.exec python3
ln -s singularity.exec prodigal
ln -s singularity.exec portcullis
ln -s singularity.exec plot-bamstats
ln -s singularity.exec junctools
ln -s singularity.exec get_top_longest_fasta_entries.pl
ln -s singularity.exec get_longest_ORF_per_transcript.pl
ln -s singularity.exec diamond
ln -s singularity.exec blastx
ln -s singularity.exec blastp
ln -s singularity.exec blastn
ln -s singularity.exec makeblastdb

popd