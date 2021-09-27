import os


def symlink(path, out_file):
    if os.path.exists(os.path.join(path, os.path.basename(out_file))):
        os.unlink(os.path.join(path, os.path.basename(out_file)))
    os.symlink(out_file, os.path.join(path, os.path.basename(out_file)))


def link_mikado(outputs, outputs_path, mikado):
    for suffix in ('loci', 'scores', 'metrics', 'stats'):
        mikado_name = mikado + suffix
        if outputs[mikado_name]:
            symlink(outputs_path, outputs[mikado_name])


def link_assemblies(assemblies, assembly_path, assembly_stats, outputs):
    if outputs[assemblies]:
        if not os.path.exists(assembly_path):
            os.mkdir(assembly_path)
        if outputs[assembly_stats]:
            for stats in outputs[assembly_stats]:
                symlink(assembly_path, stats)
        for sr_asm in outputs[assemblies]:
            symlink(assembly_path, sr_asm['assembly'])


def link_bams(outputs, outputs_path, bams_array, stats_array, stats_table):
    alignments_path = os.path.join(outputs_path, 'alignments')
    if outputs[bams_array]:
        if not os.path.exists(alignments_path):
            os.mkdir(alignments_path)
        for aligned_sample in outputs[bams_array]:
            for bam in aligned_sample['bam']:
                symlink(alignments_path, bam)
    if outputs[stats_array]:
        for stats in outputs[stats_array]:
            symlink(alignments_path, stats)
    if outputs[stats_table]:
        symlink(alignments_path, outputs[stats_table])