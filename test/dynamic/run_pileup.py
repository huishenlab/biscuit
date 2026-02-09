import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

EXTS = ['.vcf', '.vcf_meth_average.tsv']

def run_pileup(biscuit_path, dir, tag, ref, align_dir, force):
    if not os.path.exists(f'{align_dir}/{tag}.bam') or not os.path.exists(f'{align_dir}/{tag}.bam.csi'):
        logger.error('Missing BAM or BAM index. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if all([os.path.exists(f'{dir}/{tag}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {tag} pileup files in {dir}. Rerun with `force.pileup = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} pileup files in {dir}, but `force.pileup = true` - REGENERATING pileup files')

    logger.info(f'Running {tag} BISCUIT pileup')

    cmd = f'{biscuit_path}/biscuit pileup -o {dir}/{tag}.vcf {ref} {align_dir}/{tag}.bam'
    subprocess.run(cmd.split(' '))

    return None

def main(dir, ref, old, new, align_dir, force):
    logger.info('Starting pileup testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_pileup(old, f'{dir}', 'old', ref, align_dir, force)
    run_pileup(new, f'{dir}', 'new', ref, align_dir, force)

    for ext in EXTS:
        if compare_files.compare_files(ext, f'{dir}/old', f'{dir}/new'):
            logger.info(f'*{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'{dir}/old', f'{dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_old, l_new = diff

                if ext == '.vcf':
                    if l_old.startswith('##program') and l_new.startswith('##program'):
                        continue
                    elif l_old.startswith('#CHROM') and l_new.startswith('#CHROM'):
                        continue
                    else:
                        print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')
                        n_diffs += 1
                elif ext == '.vcf_meth_average.tsv':
                    r_old = l_old.replace('old', 'sample')
                    r_new = l_new.replace('new', 'sample')

                    if r_old == r_new:
                        continue
                    else:
                        print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')
                        n_diffs += 1


            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            elif ext == '.vcf_meth_average.tsv':
                logger.warning(f'Mismatch only in sample name(s) in files: *{ext}')
            elif ext == '.vcf':
                logger.warning(f'Mismatch only in program line and sample name(s) in files: *{ext}')

    return None
