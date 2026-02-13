import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

EXTS = ['.vcf', '.vcf_meth_average.tsv']

def run_pileup(biscuit_dir, out_dir, prefix, ref_path, align_dir, force):
    if not os.path.exists(f'{align_dir}/{prefix}.bam') or not os.path.exists(f'{align_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} pileup files in {out_dir}. Rerun with `force.pileup = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} pileup files in {out_dir}, but `force.pileup = true` - REGENERATING pileup files')

    logger.info(f'Running {prefix} BISCUIT pileup')

    cmd = f'{biscuit_dir}/biscuit pileup -o {out_dir}/{prefix}.vcf {ref_path} {align_dir}/{prefix}.bam'
    subprocess.run(cmd.split(' '), stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, align_dir, force):
    logger.info('Starting pileup testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_pileup(biscuit_dir, out_dir, 'new', ref_path, align_dir, force)

    for ext in EXTS:
        if compare_files.compare_files(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
            logger.info(f'*{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_current, l_new = diff

                if ext == '.vcf':
                    if l_current.startswith('##program') and l_new.startswith('##program'):
                        continue
                    elif l_current.startswith('#CHROM') and l_new.startswith('#CHROM'):
                        continue
                    else:
                        print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')
                        n_diffs += 1
                elif ext == '.vcf_meth_average.tsv':
                    r_current = l_current.replace('current', 'sample')
                    r_new = l_new.replace('new', 'sample')

                    if r_current == r_new:
                        continue
                    else:
                        print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')
                        n_diffs += 1


            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            elif ext == '.vcf_meth_average.tsv':
                logger.warning(f'Mismatch only in sample name(s) in files: *{ext}')
            elif ext == '.vcf':
                logger.warning(f'Mismatch only in program line and sample name(s) in files: *{ext}')

    return None
