import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

EXTS = [
    '_covdist_all_base_botgc_table.txt',
    '_covdist_all_base_table.txt',
    '_covdist_all_base_topgc_table.txt',
    '_covdist_all_cpg_botgc_table.txt',
    '_covdist_all_cpg_table.txt',
    '_covdist_all_cpg_topgc_table.txt',
    '_covdist_q40_base_botgc_table.txt',
    '_covdist_q40_base_table.txt',
    '_covdist_q40_base_topgc_table.txt',
    '_covdist_q40_cpg_botgc_table.txt',
    '_covdist_q40_cpg_table.txt',
    '_covdist_q40_cpg_topgc_table.txt',
    '_cv_table.txt',
    '_filtered_read_counts.txt',
]

def run_qc_coverage(biscuit_dir, out_dir, prefix, ref_path, bot_path, top_path, cpg_path, bam_dir, force):
    """Generate qc_coverage files."""
    if not os.path.exists(f'{bam_dir}/{prefix}.bam') or not os.path.exists(f'{bam_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} qc_coverage files in {out_dir}. Rerun with `force.qc_coverage = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} qc_coverage files in {out_dir}, but `force.qc_coverage = true` - REGENERATING qc_coverage files')

    logger.info(f'Running {prefix} BISCUIT qc_coverage')

    cmd = f'{biscuit_dir}/biscuit qc_coverage -P {out_dir}/{prefix} -B {bot_path} -T {top_path} {ref_path} {cpg_path} {bam_dir}/{prefix}.bam'
    logger.debug(f'{prefix} BISCUIT qc_coverage command: {cmd}')
    with open(f'{out_dir}/{prefix}.err', 'w') as err:
        subprocess.run(cmd.split(' '), stderr=err)

    return None

def main(biscuit_dir, out_dir, ref_path, bot_path, top_path, cpg_path, bam_dir, force):
    logger.info('Starting qc_coverage testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_qc_coverage(biscuit_dir, out_dir, 'new', ref_path, bot_path, top_path, cpg_path, bam_dir, force)

    for ext in EXTS:
        if compare_files.compare_files(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
            logger.info(f'*{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
            for diff in diffs:
                idx, l_current, l_new = diff
                logger.debug(f'line {idx} --- OLD --- {l_current}')
                logger.debug(f'line {idx} --- NEW --- {l_new}')

            if logger.level == logging.DEBUG:
                logger.error(f'Mismatch in files: *{ext} - see above to see differences')
            else:
                logger.error(f'Mismatch in files: *{ext} - set `verbose = true` in config.toml to see differences')
            sys.exit(1)

    return None
