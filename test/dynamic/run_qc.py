import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

EXTS = [
    '_CpGRetentionByReadPos.txt',
    '_CpHRetentionByReadPos.txt',
    '_dup_report.txt',
    '_isize_table.txt',
    '_mapq_table.txt',
    '_strand_table.txt',
    '_totalReadConversionRate.txt',
]

def run_qc(biscuit_dir, out_dir, prefix, ref_path, bam_dir, single_end, force):
    """Generate qc files."""
    if not os.path.exists(f'{bam_dir}/{prefix}.bam') or not os.path.exists(f'{bam_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} qc files in {out_dir}. Rerun with `force.qc = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} qc files in {out_dir}, but `force.qc = true` - REGENERATING qc files')

    logger.info(f'Running {prefix} BISCUIT qc')

    cmd = f'{biscuit_dir}/biscuit qc {ref_path} {bam_dir}/{prefix}.bam {out_dir}/{prefix}'
    if single_end:
        cmd = f'{biscuit_dir}/biscuit qc -s {ref_path} {bam_dir}/{prefix}.bam {out_dir}/{prefix}'
    logger.debug(f'{prefix} BISCUIT qc command: {cmd}')
    with open(f'{out_dir}/{prefix}.err', 'w') as err:
        subprocess.run(cmd.split(' '), stderr=err)

    return None

def main(biscuit_dir, out_dir, ref_path, bam_dir, single_end, force):
    logger.info('Starting qc testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_qc(biscuit_dir, out_dir, 'new', ref_path, bam_dir, single_end, force)

    for ext in EXTS:
        if ext == '_isize_table.txt' and single_end:
            continue

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
