import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def check_index_exists(idx_path, prefix):
    """Checks for necessary input files"""
    for ext in ['.bis.amb', '.bis.ann', '.bis.pac', '.dau.bwt', '.dau.sa', '.par.bwt', '.par.sa']:
        if not os.path.exists(f'{idx_path}/{prefix}{ext}'):
            return False

    return True

def run_align(biscuit_dir, out_dir, prefix, idx_path, fq_dir, force):
    if not check_index_exists(idx_path, prefix):
        logger.error('Missing index files. Rerun with `run.index = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    EXTS = ['.sam', '.bam', '.bam.csi', '.debug']
    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} alignment files in {out_dir}. Rerun with `force.align = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} alignment files in {out_dir}, but `force.align = true` - REGENERATING alignment files')

    logger.info(f'Running {prefix} BISCUIT alignment')

    # Basic BISCUIT alignment output
    cmd1 = f'{biscuit_dir}/biscuit align {idx_path}/{prefix} {fq_dir}/simulated_1.fastq.gz {fq_dir}/simulated_2.fastq.gz'
    logger.debug(f'{prefix} BISCUIT alignment command: {cmd1}')
    with open(f'{out_dir}/{prefix}.sam', 'w') as f:
        subprocess.run(cmd1.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    # Verbose BISCUIT alignment output
    cmd2 = f'{biscuit_dir}/biscuit align -v 4 {idx_path}/{prefix} {fq_dir}/simulated_1.fastq.gz {fq_dir}/simulated_2.fastq.gz'
    logger.debug(f'{prefix} BISCUIT alignment command: {cmd2}')
    with open(f'{out_dir}/{prefix}.debug', 'w') as f:
        subprocess.run(cmd2.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    # Sort and index basic BISCUIT alignment output (used for downstream tests)
    cmd3 = f'samtools sort -o {out_dir}/{prefix}.bam -O BAM --write-index {out_dir}/{prefix}.sam'
    logger.debug(f'{prefix} BISCUIT alignment command: {cmd3}')
    subprocess.run(cmd3.split(' '), stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, idx_path, force):
    logger.info('Starting align testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_align(biscuit_dir, out_dir, 'new', idx_path, '../data', force)

    for ext in ['.sam', '.debug']:
        if compare_files.compare_files(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
            logger.info(f'*{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_current, l_new = diff
                if l_current.startswith('@PG') and l_new.startswith('@PG'):
                    continue
                else:
                    n_diffs += 1
                    print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            else:
                logger.warning(f'Mismatch only in @PG tag(s) in files: *{ext}')

    return None
