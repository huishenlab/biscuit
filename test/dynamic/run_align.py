import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def check_index_exists(idx_path, tag):
    """Checks for necessary input files"""
    for ext in ['.bis.amb', '.bis.ann', '.bis.pac', '.dau.bwt', '.dau.sa', '.par.bwt', '.par.sa']:
        if not os.path.exists(f'{idx_path}/{tag}{ext}'):
            return False

    return True

def run_align(biscuit_path, dir, tag, idx_path, fq_dir, force):
    if not check_index_exists(idx_path, tag):
        logger.error('Missing index files. Rerun with `run.index = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    EXTS = ['.sam', '.bam', '.bam.csi', '.debug']
    if all([os.path.exists(f'{dir}/{tag}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {tag} alignment files in {dir}. Rerun with `force.align = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} alignment files in {dir}, but `force.align = true` - REGENERATING alignment files')

    logger.info(f'Running {tag} BISCUIT alignment')

    # Basic BISCUIT alignment output
    cmd1 = f'{biscuit_path}/biscuit align {idx_path}/{tag} {fq_dir}/simulated_1.fastq.gz {fq_dir}/simulated_2.fastq.gz'
    logger.debug(f'{tag} BISCUIT alignment command: {cmd1}')
    with open(f'{dir}/{tag}.sam', 'w') as f:
        subprocess.run(cmd1.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    # Verbose BISCUIT alignment output
    cmd2 = f'{biscuit_path}/biscuit align -v 4 {idx_path}/{tag} {fq_dir}/simulated_1.fastq.gz {fq_dir}/simulated_2.fastq.gz'
    logger.debug(f'{tag} BISCUIT alignment command: {cmd2}')
    with open(f'{dir}/{tag}.debug', 'w') as f:
        subprocess.run(cmd2.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    # Sort and index basic BISCUIT alignment output (used for downstream tests)
    cmd3 = f'samtools sort -o {dir}/{tag}.bam -O BAM --write-index {dir}/{tag}.sam'
    logger.debug(f'{tag} BISCUIT alignment command: {cmd3}')
    subprocess.run(cmd3.split(' '), stderr=subprocess.DEVNULL)

    return None

def main(dir, idx_path, old, new, force):
    logger.info('Starting align testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_align(old, f'{dir}', 'old', idx_path, '../data', force)
    run_align(new, f'{dir}', 'new', idx_path, '../data', force)

    for ext in ['.sam', '.debug']:
        if compare_files.compare_files(ext, f'{dir}/old', f'{dir}/new'):
            logger.info(f'*{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'{dir}/old', f'{dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_old, l_new = diff
                if l_old.startswith('@PG') and l_new.startswith('@PG'):
                    continue
                else:
                    n_diffs += 1
                    print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            else:
                logger.warning(f'Mismatch only in @PG tag(s) in files: *{ext}')

    return None
