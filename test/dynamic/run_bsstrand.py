import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_bsstrand(biscuit_path, dir, tag, bam_dir, ref, force):
    """Generate bsstrand files."""
    if not os.path.exists(f'{bam_dir}/{tag}.bam') or not os.path.exists(f'{bam_dir}/{tag}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    EXTS = ['.bss.bam', '.bss.sam', '.bss']
    if all([os.path.exists(f'{dir}/{tag}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {tag} bsstrand files in {dir}. Rerun with `force.bsstrand = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} bsstrand files in {dir}, but `force.bsstrand = true` - REGENERATING bsstrand files')

    logger.info(f'Running {tag} BISCUIT bsstrand')

    cmd1 = f'{biscuit_path}/biscuit bsstrand {ref} {bam_dir}/{tag}.bam {dir}/{tag}.bss.bam'
    logger.debug(f'{tag} BISCUIT bsstrand command: {cmd1}')
    with open(f'{dir}/{tag}.bss', 'w') as f:
        subprocess.run(cmd1.split(' '), stdout=subprocess.DEVNULL, stderr=f)

    cmd2 = f'samtools view -h -o {dir}/{tag}.bss.sam -O SAM {dir}/{tag}.bss.bam'
    logger.debug(f'{tag} BISCUIT bsstrand command: {cmd2}')
    subprocess.run(cmd2.split(' '), stderr=subprocess.DEVNULL)

    return None

def main(dir, bam_dir, old, new, ref, force):
    logger.info('Starting bsstrand testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_bsstrand(old, f'{dir}', 'old', bam_dir, ref, force)
    run_bsstrand(new, f'{dir}', 'new', bam_dir, ref, force)

    for ext in ['.bss.sam', '.bss']:
        if compare_files.compare_files(ext, f'{dir}/old', f'{dir}/new'):
            logger.info(f'*.{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'{dir}/old', f'{dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_old, l_new = diff
                if ext == '.bss.sam' and l_old.startswith('@PG') and l_new.startswith('@PG'):
                    continue
                elif ext == '.bss' and l_old.startswith('[main]') and l_new.startswith('[main]'):
                    continue
                else:
                    n_diffs += 1
                    print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            elif n_diffs == 0 and ext == '.bss.sam':
                logger.warning(f'Mismatch only in @PG tag(s) in files: *{ext}')
            elif n_diffs == 0 and ext == '.bss':
                logger.warning(f'Mismatch only in [main] line(s) in files: *{ext}')

    return None
