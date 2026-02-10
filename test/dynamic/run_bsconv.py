import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_bsconv(biscuit_path, dir, tag, bam_dir, ref, force):
    """Generate bisulfite conversion files."""
    if not os.path.exists(f'{bam_dir}/{tag}.bam') or not os.path.exists(f'{bam_dir}/{tag}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{dir}/{tag}.bsconv'):
        if not force:
            logger.info(f'Found {tag} bsconv files in {dir}. Rerun with `force.bsconv = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} bsconv files in {dir}, but `force.bsconv = true` - REGENERATING bsconv files')

    logger.info(f'Running {tag} BISCUIT bsconv')

    cmd = f'{biscuit_path}/biscuit bsconv {ref} {bam_dir}/{tag}.bam'
    logger.debug(f'{tag} BISCUIT bsconv command: {cmd}')
    with open(f'{dir}/{tag}.bsconv', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(dir, bam_dir, old, new, ref, force):
    logger.info('Starting bsconv testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_bsconv(old, f'{dir}', 'old', bam_dir, ref, force)
    run_bsconv(new, f'{dir}', 'new', bam_dir, ref, force)

    if compare_files.compare_files('.bsconv', f'{dir}/old', f'{dir}/new'):
        logger.info(f'*.bsconv match')
    else:
        diffs = compare_files.compare_line_by_line('.bsconv', f'{dir}/old', f'{dir}/new')
        n_diffs = 0
        for diff in diffs:
            idx, l_old, l_new = diff
            if l_old.startswith('@PG') and l_new.startswith('@PG'):
                continue
            else:
                n_diffs += 1
                print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

        if n_diffs > 0:
            logger.error(f'Mismatch in files: *.bsconv - see above for differences')
            sys.exit(1)
        else:
            logger.warning(f'Mismatch only in @PG tag(s) in files: *.bsconv')

    return None
