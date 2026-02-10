import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_cinread(biscuit_path, dir, tag, bam_dir, ref, force):
    """Generate cinread files."""
    if not os.path.exists(f'{bam_dir}/{tag}.bam') or not os.path.exists(f'{bam_dir}/{tag}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{dir}/{tag}.cinread'):
        if not force:
            logger.info(f'Found {tag} cinread files in {dir}. Rerun with `force.cinread = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} cinread files in {dir}, but `force.cinread = true` - REGENERATING cinread files')

    logger.info(f'Running {tag} BISCUIT cinread')

    cmd = f'{biscuit_path}/biscuit cinread {ref} {bam_dir}/{tag}.bam'
    logger.debug(f'{tag} BISCUIT cinread command: {cmd}')
    with open(f'{dir}/{tag}.cinread', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(dir, bed_dir, old, new, ref, force):
    logger.info('Starting cinread testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_cinread(old, f'{dir}', 'old', bed_dir, ref, force)
    run_cinread(new, f'{dir}', 'new', bed_dir, ref, force)

    if compare_files.compare_files('.cinread', f'{dir}/old', f'{dir}/new'):
        logger.info(f'*.cinread match')
    else:
        diffs = compare_files.compare_line_by_line('.cinread', f'{dir}/old', f'{dir}/new')
        for diff in diffs:
            idx, l_old, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.cinread - see above for differences')
            sys.exit(1)

    return None
