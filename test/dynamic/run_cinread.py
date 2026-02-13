import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_cinread(biscuit_dir, out_dir, prefix, ref_path, bam_dir, force):
    """Generate cinread files."""
    if not os.path.exists(f'{bam_dir}/{prefix}.bam') or not os.path.exists(f'{bam_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{out_dir}/{prefix}.cinread'):
        if not force:
            logger.info(f'Found {prefix} cinread files in {out_dir}. Rerun with `force.cinread = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} cinread files in {out_dir}, but `force.cinread = true` - REGENERATING cinread files')

    logger.info(f'Running {prefix} BISCUIT cinread')

    cmd = f'{biscuit_dir}/biscuit cinread {ref_path} {bam_dir}/{prefix}.bam'
    logger.debug(f'{prefix} BISCUIT cinread command: {cmd}')
    with open(f'{out_dir}/{prefix}.cinread', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, bam_dir, force):
    logger.info('Starting cinread testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_cinread(biscuit_dir, f'{out_dir}', 'new', ref_path, bam_dir, force)

    if compare_files.compare_files('.cinread', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
        logger.info(f'*.cinread match')
    else:
        diffs = compare_files.compare_line_by_line('.cinread', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
        for diff in diffs:
            idx, l_current, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.cinread - see above for differences')
            sys.exit(1)

    return None
