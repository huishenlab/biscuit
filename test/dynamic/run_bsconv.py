import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_bsconv(biscuit_dir, out_dir, prefix, ref_path, bam_dir, force):
    """Generate bisulfite conversion files."""
    if not os.path.exists(f'{bam_dir}/{prefix}.bam') or not os.path.exists(f'{bam_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{out_dir}/{prefix}.bsconv'):
        if not force:
            logger.info(f'Found {prefix} bsconv files in {out_dir}. Rerun with `force.bsconv = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} bsconv files in {out_dir}, but `force.bsconv = true` - REGENERATING bsconv files')

    logger.info(f'Running {prefix} BISCUIT bsconv')

    cmd = f'{biscuit_dir}/biscuit bsconv {ref_path} {bam_dir}/{prefix}.bam'
    logger.debug(f'{prefix} BISCUIT bsconv command: {cmd}')
    with open(f'{out_dir}/{prefix}.bsconv', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, bam_dir, force):
    logger.info('Starting bsconv testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_bsconv(biscuit_dir, out_dir, 'new', ref_path, bam_dir, force)

    if compare_files.compare_files('.bsconv', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
        logger.info(f'*.bsconv match')
    else:
        diffs = compare_files.compare_line_by_line('.bsconv', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
        n_diffs = 0
        for diff in diffs:
            idx, l_current, l_new = diff
            if l_current.startswith('@PG') and l_new.startswith('@PG'):
                continue
            else:
                n_diffs += 1
                print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

        if n_diffs > 0:
            logger.error(f'Mismatch in files: *.bsconv - see above for differences')
            sys.exit(1)
        else:
            logger.warning(f'Mismatch only in @PG tag(s) in files: *.bsconv')

    return None
