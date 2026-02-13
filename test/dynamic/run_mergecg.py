import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_mergecg(biscuit_dir, out_dir, prefix, ref_path, bed_dir, force):
    """Generate BISCUIT BED files."""
    if not os.path.exists(f'{bed_dir}/{prefix}.bed'):
        logger.error('Missing BISCUIT BED file. Rerun with `run.vcf2bed = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{out_dir}/{prefix}.mergecg.bed'):
        if not force:
            logger.info(f'Found {prefix} merged CG BED files in {out_dir}. Rerun with `force.mergecg = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} merged CG BED files in {out_dir}, but `force.mergecg = true` - REGENERATING merged CG BED files')

    logger.info(f'Running {prefix} BISCUIT mergecg')

    cmd = f'{biscuit_dir}/biscuit mergecg {ref_path} {bed_dir}/{prefix}.bed'
    logger.debug(f'{prefix} BISCUIT mergecg command: {cmd}')
    with open(f'{out_dir}/{prefix}.mergecg.bed', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, bed_dir, force):
    logger.info('Starting mergecg testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_mergecg(biscuit_dir, out_dir, 'new', ref_path, bed_dir, force)

    if compare_files.compare_files('.mergecg.bed', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
        logger.info(f'*.mergecg.bed match')
    else:
        diffs = compare_files.compare_line_by_line('.mergecg.bed', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
        for diff in diffs:
            idx, l_current, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.mergecg.bed - see above for differences')
            sys.exit(1)

    return None
