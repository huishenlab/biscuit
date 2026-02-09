import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_mergecg(biscuit_path, dir, tag, bed_dir, ref, force):
    """Generate BISCUIT BED files."""
    if not os.path.exists(f'{bed_dir}/{tag}.bed'):
        logger.error('Missing BISCUIT BED file. Rerun with `run.vcf2bed = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{dir}/{tag}.mergecg.bed'):
        if not force:
            logger.info(f'Found {tag} merged CG BED files in {dir}. Rerun with `force.mergecg = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} merged CG BED files in {dir}, but `force.mergecg = true` - REGENERATING merged CG BED files')

    logger.info(f'Running {tag} BISCUIT mergecg')

    cmd = f'{biscuit_path}/biscuit mergecg {ref} {bed_dir}/{tag}.bed'
    logger.debug(f'{tag} BISCUIT mergecg command: {cmd}')
    with open(f'{dir}/{tag}.mergecg.bed', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(dir, bed_dir, old, new, ref, force):
    logger.info('Starting mergecg testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_mergecg(old, f'{dir}', 'old', bed_dir, ref, force)
    run_mergecg(new, f'{dir}', 'new', bed_dir, ref, force)

    if compare_files.compare_files('.mergecg.bed', f'{dir}/old', f'{dir}/new'):
        logger.info(f'*.mergecg.bed match')
    else:
        diffs = compare_files.compare_line_by_line('.mergecg.bed', f'{dir}/old', f'{dir}/new')
        for diff in diffs:
            idx, l_old, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.mergecg.bed - see above for differences')
            sys.exit(1)

    return None
