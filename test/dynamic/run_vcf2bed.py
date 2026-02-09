import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_vcf2bed(biscuit_path, dir, tag, vcf_dir, force):
    """Generate BISCUIT BED files."""
    if not os.path.exists(f'{vcf_dir}/{tag}.vcf'):
        logger.error('Missing VCF file. Rerun with `run.pileup = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{dir}/{tag}.bed'):
        if not force:
            logger.info(f'Found {tag} BISCUIT BED files in {dir}. Rerun with `force.vcf2bed = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} BISCUIT BED files in {dir}, but `force.vcf2bed = true` - REGENERATING BISCUIT BED files')

    logger.info(f'Running {tag} BISCUIT vcf2bed')

    cmd = f'{biscuit_path}/biscuit vcf2bed {vcf_dir}/{tag}.vcf'
    logger.debug(f'{tag} BISCUIT vcf2bed command: {cmd}')
    with open(f'{dir}/{tag}.bed', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(dir, vcf_dir, old, new, force):
    logger.info('Starting vcf2bed testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_vcf2bed(old, f'{dir}', 'old', vcf_dir, force)
    run_vcf2bed(new, f'{dir}', 'new', vcf_dir, force)

    if compare_files.compare_files('.bed', f'{dir}/old', f'{dir}/new'):
        logger.info(f'*.bed match')
    else:
        diffs = compare_files.compare_line_by_line('.bed', f'{dir}/old', f'{dir}/new')
        for diff in diffs:
            idx, l_old, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_old}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.bed - see above for differences')
            sys.exit(1)

    return None
