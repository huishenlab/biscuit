import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_vcf2bed(biscuit_dir, out_dir, prefix, vcf_dir, force):
    """Generate BISCUIT BED files."""
    if not os.path.exists(f'{vcf_dir}/{prefix}.vcf'):
        logger.error('Missing VCF file. Rerun with `run.pileup = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    if os.path.exists(f'{out_dir}/{prefix}.bed'):
        if not force:
            logger.info(f'Found {prefix} BISCUIT BED files in {out_dir}. Rerun with `force.vcf2bed = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} BISCUIT BED files in {out_dir}, but `force.vcf2bed = true` - REGENERATING BISCUIT BED files')

    logger.info(f'Running {prefix} BISCUIT vcf2bed')

    cmd = f'{biscuit_dir}/biscuit vcf2bed {vcf_dir}/{prefix}.vcf'
    logger.debug(f'{prefix} BISCUIT vcf2bed command: {cmd}')
    with open(f'{out_dir}/{prefix}.bed', 'w') as f:
        subprocess.run(cmd.split(' '), stdout=f, stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, vcf_dir, force):
    logger.info('Starting vcf2bed testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_vcf2bed(biscuit_dir, out_dir, 'new', vcf_dir, force)

    if compare_files.compare_files('.bed', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
        logger.info(f'*.bed match')
    else:
        diffs = compare_files.compare_line_by_line('.bed', f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
        for diff in diffs:
            idx, l_current, l_new = diff
            print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

            logger.error(f'Mismatch in files: *.bed - see above for differences')
            sys.exit(1)

    return None
