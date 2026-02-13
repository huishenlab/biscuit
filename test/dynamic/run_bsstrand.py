import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

def run_bsstrand(biscuit_dir, out_dir, prefix, ref_path, bam_dir, force):
    """Generate bsstrand files."""
    if not os.path.exists(f'{bam_dir}/{prefix}.bam') or not os.path.exists(f'{bam_dir}/{prefix}.bam.csi'):
        logger.error('Missing BAM or BAM index file. Rerun with `run.align = true`')
        sys.exit(1)

    # If files exist and user doesn't force regeneration, skip processing
    EXTS = ['.bss.bam', '.bss.sam', '.bss']
    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} bsstrand files in {out_dir}. Rerun with `force.bsstrand = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} bsstrand files in {out_dir}, but `force.bsstrand = true` - REGENERATING bsstrand files')

    logger.info(f'Running {prefix} BISCUIT bsstrand')

    cmd1 = f'{biscuit_dir}/biscuit bsstrand {ref_path} {bam_dir}/{prefix}.bam {out_dir}/{prefix}.bss.bam'
    logger.debug(f'{prefix} BISCUIT bsstrand command: {cmd1}')
    with open(f'{out_dir}/{prefix}.bss', 'w') as f:
        subprocess.run(cmd1.split(' '), stdout=subprocess.DEVNULL, stderr=f)

    cmd2 = f'samtools view -h -o {out_dir}/{prefix}.bss.sam -O SAM {out_dir}/{prefix}.bss.bam'
    logger.debug(f'{prefix} BISCUIT bsstrand command: {cmd2}')
    subprocess.run(cmd2.split(' '), stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, bam_dir, force):
    logger.info('Starting bsstrand testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_bsstrand(biscuit_dir, f'{out_dir}', 'new', ref_path, bam_dir, force)

    for ext in ['.bss.sam', '.bss']:
        if compare_files.compare_files(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
            logger.info(f'*.{ext} match')
        else:
            diffs = compare_files.compare_line_by_line(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new')
            n_diffs = 0
            for diff in diffs:
                idx, l_current, l_new = diff
                if ext == '.bss.sam' and l_current.startswith('@PG') and l_new.startswith('@PG'):
                    continue
                elif ext == '.bss' and l_current.startswith('[main]') and l_new.startswith('[main]'):
                    continue
                else:
                    n_diffs += 1
                    print(f'line {idx}\n\tOLD -- {l_current}\n\tNEW -- {l_new}')

            if n_diffs > 0:
                logger.error(f'Mismatch in files: *{ext} - see above for differences')
                sys.exit(1)
            elif n_diffs == 0 and ext == '.bss.sam':
                logger.warning(f'Mismatch only in @PG tag(s) in files: *{ext}')
            elif n_diffs == 0 and ext == '.bss':
                logger.warning(f'Mismatch only in [main] line(s) in files: *{ext}')

    return None
