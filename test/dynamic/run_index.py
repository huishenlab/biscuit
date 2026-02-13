import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

# Index file extensions to compare
EXTS = ['.bis.amb', '.bis.ann', '.bis.pac', '.dau.bwt', '.dau.sa', '.par.bwt', '.par.sa']

def run_index(biscuit_dir, out_dir, prefix, ref_path, force):
    """Generate index files."""
    # If files exist and user doesn't force regeneration, skip processing
    if all([os.path.exists(f'{out_dir}/{prefix}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {prefix} index files in {out_dir}. Rerun with `force.index = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {prefix} index files in {out_dir}, but `force.index = true` - REGENERATING index files')

    logger.info(f'Running {prefix} BISCUIT indexing - this will take a while!')

    cmd = f'{biscuit_dir}/biscuit index -p {out_dir}/{prefix} {ref_path}'
    logger.debug(f'{prefix} BISCUIT indexing command: {cmd}')
    subprocess.run(cmd.split(' '), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return None

def main(biscuit_dir, out_dir, ref_path, force):
    logger.info('Starting index testing')

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run_index(biscuit_dir, out_dir, 'new', ref_path, force)

    for ext in EXTS:
        if compare_files.compare_files(ext, f'../data/dynamic/{out_dir}/current', f'{out_dir}/new'):
            logger.info(f'*{ext} match')
        else:
            logger.error(f'Mismatch in binary files: *{ext}')
            sys.exit(1)

    return None
