import subprocess
import logging
import sys
import os

import compare_files

logger = logging.getLogger(__name__)

# Index file extensions to compare
EXTS = ['.bis.amb', '.bis.ann', '.bis.pac', '.dau.bwt', '.dau.sa', '.par.bwt', '.par.sa']

def run_index(biscuit_path, dir, tag, ref, force):
    """Generate index files."""
    # If files exist and user doesn't force regeneration, skip processing
    if all([os.path.exists(f'{dir}/{tag}{ext}') for ext in EXTS]):
        if not force:
            logger.info(f'Found {tag} index files in {dir}. Rerun with `force.index = true` to regenerate these files')
            return None
        else:
            logger.info(f'Found {tag} index files in {dir}, but `force.index = true` - REGENERATING index files')

    logger.info(f'Running {tag} BISCUIT indexing - this will take a while!')

    cmd = f'{biscuit_path}/biscuit index -p {dir}/{tag} {ref}'
    logger.debug(f'{tag} BISCUIT indexing command: {cmd}')
    subprocess.run(cmd.split(' '), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    return None

def main(dir, ref, old, new, force):
    logger.info('Starting index testing')

    if not os.path.exists(dir):
        os.makedirs(dir)

    run_index(old, f'{dir}', 'old', ref, force)
    run_index(new, f'{dir}', 'new', ref, force)

    for ext in EXTS:
        if compare_files.compare_files(ext, f'{dir}/old', f'{dir}/new'):
            logger.info(f'*{ext} match')
        else:
            logger.error(f'Mismatch in binary files: *{ext}')
            sys.exit(1)

    return None
