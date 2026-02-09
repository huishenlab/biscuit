import logging
import tomllib
import sys
import os

import run_index
import run_align
import run_pileup
import run_vcf2bed
import run_mergecg

# Define logging for entire program here
def setup_logger():
    """
    Setup logging for all tests. All modules should run

    import logging
    logger = logging.getLogger(__name__)

    at the top of its file to become children of this main logger
    """
    FORMAT = "[{levelname:<7}] {asctime} - {name:<15} - {message}"
    logging.basicConfig(format=FORMAT, style="{", level=logging.INFO)

    return logging.getLogger(__name__)

logger = setup_logger()

def check_path(path):
    """Wrapper to check existence of file or directory."""
    return os.path.exists(path)

def read_config():
    """Read configuration"""
    with open('config.toml', 'rb') as f:
        data = tomllib.load(f)

    return data

def main():
    # Reference FASTA
    REF = '../data/ref/chr1.fa.gz'
    if not check_path(REF) or not check_path(f'{REF}.fai'):
        print('Reference FASTA missing. Please run `get_fasta.sh` in ../data/ to retrieve.')
        sys.exit(1)

    # Previous BISCUIT version
    OLD = '../data/biscuit_master/bin'
    if not check_path(OLD) and not check_path(f'{OLD}/biscuit'):
        print('Original BISCUIT not downloaded and/or compiled. Please run `get_biscuit.sh` in ../data to retrieve.')
        sys.exit(1)

    # New BISCUIT version
    NEW = None
    if check_path('../../bin') and check_path('../../bin/biscuit'):
        NEW = '../../bin'
    elif check_path('../../build/src') and check_path('../../build/src/biscuit'):
        NEW = '../../build/src'
    else:
        print('Have you compiled BISCUIT yet?')
        sys.exit(1)

    # Runtime configuration
    conf = read_config()

    logger.info(f'Reference path: {REF}')
    logger.info(f'Old BISCUIT path: {OLD}')
    logger.info(f'New BISCUIT path: {NEW}')
    for outer_key, dic in conf.items():
        for inner_key, value in dic.items():
            logger.info(f'Runtime configuration: {outer_key}.{inner_key} = {value}')
    #logger.info(f'Runtime configuration: {conf}')

    if conf['run']['index']:
        run_index.main('00_index', REF, OLD, NEW, conf['force']['index'])
    if conf['run']['align']:
        run_align.main('01_align', '00_index', OLD, NEW, conf['force']['align'])
    if conf['run']['pileup']:
        run_pileup.main('02_pileup', REF, OLD, NEW, '01_align', conf['force']['pileup'])
    if conf['run']['vcf2bed']:
        run_vcf2bed.main('03_vcf2bed', '02_pileup', OLD, NEW, conf['force']['vcf2bed'])
    if conf['run']['mergecg']:
        run_mergecg.main('04_mergecg', '03_vcf2bed', OLD, NEW, REF, conf['force']['mergecg'])

    return None

if __name__ == '__main__':
    main()
