import logging
import tomllib
import sys
import os

import run_index
import run_align
import run_pileup
import run_vcf2bed
import run_mergecg
import run_bsconv
import run_bsstrand
import run_cinread
import run_tview

# Define logging for entire program here
def setup_logger():
    """
    Setup logging for all tests. All modules should run

    import logging
    logger = logging.getLogger(__name__)

    at the top of its file to become children of this main logger
    """
    FORMAT = "[{levelname:<7}] {asctime} - {name:<12} :: {funcName:<15} - {message}"
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
    REF = '../data/ref/chr22.fa.gz'
    if not check_path(REF) or not check_path(f'{REF}.fai'):
        print('Reference FASTA missing. Please move up a directory and run `setup_tests.py` to retrieve.')
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
    logger.info(f'New BISCUIT path: {NEW}')
    for outer_key, dic in conf.items():
        for inner_key, value in dic.items():
            logger.info(f'Runtime configuration: {outer_key}.{inner_key} = {value}')

    if conf['run']['index']:
        run_index.main(NEW, '00_index', REF, conf['force']['index'])
    if conf['run']['align']:
        run_align.main(NEW, '01_align', '00_index', conf['force']['align'])
    if conf['run']['pileup']:
        run_pileup.main(NEW, '02_pileup', REF, '01_align', conf['force']['pileup'])
    if conf['run']['vcf2bed']:
        run_vcf2bed.main(NEW, '03_vcf2bed', '02_pileup', conf['force']['vcf2bed'])
    if conf['run']['mergecg']:
        run_mergecg.main(NEW, '04_mergecg', REF, '03_vcf2bed', conf['force']['mergecg'])
    if conf['run']['bsconv']:
        run_bsconv.main(NEW, '05_bsconv', REF, '01_align', conf['force']['bsconv'])
    if conf['run']['bsstrand']:
        run_bsstrand.main(NEW, '06_bsstrand', REF, '01_align', conf['force']['bsstrand'])
    if conf['run']['cinread']:
        run_cinread.main(NEW, '07_cinread', REF, '01_align', conf['force']['cinread'])
    if conf['run']['tview']:
        run_tview.main(NEW, 'XX_tview', REF, '01_align', conf['force']['tview'])

    logger.info('Finished testing')

    return None

if __name__ == '__main__':
    main()
