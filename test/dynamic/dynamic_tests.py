import logging
import tomllib
import sys
import os

import run_index
import run_align
import check_alignments
import run_pileup
import run_vcf2bed
import run_mergecg
import run_bsconv
import run_bsstrand
import run_cinread
import run_qc
import run_qc_coverage
import run_tview
import run_se_align

# Define logging for entire program here
def setup_logger():
    """
    Setup logging for all tests. All modules should run

    import logging
    logger = logging.getLogger(__name__)

    at the top of its file to become children of this main logger
    """
    FORMAT = "[{levelname:<7}] {asctime} - {name:<16} :: {funcName:<15} - {message}"
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

def get_reference_files():
    REF = '../data/ref/chr22.fa.gz'
    if not check_path(REF) or not check_path(f'{REF}.fai'):
        logger.error('Reference FASTA missing. Please move up a directory and run `setup_tests.py` to retrieve.')
        sys.exit(1)

    BOT = '../data/ref/chr22.gc.bot10p.bed.gz'
    if not check_path(BOT) or not check_path(f'{BOT}.tbi'):
        logger.error('Bottom 10% GC content file missing. Please move up a directory and run `setup_tests.py` to retrieve.')
        sys.exit(1)

    TOP = '../data/ref/chr22.gc.top10p.bed.gz'
    if not check_path(TOP) or not check_path(f'{TOP}.tbi'):
        logger.error('Top 10% GC content file missing. Please move up a directory and run `setup_tests.py` to retrieve.')
        sys.exit(1)

    CPG = '../data/ref/chr22.cpg.bed.gz'
    if not check_path(CPG) or not check_path(f'{CPG}.tbi'):
        logger.error('CpG BED file missing. Please move up a directory and run `setup_tests.py` to retrieve.')
        sys.exit(1)

    return REF, BOT, TOP, CPG

def main():
    # Runtime configuration
    conf = read_config()
    if conf['verbose']:
        # Set top logging level
        logger.setLevel(logging.DEBUG)
        logger.debug('Exact mismatches will be shown')

        # Set imported module logging levels
        logging.getLogger('run_index').setLevel(logging.DEBUG)
        logging.getLogger('run_align').setLevel(logging.DEBUG)
        logging.getLogger('check_alignments').setLevel(logging.DEBUG)
        logging.getLogger('run_pileup').setLevel(logging.DEBUG)
        logging.getLogger('run_vcf2bed').setLevel(logging.DEBUG)
        logging.getLogger('run_mergecg').setLevel(logging.DEBUG)
        logging.getLogger('run_bsconv').setLevel(logging.DEBUG)
        logging.getLogger('run_bsstrand').setLevel(logging.DEBUG)
        logging.getLogger('run_cinread').setLevel(logging.DEBUG)
        logging.getLogger('run_qc').setLevel(logging.DEBUG)
        logging.getLogger('run_qc_coverage').setLevel(logging.DEBUG)
        logging.getLogger('run_tview').setLevel(logging.DEBUG)
        logging.getLogger('run_se_align').setLevel(logging.DEBUG)

    # Reference FASTA
    REF, BOT, TOP, CPG = get_reference_files()

    # New BISCUIT version
    NEW = None
    if check_path('../../bin') and check_path('../../bin/biscuit'):
        NEW = '../../bin'
    elif check_path('../../build/src') and check_path('../../build/src/biscuit'):
        NEW = '../../build/src'
    else:
        logger.error('Have you compiled BISCUIT yet?')
        sys.exit(1)

    logger.info(f'Reference path: {REF}')
    logger.info(f'New BISCUIT path: {NEW}')
    for outer_key, dic in conf.items():
        try:
            for inner_key, value in dic.items():
                logger.info(f'Runtime configuration: {outer_key}.{inner_key} = {value}')
        except AttributeError: # not a dictionary
            logger.info(f'Runtime configuration: {outer_key} = {dic}')

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
    if conf['run']['qc']:
        run_qc.main(NEW, '08_qc', REF, '01_align', conf['force']['qc'])
    if conf['run']['qc_coverage']:
        run_qc_coverage.main(NEW, '09_qc_coverage', REF, BOT, TOP, CPG, '01_align', conf['force']['qc_coverage'])
    if conf['run']['se_align']:
        run_se_align.main(NEW, '10_se_align', '00_index', conf['force']['se_align'])
    if conf['run']['se_pileup']:
        run_pileup.main(NEW, '11_se_pileup', REF, '10_se_align', conf['force']['se_pileup'])
    if conf['run']['se_vcf2bed']:
        run_vcf2bed.main(NEW, '12_se_vcf2bed', '11_se_pileup', conf['force']['se_vcf2bed'])
    if conf['run']['se_mergecg']:
        run_mergecg.main(NEW, '13_se_mergecg', REF, '12_se_vcf2bed', conf['force']['se_mergecg'])
    if conf['run']['se_bsconv']:
        run_bsconv.main(NEW, '14_se_bsconv', REF, '10_se_align', conf['force']['se_bsconv'])
    if conf['run']['se_bsstrand']:
        run_bsstrand.main(NEW, '15_se_bsstrand', REF, '10_se_align', conf['force']['se_bsstrand'])
    if conf['run']['se_cinread']:
        run_cinread.main(NEW, '16_se_cinread', REF, '10_se_align', conf['force']['se_cinread'])
    if conf['run']['se_qc_coverage']:
        run_qc_coverage.main(NEW, '18_se_qc_coverage', REF, BOT, TOP, CPG, '10_se_align', conf['force']['se_qc_coverage'])
    if conf['run']['tview']:
        run_tview.main(NEW, 'XX_tview', REF, '01_align', conf['force']['tview'])

    logger.info('Finished testing')

    return None

if __name__ == '__main__':
    main()
