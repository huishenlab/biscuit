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
        MODULES = [
            'run_index',
            'run_align',
            'check_alignments',
            'run_pileup',
            'run_vcf2bed',
            'run_mergecg',
            'run_bsconv',
            'run_bsstrand',
            'run_cinread',
            'run_qc',
            'run_qc_coverage',
            'run_tview',
        ]
        for mod in MODULES:
            logging.getLogger(mod).setLevel(logging.DEBUG)

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

    # Index only needs to be run one time
    if conf['run']['index']:
        run_index.main(NEW, '000_index', REF, conf['force']['index'])

    # Loop through test types and run any requested tests
    type_to_fqs = {
        'wgbs_pe': '../data/sim_pe_wgbs_1.fastq.gz ../data/sim_pe_wgbs_2.fastq.gz',
        'wgbs_se': '../data/sim_se_wgbs.fastq.gz',
        'wgs_pe': '../data/sim_pe_wgs_1.fastq.gz ../data/sim_pe_wgs_2.fastq.gz',
        'wgs_se': '../data/sim_se_wgs.fastq.gz',
    }
    for key in conf['type']:
        if conf['type'][key]:
            logger.info(f'Starting tests for test type: {key}')
            if conf['run']['align']:
                run_align.main(
                    NEW, f'{key}_00_align', '000_index', type_to_fqs[key], 'se' in key, conf['force']['align']
                )
            if conf['run']['pileup']:
                run_pileup.main(
                    NEW, f'{key}_01_pileup', REF, f'{key}_00_align', conf['force']['align']
                )
            if conf['run']['vcf2bed']:
                run_vcf2bed.main(
                    NEW, f'{key}_02_vcf2bed', f'{key}_01_pileup', conf['force']['vcf2bed']
                )
            if conf['run']['mergecg']:
                run_mergecg.main(
                    NEW, f'{key}_03_mergecg', REF, f'{key}_02_vcf2bed', conf['force']['mergecg']
                )
            if conf['run']['bsconv']:
                run_bsconv.main(
                    NEW, f'{key}_04_bsconv', REF, f'{key}_00_align', conf['force']['bsconv']
                )
            if conf['run']['bsstrand']:
                run_bsstrand.main(
                    NEW, f'{key}_05_bsstrand', REF, f'{key}_00_align', conf['force']['bsstrand']
                )
            if conf['run']['cinread']:
                run_cinread.main(
                    NEW, f'{key}_06_cinread', REF, f'{key}_00_align', conf['force']['cinread']
                )
            if conf['run']['qc']:
                run_qc.main(
                    NEW, f'{key}_07_qc', REF, f'{key}_00_align', 'se' in key, conf['force']['qc']
                )
            if conf['run']['qc_coverage']:
                run_qc_coverage.main(
                    NEW, f'{key}_08_qc_coverage', REF, BOT, TOP, CPG, f'{key}_00_align', conf['force']['qc_coverage']
                )
            if conf['run']['tview']:
                run_tview.main(
                    NEW, f'{key}_XX_tview', REF, f'{key}_00_align', conf['force']['tview']
                )
        else:
            logger.info(f'Skipping tests for test type: {key}')

    logger.info('Finished testing')

    return None

if __name__ == '__main__':
    main()
