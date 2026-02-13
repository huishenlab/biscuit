import subprocess
import logging
import tarfile
import urllib.request
import sys
import os

# Define logging for entire program here
def setup_logger():
    """
    Setup logging for all tests. All modules should run

    import logging
    logger = logging.getLogger(__name__)

    at the top of its file to become children of this main logger
    """
    FORMAT = "[{levelname:<7}] {asctime} - {name}::{funcName:<15} - {message}"
    logging.basicConfig(format=FORMAT, style="{", level=logging.INFO)

    return logging.getLogger(__name__)

logger = setup_logger()

def get_test_files(version):
    test_fname = 'biscuit_test_files.tar.gz'

    logger.info('Downloading BISCUIT test files')
    try:
        urllib.request.urlretrieve(
            url = f'https://github.com/huishenlab/biscuit/releases/download/{version}/{test_fname}',
            filename = test_fname
        )
    except Exception as e:
        logger.error(f'Problem download biscuit test files ({e})')
        sys.exit(1)
    logger.info('Finished downloading')

    logger.info('Extracting tarball')
    tar = tarfile.open(test_fname, mode='r:gz')
    tar.extractall()
    logger.info('Extraction finished')

    logger.info('Cleaning up tarball')
    try:
        os.remove(test_fname)
    except Exception as e:
        logger.error(f'Error deleting tarball ({e})')
        sys.exit(1)
    logger.info('Successfully cleaned up tarball')

    return None

def main():
    get_test_files('v1.7.1.20250908')

    return None

if __name__ == '__main__':
    main()
