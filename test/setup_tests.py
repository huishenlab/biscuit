import subprocess
import logging
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

def get_fasta():
    # Make directory where the reference FASTA will live
    os.makedirs('data/ref', exist_ok=True)

    CHR_GZ = 'chr22.fa.gz'
    CHR = CHR_GZ.replace('.gz', '')

    # Retrieve FASTA
    logger.info(f'Downloading {CHR_GZ}')
    try:
        urllib.request.urlretrieve(
            url = f'https://hgdownload.soe.ucsc.edu/goldenpath/hg38/chromosomes/{CHR_GZ}',
            filename = f'data/ref/{CHR_GZ}'
        )
    except Exception as e:
        logging.error(f'Problem downloading {CHR_GZ} ({e})')
        sys.exit(1)
    logger.info('Finished downloading')

    # Decompress so that bgzip compression can happen
    logger.info('Decompressing for eventual bgzip compression')
    cmd1 = f'gunzip data/ref/{CHR_GZ}'
    subprocess.run(cmd1.split(' '), stderr=subprocess.DEVNULL)
    logger.info('Done with decompression')

    # bgzip compression
    logger.info('bgzip compressing FASTA')
    cmd2 = f'bgzip data/ref/{CHR}'
    subprocess.run(cmd2.split(' '), stderr=subprocess.DEVNULL)
    logger.info('Done with compression')

    # FASTA index
    logger.info('Indexing reference FASTA')
    cmd3 = f'samtools faidx data/ref/{CHR_GZ}'
    subprocess.run(cmd3.split(' '), stderr=subprocess.DEVNULL)
    logger.info('Done with index creation')

    return None

def main():
    # Get reference FASTA
    get_fasta()

    return None

if __name__ == '__main__':
    main()
