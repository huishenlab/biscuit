import itertools
import hashlib
import logging

logger = logging.getLogger(__name__)

def compare_files(ext, prefix_1, prefix_2):
    """Compare two output files."""
    # sha256 hash file 1 contents
    with open(f'{prefix_1}{ext}', 'rb') as f:
        digest_1 = hashlib.file_digest(f, 'sha256').hexdigest()

    # sha256 hash file 2 contents
    with open(f'{prefix_2}{ext}', 'rb') as f:
        digest_2 = hashlib.file_digest(f, 'sha256').hexdigest()

    logging.debug(f'{prefix_1}{ext} sha256 hash: {digest_1}')
    logging.debug(f'{prefix_2}{ext} sha256 hash: {digest_2}')

    return digest_1 == digest_2

def compare_line_by_line(ext, prefix_1, prefix_2):
    """Compare two output files line by line."""
    out = []
    with open(f'{prefix_1}{ext}', 'r') as f1, open(f'{prefix_2}{ext}', 'r') as f2:
        # n_lines counts the number of shared lines so we can start at the right line number for extra lines in one file
        # but not the other
        n_lines = 0
        for idx, (line1, line2) in enumerate(itertools.zip_longest(f1, f2, fillvalue='(empty)'), start=1):
            n_lines = idx
            l1 = line1.strip()
            l2 = line2.strip()

            # Collect lines that differ for printing downstream
            if l1 != l2:
                out.append( (idx, l1, l2) )

    return out
