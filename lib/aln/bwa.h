/* bwa alignment base
 *
 * Newly added copyright in 2022
 * Copyright (c) 2022-2025 Jacob.Morrison@vai.org
 *
 * The MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef BWA_H_
#define BWA_H_

#include <stdint.h>
#include "bntseq.h"
#include "bwt.h"

#define BWA_IDX_BWT 0x1
#define BWA_IDX_BNS 0x2
#define BWA_IDX_PAC 0x4
#define BWA_IDX_ALL 0x7

#define BWA_CTL_SIZE 0x10000

typedef struct {
    bwt_t    bwt[2]; // FM-index
    bntseq_t *bns; // information on the reference sequences
    uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base

    int    is_shm;
    int64_t l_mem;
    uint8_t  *mem;
} bwaidx_t;

typedef struct {
    int l_seq, id;                /* check if l_seq can be unsigned? */
    char *name, *comment, *barcode, *umi, *qual, *sam; /* sam stored the end output of sam record */
    uint8_t *seq, *bisseq[2];
    uint8_t *seq0;              /* pointer to sequence beginning before clipping */
    int l_seq0;                 /* the original l_seq before clipping */
    int l_adaptor;              /* length of adaptor sequence from the 3' end */
    int clip5;                  /* actual length of sequence to clip from 5' */
    int clip3; /* actual length of sequence to clip from 3', should include l_adaptor */
} bseq1_t;

extern int bwa_verbose;
extern char bwa_rg_id[256];


#ifdef __cplusplus
extern "C" {
#endif

    void bseq1_code_nt4(bseq1_t *s);
    bseq1_t *bis_create_bseq1(char *seq1, char *seq2, int *n);

    bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);
    bseq1_t *bis_bseq_read(int chunk_size, uint8_t has_bc, int *n_, void *ks1_, void *ks2_);
    void bseq_classify(int n, bseq1_t *seqs, int m[2], bseq1_t *sep[2]);

    void bwa_fill_scmat(int a, int b, int8_t mat[25]);
    void bwa_fill_scmat_ct(int a, int b, int8_t mat[25]); /* WZBS */
    void bwa_fill_scmat_ga(int a, int b, int8_t mat[25]); /* WZBS */

    uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);
    uint32_t *bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM);
    uint32_t *bis_bwa_gen_cigar2(const int8_t mat[25], int o_del, int e_del, int o_ins, int e_ins, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM, uint32_t *ZC, uint32_t *ZR, int *bss_u, uint8_t parent); /* WZBS */

    char *bwa_idx_infer_prefix(const char *hint);
    void bwa_idx_load_bwt(const char *hint, uint8_t parent, bwt_t *bwt);

    bwaidx_t *bwa_idx_load_from_shm(const char *hint);
    bwaidx_t *bwa_idx_load_from_disk(const char *hint, int which);
    bwaidx_t *bwa_idx_load(const char *hint, int which);
    void bwa_idx_destroy(bwaidx_t *idx);
    int bwa_idx2mem(bwaidx_t *idx);
    int bwa_mem2idx(int64_t l_mem, uint8_t *mem, bwaidx_t *idx);

    void bwa_print_sam_hdr(const bntseq_t *bns, const char *hdr_line);
    char *bwa_set_rg(const char *s);
    char *bwa_insert_header(const char *s, char *hdr);

#ifdef __cplusplus
}
#endif

#endif
