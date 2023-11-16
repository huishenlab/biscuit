/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2023 Jacob.Morrison@vai.org
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
**/

#ifndef _PILEUP_H_
#define _PILEUP_H_

#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <inttypes.h>
#include <libgen.h>
#include "wqueue.h"
#include "wzmisc.h"
#include "encode.h"
#include "sam.h"
#include "hts.h"
#include "refcache.h"
#include "kstring.h"
#include "wvec.h"
#include "stats.h"
#include "biscuit.h"

// Common parameters used in biscuit subcommands
typedef struct {
    uint8_t is_nome:1; /* input data is NOMe-seq */
    uint8_t verbose:1; /* print out extra information during processing */
} bisc_common_t;

static inline bisc_common_t bisc_common_init() {
    bisc_common_t out;

    out.is_nome = 0;
    out.verbose = 0;

    return out;
}

// Multithreading specific parameters
// TODO: can probably merge this into bisc_common_t
typedef struct {
    int step;      /* step size of window dispatching */
    int n_threads; /* number of processing threads */
} bisc_threads_t;

static inline bisc_threads_t bisc_threads_init() {
    bisc_threads_t out;

    out.step = 100000;
    out.n_threads = 3;

    return out;
}

// parameters affecting methylation extraction
typedef struct {
    // minimum values
    uint32_t min_base_qual;      /* minimum base quality */
    uint32_t min_read_len;       /* minimum read length */
    uint8_t  min_dist_end_5p;    /* minimum distance to 5' end of read to be included in methylation count */
    uint8_t  min_dist_end_3p;    /* minimum distance to 3' end of read to be included in methylation count */
    uint8_t  min_mapq;           /* minimum MAPQ score to include read */
    int      min_score;          /* minimum alignment score (from AS tag) */

    // maximum values
    int      max_nm;             /* maximum number of non-cytosine-conversion mismatches (from NM tag) */
    uint32_t max_retention;      /* maximum cytosine retention in read */

    // filter flags
    uint8_t  filter_ppair:1;     /* only use BAM_FPROPER_PAIR reads (hidden CLI option) */
    uint8_t  filter_secondary:1; /* don't include BAM_FSECONDARY reads (hidden CLI option) */
    uint8_t  filter_duplicate:1; /* don't include BAM_FDUP reads */
    uint8_t  filter_qcfail:1;    /* don't include BAM_FQCFAIL reads */
    uint8_t  filter_doublecnt:1; /* stop double-counting cytosine in overlapping mate reads */
} meth_filter_t;

static inline meth_filter_t meth_filter_init() {
    meth_filter_t out;

    out.min_base_qual = 20;
    out.min_read_len = 10;
    out.min_dist_end_5p = 3;
    out.min_dist_end_3p = 3;
    out.min_mapq = 40;
    out.min_score = 40;
    out.max_nm = 999999;
    out.max_retention = 999999;
    out.filter_ppair = 1;
    out.filter_secondary = 1;
    out.filter_duplicate = 1;
    out.filter_qcfail = 1;
    out.filter_doublecnt = 1;

    return out;
}

typedef struct {
    bisc_common_t comm; /* common parameters across subcommands */
    bisc_threads_t bt;  /* multithreading parameters */
    meth_filter_t filt; /* methylation extraction filters */

    // pileup only
    int somatic;                    /* call somatic mutation by assuming sample 1 is tumor and sample 2 is normal */
    uint8_t noheader:1;
    uint8_t ambi_redist:1;
    double error;
    double mu;
    double mu_somatic;
    double contam;
    double prior0;
    double prior1;
    double prior2;

    // epiread only
    int epiread_old;         /* print old BISCUIT epiread format */
    int print_all_locations; /* print all CpG and SNP locations in location column of epiread format */
    uint8_t is_long_read;           /* data is from long read sequencing */
    int epiread_pair;               /* pair output mode in epireads, doesn't mean "paired-end" */
    uint32_t epiread_reg_start;     /* first location of region provided to epiread */
    uint32_t epiread_reg_end;       /* final location of region provided to epiread */
    uint8_t filter_empty_epiread:1; /* remove epireads that only have F/x/P */
    uint32_t n_g_first;             /* number of OB/CTOB reads with a G in a CG context as first base, only incremented if not already filtered */
} conf_t;

void pileup_conf_init(conf_t *conf);

typedef struct {
    int64_t block_id;
    int32_t tid;
    uint32_t beg, end;
} window_t;

DEFINE_WQUEUE(window, window_t)

/* mutation-methylation code */
extern const char nt256int8_to_methcode[3];
extern const char nt256int8_to_basecode[7];
#define NSTATUS_METH 3
#define NSTATUS_BASE 7
typedef enum {METH_RETENTION, METH_CONVERSION, METH_NA} status_meth_t;
typedef enum {BASE_A, BASE_C, BASE_G, BASE_T, BASE_N, BASE_Y, BASE_R} status_base_t;

/* cytosine context code */
#define NCONTXTS 6 /* not including NA */
typedef enum {CTXT_HCG, CTXT_HCHG, CTXT_HCHH, CTXT_GCG, CTXT_GCHG, CTXT_GCHH, CTXT_NA} cytosine_context_t;
extern const char *cytosine_context[];

typedef struct {
    uint8_t sid;        /* which sample */
    uint8_t bsstrand:1; /* bisulfite strand */
    uint8_t qual:7;     /* base quality */
    uint8_t strand:1;   /* read stand */
    uint16_t qpos;      /* position on read */
    uint8_t cnt_ret;    /* count of retention of entire read */
    uint16_t rlen;      /* read length */
    char qb;            /* query base */
    uint8_t stat;       /* (status_base_t << 4) | status_meth_t */
} __attribute__((__packed__)) pileup_data_t;

DEFINE_VECTOR(pileup_data_v, pileup_data_t)

#define RECORD_QUEUE_END -2
#define RECORD_SLOT_OBSOLETE -1

typedef struct {
    int64_t block_id;
    kstring_t s; /* vcf record */

    /* coverage */
    int tid;

    double *betasum_context; /* CG, CHG, CHH */
    int64_t *cnt_context;    /* number of C's in each context */
} record_t;

DEFINE_VECTOR(record_v, record_t)

DEFINE_WQUEUE(record, record_t)

void pop_record_by_block_id(record_v *records, int64_t block_id, record_t *record);
void put_into_record_v(record_v *records, record_t rec);

uint32_t cnt_retention(refcache_t *rs, bam1_t *b, uint8_t bsstrand);

uint8_t infer_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual);

uint8_t get_bsstrand(refcache_t *rs, bam1_t *b, uint32_t min_base_qual, int allow_u);

uint32_t get_mate_length(char *m_cigar);

typedef struct {
    int32_t tid;
    char *name;
    uint32_t len;
} target_t;

DEFINE_VECTOR(target_v, target_t);

static inline int compare_targets(const void *a, const void *b) {
    return strcmp(((target_t*)a)->name, ((target_t*)b)->name);
}

cytosine_context_t fivenuc_context(refcache_t *rs, uint32_t rpos, char rb, char *fivenuc);

#define max(a,b)                \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b);   \
     _a > _b ? _a : _b; })

#define min(a,b)                \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b);   \
     _a > _b ? _b : _a; })

void pileup_genotype(int cref, int altsupp, conf_t *conf, char gt[4], double *_gl0, double *_gl1, double *_gl2, double *_gq);
int reference_supp(int cnts[9]);

static inline int compare_supp(const void *a, const void *b) {
    return ((*(uint32_t*)b)>>4) - ((*(uint32_t*)a)>>4);
}

typedef struct {
    wqueue_t(record) *q;
    int n_bams;
    char **bam_fns;
    char *outfn;
    char *statsfn;
    char *header;
    target_v *targets;
    conf_t *conf;
} writer_conf_t;

void *write_func(void *data);

static inline void pileup_parse_region(const char *reg, void *hdr, int *tid, int *beg, int *end) {
    const char *q = hts_parse_reg(reg, beg, end);
    if (q) {
        char *tmp = (char*)malloc(q - reg + 1);
        strncpy(tmp, reg, q - reg);
        tmp[q - reg] = 0;
        *tid = bam_name2id(hdr, tmp);
        free(tmp);
    }
    else {
        // not parsable as a region, but possibly a sequence named "foo:a"
        *tid = bam_name2id(hdr, reg);
        *beg = 0; *end = INT_MAX;
    }
}

#ifndef kroundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#endif /* _PILEUP_H_ */
