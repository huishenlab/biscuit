#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "wqueue.h"
#include "wvec.h"
#include "wzmisc.h"

#include "khashl.h"

#define N_THREADS_DEFAULT 1
#define STEP_SIZE_DEFAULT 100000

// Command line interface configuration
typedef struct {
    int step;
    int n_threads;
} covg_conf_t;

static void covg_conf_init(covg_conf_t *conf) {
    conf->step = STEP_SIZE_DEFAULT;
    conf->n_threads = N_THREADS_DEFAULT;
}

// Coverage table data
typedef struct {
    uint32_t covg; /* coverage */
    uint32_t base; /* number of bases with coverage */
} covg_data_t;

static inline covg_data_t init_covg_data(uint32_t covg, uint32_t base) {
    covg_data_t out = { .covg = covg, .base = base };

    return out;
}

static inline int compare_covg_data(const void *a, const void *b) {
    covg_data_t *a1 = (covg_data_t *)a;
    covg_data_t *b1 = (covg_data_t *)b;

    if (a1->covg < b1->covg) return -1;
    if (a1->covg > b1->covg) return  1;

    return 0;
}

DEFINE_VECTOR(covg_data_v, covg_data_t);

// Store two components of a fraction separately
// Allows for modifying components before calculating decimal value
typedef struct {
    uint32_t num; /* fraction numerator */
    uint32_t den; /* fraction denominator */
} fraction_t;

static inline fraction_t init_fraction() {
    fraction_t out;

    out.num = 0;
    out.den = 0;

    return out;
}

static inline double divide(fraction_t frac) {
    return (double)frac.num / (double)frac.den;
}

// Fraction components for calculating mean and standard deviation
// Keep code cleaner by storing all outputs together
typedef struct {
    fraction_t all_base;     /* all base coverage */
    fraction_t q40_base;     /* q40 base coverage */
    fraction_t all_cpg;      /* all cpg coverage */
    fraction_t q40_cpg;      /* q40 cpg coverage */
    fraction_t all_base_top; /* all base top gc content coverage */
    fraction_t q40_base_top; /* q40 base top gc content coverage */
    fraction_t all_cpg_top;  /* all cpg top gc content coverage */
    fraction_t q40_cpg_top;  /* q40 cpg top gc content coverage */
    fraction_t all_base_bot; /* all base bottom gc content coverage */
    fraction_t q40_base_bot; /* q40 base bottom gc content coverage */
    fraction_t all_cpg_bot;  /* all cpg bottom gc content coverage */
    fraction_t q40_cpg_bot;  /* q40 cpg bottom gc content coverage */
} total_coverage_t;

static inline total_coverage_t init_total_coverage() {
    total_coverage_t out;

    out.all_base     = init_fraction();
    out.q40_base     = init_fraction();
    out.all_cpg      = init_fraction();
    out.q40_cpg      = init_fraction();
    out.all_base_top = init_fraction();
    out.q40_base_top = init_fraction();
    out.all_cpg_top  = init_fraction();
    out.q40_cpg_top  = init_fraction();
    out.all_base_bot = init_fraction();
    out.q40_base_bot = init_fraction();
    out.all_cpg_bot  = init_fraction();
    out.q40_cpg_bot  = init_fraction();

    return out;
}

// Initialize hashmap for storing coverages and corresponding number of bases
KHASHL_MAP_INIT(KH_LOCAL, covg_map, cm, uint32_t, uint32_t, kh_hash_uint32, kh_eq_generic)

// Add key-value to bucket if it doesn't exist, otherwise increment value
static inline void increment_map(covg_map *map, khint_t bucket, int absent) {
    if (absent) {
        kh_val(map, bucket) = 1;
    } else {
        kh_val(map, bucket) += 1;
    }
}

// Keep code cleaner by storing hashmaps together for each output
typedef struct {
    covg_map *all_base;     /* all reads base coverage */
    covg_map *q40_base;     /* q40 reads cpg coverage */
    covg_map *all_cpg;      /* all reads base coverage */
    covg_map *q40_cpg;      /* q40 reads cpg coverage */
    covg_map *all_base_top; /* all reads base top GC content coverage */
    covg_map *q40_base_top; /* q40 reads base top GC content coverage */
    covg_map *all_cpg_top;  /* all reads cpg top GC content coverage */
    covg_map *q40_cpg_top;  /* q40 reads cpg top GC content coverage */
    covg_map *all_base_bot; /* all reads base bottom GC content coverage */
    covg_map *q40_base_bot; /* q40 reads base bottom GC content coverage */
    covg_map *all_cpg_bot;  /* all reads cpg bottom GC content coverage */
    covg_map *q40_cpg_bot;  /* q40 reads cpg bottom GC content coverage */
} maps_t;

static inline maps_t *init_maps() {
    maps_t *out = calloc(1, sizeof(maps_t));

    out->all_base     = cm_init();
    out->q40_base     = cm_init();
    out->all_cpg      = cm_init();
    out->q40_cpg      = cm_init();
    out->all_base_top = cm_init();
    out->q40_base_top = cm_init();
    out->all_cpg_top  = cm_init();
    out->q40_cpg_top  = cm_init();
    out->all_base_bot = cm_init();
    out->q40_base_bot = cm_init();
    out->all_cpg_bot  = cm_init();
    out->q40_cpg_bot  = cm_init();

    return out;
}

static inline void destroy_maps(maps_t *maps) {
    cm_destroy(maps->q40_cpg_bot);
    cm_destroy(maps->all_cpg_bot);
    cm_destroy(maps->q40_base_bot);
    cm_destroy(maps->all_base_bot);
    cm_destroy(maps->q40_cpg_top);
    cm_destroy(maps->all_cpg_top);
    cm_destroy(maps->q40_base_top);
    cm_destroy(maps->all_base_top);
    cm_destroy(maps->q40_cpg);
    cm_destroy(maps->all_cpg);
    cm_destroy(maps->q40_base);
    cm_destroy(maps->all_base);

    free(maps);
}

// Test and set bits for whether a region covers a certain base
#define region_test(regs, i) regs[(i)>>3]&(1<<((i)&0x7))
#define region_set(regs, i) regs[(i)>>3] |= 1<<((i)&0x7)

// Information stored for each window
typedef struct {
    int64_t  block_id; /* ID of block processed by thread */
    maps_t  *maps;     /* coverage hash maps */
} qc_cov_record_t;

// Queue for processing records as they come in from various threads
DEFINE_WQUEUE(qc_cov_record, qc_cov_record_t)

// Window blocks for processing regions
typedef struct {
    int64_t   block_id; /* ID number for window */
    int32_t   tid;      /* contig ID number of region */
    uint32_t  beg;      /* beginning of region for window */
    uint32_t  end;      /* end of region for window */
    uint8_t  *cpg;      /* bit array of CpG locations */
    uint8_t  *top;      /* bit array of top GC content regions */
    uint8_t  *bot;      /* bit array of bottom GC content regions */
} qc_cov_window_t;

// Queue for pulling windows for multithreaded processing
DEFINE_WQUEUE(qc_cov_window, qc_cov_window_t)

// Information shared across threads
typedef struct {
    char             *bam_fn;      /* BAM filename */
    wqueue_t(qc_cov_window) *q;           /* window queue */
    wqueue_t(qc_cov_record) *rq;          /* records queue */
} result_t;

// Information needed for writing outputs
typedef struct {
    uint8_t has_top;
    uint8_t has_bot;
    wqueue_t(qc_cov_record) *q;
    char *prefix;
} writer_conf_t;

// Collect output file names together
typedef struct {
    char *cv_table;     /* coefficient of variation table */
    char *all_base;     /* all base coverage */
    char *q40_base;     /* Q40 base coverage */
    char *all_cpg;      /* all cpg coverage */
    char *q40_cpg;      /* Q40 cpg coverage */
    char *all_base_top; /* all Top GC content base coverage */
    char *q40_base_top; /* Q40 Top GC content base coverage */
    char *all_cpg_top;  /* All Top GC content cpg coverage */
    char *q40_cpg_top;  /* Q40 Top GC content cpg coverage */
    char *all_base_bot; /* all Bottom GC content base coverage */
    char *q40_base_bot; /* Q40 Bottom GC content base coverage */
    char *all_cpg_bot;  /* All Bottom GC content cpg coverage */
    char *q40_cpg_bot;  /* Q40 Bottom GC content cpg coverage */
} output_names_t;

static inline output_names_t *init_output_names(char *prefix) {
    output_names_t *out = calloc(1, sizeof(output_names_t));

    // If we have a prefix to prepend, then add 1 to its length to account for
    // underscore (_) connecting prefix to default file name
    size_t len_prefix = (prefix == NULL) ? 0 : strlen(prefix)+1;

    out->cv_table     = calloc(len_prefix + 15, sizeof(char));
    out->all_base     = calloc(len_prefix + 30, sizeof(char));
    out->q40_base     = calloc(len_prefix + 30, sizeof(char));
    out->all_cpg      = calloc(len_prefix + 30, sizeof(char));
    out->q40_cpg      = calloc(len_prefix + 30, sizeof(char));
    out->all_base_top = calloc(len_prefix + 35, sizeof(char));
    out->q40_base_top = calloc(len_prefix + 35, sizeof(char));
    out->all_cpg_top  = calloc(len_prefix + 35, sizeof(char));
    out->q40_cpg_top  = calloc(len_prefix + 35, sizeof(char));
    out->all_base_bot = calloc(len_prefix + 35, sizeof(char));
    out->q40_base_bot = calloc(len_prefix + 35, sizeof(char));
    out->all_cpg_bot  = calloc(len_prefix + 35, sizeof(char));
    out->q40_cpg_bot  = calloc(len_prefix + 35, sizeof(char));

    if (prefix != NULL) {
        strcat(out->cv_table, prefix);
        strcat(out->all_base, prefix);
        strcat(out->q40_base, prefix);
        strcat(out->all_cpg, prefix);
        strcat(out->q40_cpg, prefix);
        strcat(out->all_base_top, prefix);
        strcat(out->q40_base_top, prefix);
        strcat(out->all_cpg_top, prefix);
        strcat(out->q40_cpg_top, prefix);
        strcat(out->all_base_bot, prefix);
        strcat(out->q40_base_bot, prefix);
        strcat(out->all_cpg_bot, prefix);
        strcat(out->q40_cpg_bot, prefix);

        strcat(out->cv_table, "_");
        strcat(out->all_base, "_");
        strcat(out->q40_base, "_");
        strcat(out->all_cpg, "_");
        strcat(out->q40_cpg, "_");
        strcat(out->all_base_top, "_");
        strcat(out->q40_base_top, "_");
        strcat(out->all_cpg_top, "_");
        strcat(out->q40_cpg_top, "_");
        strcat(out->all_base_bot, "_");
        strcat(out->q40_base_bot, "_");
        strcat(out->all_cpg_bot, "_");
        strcat(out->q40_cpg_bot, "_");
    }

    strcat(out->cv_table, "cv_table.txt");
    strcat(out->all_base, "covdist_all_base_table.txt");
    strcat(out->q40_base, "covdist_q40_base_table.txt");
    strcat(out->all_cpg, "covdist_all_cpg_table.txt");
    strcat(out->q40_cpg, "covdist_q40_cpg_table.txt");
    strcat(out->all_base_top, "covdist_all_base_topgc_table.txt");
    strcat(out->q40_base_top, "covdist_q40_base_topgc_table.txt");
    strcat(out->all_cpg_top, "covdist_all_cpg_topgc_table.txt");
    strcat(out->q40_cpg_top, "covdist_q40_cpg_topgc_table.txt");
    strcat(out->all_base_bot, "covdist_all_base_botgc_table.txt");
    strcat(out->q40_base_bot, "covdist_q40_base_botgc_table.txt");
    strcat(out->all_cpg_bot, "covdist_all_cpg_botgc_table.txt");
    strcat(out->q40_cpg_bot, "covdist_q40_cpg_botgc_table.txt");

    return out;
}

static inline void destroy_output_names(output_names_t *get_wrecked) {
    free(get_wrecked->q40_cpg_bot);
    free(get_wrecked->all_cpg_bot);
    free(get_wrecked->q40_base_bot);
    free(get_wrecked->all_base_bot);
    free(get_wrecked->q40_cpg_top);
    free(get_wrecked->all_cpg_top);
    free(get_wrecked->q40_base_top);
    free(get_wrecked->all_base_top);
    free(get_wrecked->q40_cpg);
    free(get_wrecked->all_cpg);
    free(get_wrecked->q40_base);
    free(get_wrecked->all_base);
    free(get_wrecked->cv_table);

    free(get_wrecked);
}
