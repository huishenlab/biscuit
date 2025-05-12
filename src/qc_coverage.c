/* find coverage metrics from a bam file
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2024 Jacob.Morrison@vai.org
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
 *
 */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "zlib.h"

#include "wqueue.h"
#include "wvec.h"
#include "wzmisc.h"

#include "hts.h"
#include "sam.h"
#include "tbx.h"

#include "khashl.h"

#include "qc_coverage.h"
#include "bisc_utils.h"

// Merge contents of one hashmap into another
// Tabulate mean coverage fraction components at the same time
static void merge(covg_map *merge_from, covg_map *merge_into, fraction_t *frac) {
    khint_t k;
    int absent;
    kh_foreach(merge_from, k) {
        uint32_t covg = kh_key(merge_from, k); // coverage level
        uint32_t base = kh_val(merge_from, k); // number of bases with coverage level

        // Tabulating weighted mean fraction components
        frac->num += covg * base;
        frac->den += base;

        khint_t key = cm_put(merge_into, covg, &absent);
        if (absent) {
            kh_val(merge_into, key) = base;
        } else {
            kh_val(merge_into, key) += base;
        }
    }
}

// Calculate mean, standard deviation, and coefficient of variation
// Write coverage values and number of bases to output files at the same time
static void process_coverage_results(covg_map *cm, fraction_t frac, char *covg_fname, char *covdist_tag, FILE *cv_file, char *cv_tag) {
    // Only calculate statistics if we have non-zero fraction components
    uint8_t is_nonzero = frac.num > 0 && frac.den > 0;

    double mean = -1.0;
    if (is_nonzero) {
        mean = divide(frac);
    }

    FILE *out = fopen(covg_fname, "w");
    fprintf(out, "BISCUITqc Depth Distribution - %s\ndepth\tcount\n", covdist_tag);

    khint_t k;
    uint32_t variance_numerator = 0;
    covg_data_v *table = init_covg_data_v(50);
    kh_foreach(cm, k) {
        uint32_t covg = kh_key(cm, k);
        uint32_t base = kh_val(cm, k);

        if (is_nonzero) {
            variance_numerator += base * (covg - mean) * (covg - mean);
        }

        push_covg_data_v(table, init_covg_data(covg, base));
    }

    qsort(table->buffer, table->size, sizeof(covg_data_t), compare_covg_data);

    size_t i;
    for (i=0; i<table->size; ++i) {
        covg_data_t *cd = ref_covg_data_v(table, i);
        fprintf(out, "%u\t%u\n", cd->covg, cd->base);
    }

    free_covg_data_v(table);
    fflush(out);
    fclose(out);

    if (is_nonzero) {
        double sigma = sqrt((double)variance_numerator / (double)frac.den);
        fprintf(cv_file, "%s\t%lf\t%lf\t%lf\n", cv_tag, mean, sigma, sigma/mean);
    }
}

// Collect results from threads and write to various output files
static void *coverage_write_func(void *data) {
    writer_conf_t *c = (writer_conf_t*) data;

    maps_t *maps = init_maps();
    filter_counts_t *counts = init_filter_counts();
    total_coverage_t covg_fracs = init_total_coverage();

    // Effectively a reduction algorithm on the hashmaps from the individual records
    while (1) {
        qc_cov_record_t rec;
        wqueue_get(qc_cov_record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;

        add_filter_counts(counts, rec.counts);

        merge(rec.maps->all_base, maps->all_base, &covg_fracs.all_base);
        merge(rec.maps->q40_base, maps->q40_base, &covg_fracs.q40_base);
        merge(rec.maps->all_cpg , maps->all_cpg , &covg_fracs.all_cpg );
        merge(rec.maps->q40_cpg , maps->q40_cpg , &covg_fracs.q40_cpg );
        if (c->has_top) {
            merge(rec.maps->all_base_top, maps->all_base_top, &covg_fracs.all_base_top);
            merge(rec.maps->q40_base_top, maps->q40_base_top, &covg_fracs.q40_base_top);
            merge(rec.maps->all_cpg_top, maps->all_cpg_top, &covg_fracs.all_cpg_top);
            merge(rec.maps->q40_cpg_top, maps->q40_cpg_top, &covg_fracs.q40_cpg_top);
        }
        if (c->has_bot) {
            merge(rec.maps->all_base_bot, maps->all_base_bot, &covg_fracs.all_base_bot);
            merge(rec.maps->q40_base_bot, maps->q40_base_bot, &covg_fracs.q40_base_bot);
            merge(rec.maps->all_cpg_bot, maps->all_cpg_bot, &covg_fracs.all_cpg_bot);
            merge(rec.maps->q40_cpg_bot, maps->q40_cpg_bot, &covg_fracs.q40_cpg_bot);
        }

        destroy_filter_counts(rec.counts);
        destroy_maps(rec.maps);
    }

    // Write collected results to output files
    output_names_t *names = init_output_names(c->prefix);

    FILE *cv_table = fopen(names->cv_table, "w");
    fprintf(cv_table, "BISCUITqc Uniformity Table\ngroup\tmu\tsigma\tcv\n");

    process_coverage_results(maps->all_base, covg_fracs.all_base, names->all_base, "All Bases", cv_table, "all_base");
    process_coverage_results(maps->q40_base, covg_fracs.q40_base, names->q40_base, "Q40 Bases", cv_table, "q40_base");
    process_coverage_results(maps->all_cpg , covg_fracs.all_cpg , names->all_cpg , "All CpGs" , cv_table, "all_cpg" );
    process_coverage_results(maps->q40_cpg , covg_fracs.q40_cpg , names->q40_cpg , "Q40 CpGs" , cv_table, "q40_cpg" );
    if (c->has_top) {
        process_coverage_results(maps->all_base_top, covg_fracs.all_base_top, names->all_base_top, "All Top GC Bases", cv_table, "all_base_topgc");
        process_coverage_results(maps->q40_base_top, covg_fracs.q40_base_top, names->q40_base_top, "Q40 Top GC Bases", cv_table, "q40_base_topgc");
        process_coverage_results(maps->all_cpg_top, covg_fracs.all_cpg_top, names->all_cpg_top, "All Top GC CpGs", cv_table, "all_cpg_topgc");
        process_coverage_results(maps->q40_cpg_top, covg_fracs.q40_cpg_top, names->q40_cpg_top, "Q40 Top GC CpGs", cv_table, "q40_cpg_topgc");
    }
    if (c->has_bot) {
        process_coverage_results(maps->all_base_bot, covg_fracs.all_base_bot, names->all_base_bot, "All Bot GC Bases", cv_table, "all_base_botgc");
        process_coverage_results(maps->q40_base_bot, covg_fracs.q40_base_bot, names->q40_base_bot, "Q40 Bot GC Bases", cv_table, "q40_base_botgc");
        process_coverage_results(maps->all_cpg_bot, covg_fracs.all_cpg_bot, names->all_cpg_bot, "All Bot GC CpGs", cv_table, "all_cpg_botgc");
        process_coverage_results(maps->q40_cpg_bot, covg_fracs.q40_cpg_bot, names->q40_cpg_bot, "Q40 Bot GC CpGs", cv_table, "q40_cpg_botgc");
    }

    fflush(cv_table);
    fclose(cv_table);

    FILE *fh_counts = fopen(names->counts, "w");
    fprintf(fh_counts, "BISCUITqc Filtered Read Counts\ncategory\tcount\n");
    fprintf(fh_counts, "n_reads_in_bam\t%u\n", counts->n_reads);
    fprintf(fh_counts, "n_reads_passing_filters_all\t%u\n", counts->n_passed);
    fprintf(fh_counts, "n_reads_passing_filters_q40\t%u\n", counts->n_q40);
    fprintf(fh_counts, "n_unmapped_reads\t%u\n", counts->n_unmapped);
    fprintf(fh_counts, "n_secondary_reads\t%u\n", counts->n_secondary);
    fprintf(fh_counts, "n_duplicate_reads\t%u\n", counts->n_duplicate);
    fprintf(fh_counts, "n_improper_pair_reads\t%u\n", counts->n_improper_pair);
    fprintf(fh_counts, "n_qc_fail_reads\t%u\n", counts->n_qc_fail);
    fprintf(fh_counts, "n_too_short_reads\t%u\n", counts->n_too_short);
    fprintf(fh_counts, "n_too_many_mismatches\t%u\n", counts->n_nm_fail);
    fprintf(fh_counts, "n_alignment_score_too_low\t%u\n", counts->n_bad_as);
    fprintf(fh_counts, "n_too_many_retained_cytosines\t%u\n", counts->n_retention);

    fflush(fh_counts);
    fclose(fh_counts);

    destroy_output_names(names);
    destroy_filter_counts(counts);
    destroy_maps(maps);

    return 0;
}

// Take array of coverages and turn those coverages into a hashmap containing the coverage and number of bases
// (the number of elements with a given coverage in the array) with that coverage
static void format_coverage_data(maps_t *maps, uint32_t *all_covgs, uint32_t *q40_covgs, uint32_t arr_len, uint8_t *cpgs, uint8_t *tops, uint8_t *bots) {
    int abs_all_base, abs_q40_base, abs_all_cpg, abs_q40_cpg;
    int abs_all_base_top, abs_q40_base_top, abs_all_cpg_top, abs_q40_cpg_top;
    int abs_all_base_bot, abs_q40_base_bot, abs_all_cpg_bot, abs_q40_cpg_bot;
    khint_t k_all_base, k_q40_base, k_all_cpg, k_q40_cpg;
    khint_t k_all_base_top, k_q40_base_top, k_all_cpg_top, k_q40_cpg_top;
    khint_t k_all_base_bot, k_q40_base_bot, k_all_cpg_bot, k_q40_cpg_bot;

    uint32_t i;
    for (i=0; i<arr_len; i++) {
        uint8_t is_match_all = 0;
        uint8_t is_match_q40 = 0;

        // If the current coverage is the same as the last one, we know the coverage has been seen before
        // and the correct bucket is already in memory, so shortcircuit by automatically incrementing value
        if (i > 0) {
            is_match_all = all_covgs[i] == all_covgs[i-1];
            is_match_q40 = q40_covgs[i] == q40_covgs[i-1];

            if (is_match_all) {
                kh_val(maps->all_base, k_all_base) += 1;
            }
            if (is_match_q40) {
                kh_val(maps->q40_base, k_q40_base) += 1;
            }
        }

        if (!is_match_all) {
            k_all_base = cm_put(maps->all_base, all_covgs[i], &abs_all_base);
            increment_map(maps->all_base, k_all_base, abs_all_base);
        }

        if (!is_match_q40) {
            k_q40_base = cm_put(maps->q40_base, q40_covgs[i], &abs_q40_base);
            increment_map(maps->q40_base, k_q40_base, abs_q40_base);
        }

        // Always check if we're in CpG location. We'll also always pull the bucket for each hash map
        // to simplify finding buckets (rather than checking the last coverage like the all_base and
        // q40_base above)
        uint8_t is_cpg = region_test(cpgs, i);
        if (is_cpg) {
            k_all_cpg = cm_put(maps->all_cpg, all_covgs[i], &abs_all_cpg);
            increment_map(maps->all_cpg, k_all_cpg, abs_all_cpg);

            k_q40_cpg = cm_put(maps->q40_cpg, q40_covgs[i], &abs_q40_cpg);
            increment_map(maps->q40_cpg, k_q40_cpg, abs_q40_cpg);
        }

        // Only check the top GC content if the file is provided
        //
        // Since these regions are larger than CpGs, I may be able to keep track of the previous bucket
        // and reuse that to cut down on look up time, this would require knowing when we transition
        // from not being in a region, to being in the region and then acting appropriately. I may be
        // able to set a flag when I find I'm in the region, and then turn it off when I find I'm no
        // longer in the region. If the flag is on, then check if the coverage is the same, otherwise
        // pull the bucket. This also goes for the bottom GC content regions as well.
        if (tops) {
            if (region_test(tops, i)) {
                k_all_base_top = cm_put(maps->all_base_top, all_covgs[i], &abs_all_base_top);
                increment_map(maps->all_base_top, k_all_base_top, abs_all_base_top);

                k_q40_base_top = cm_put(maps->q40_base_top, q40_covgs[i], &abs_q40_base_top);
                increment_map(maps->q40_base_top, k_q40_base_top, abs_q40_base_top);

                if (is_cpg) {
                    k_all_cpg_top = cm_put(maps->all_cpg_top, all_covgs[i], &abs_all_cpg_top);
                    increment_map(maps->all_cpg_top, k_all_cpg_top, abs_all_cpg_top);

                    k_q40_cpg_top = cm_put(maps->q40_cpg_top, q40_covgs[i], &abs_q40_cpg_top);
                    increment_map(maps->q40_cpg_top, k_q40_cpg_top, abs_q40_cpg_top);
                }
            }
        }

        // Only check the bottom GC content if the file is provided
        if (bots) {
            if (region_test(bots, i)) {
                k_all_base_bot = cm_put(maps->all_base_bot, all_covgs[i], &abs_all_base_bot);
                increment_map(maps->all_base_bot, k_all_base_bot, abs_all_base_bot);

                k_q40_base_bot = cm_put(maps->q40_base_bot, q40_covgs[i], &abs_q40_base_bot);
                increment_map(maps->q40_base_bot, k_q40_base_bot, abs_q40_base_bot);

                if (is_cpg) {
                    k_all_cpg_bot = cm_put(maps->all_cpg_bot, all_covgs[i], &abs_all_cpg_bot);
                    increment_map(maps->all_cpg_bot, k_all_cpg_bot, abs_all_cpg_bot);

                    k_q40_cpg_bot = cm_put(maps->q40_cpg_bot, q40_covgs[i], &abs_q40_cpg_bot);
                    increment_map(maps->q40_cpg_bot, k_q40_cpg_bot, abs_q40_cpg_bot);
                }
            }
        }
    }
}

static uint8_t *set_bit_array(htsFile *bed, tbx_t *tbx, int tid, uint32_t window_beg, uint32_t window_end, int n_cols) {
    // Must be free'd when done being used
    uint8_t *out = calloc((window_end - window_beg)/8 + 1, sizeof(uint8_t));
    if (bed == NULL) {
        return out;
    }

    hts_itr_t *iter = tbx_itr_queryi(tbx, tid, window_beg, window_end);
    kstring_t str = {0, 256, (char *)calloc(256, sizeof(char))};

    char *tok = NULL;
    while (tbx_itr_next(bed, tbx, iter, &str) >= 0) {
        if (strcount_char(str.s, '\t') == n_cols) {
            // Chromosome
            tok = strtok(str.s, "\t");

            // Start
            tok = strtok(NULL, "\t");
            ensure_number(tok);
            uint32_t start = (uint32_t)atoi(tok);

            // End
            tok = strtok(NULL, "\t");
            ensure_number(tok);
            uint32_t end = (uint32_t)atoi(tok);

            uint32_t width = end - start;
            uint32_t i;
            for (i=0; i<width; i++) {
                if (start+i >= window_beg && start+i < window_end) {
                    region_set(out, start + i - window_beg);
                }
            }
        }
    }

    free(str.s);
    hts_itr_destroy(iter);

    return out;
}

// Function that is run by each thread, collects coverages across genome (factoring in CIGAR string), and then
// turns the results into a hashmap for downstream processing
static void *process_func(void *data) {
    result_t *res = (result_t*) data;
    covg_conf_t *conf = (covg_conf_t *) res->conf;

    // Open input file and check for existing BAM index
    htsFile   *in  = hts_open(res->bam_fn, "rb");
    hts_idx_t *idx = sam_index_load(in, res->bam_fn);
    if (!idx) {
        fprintf(stderr, "[%s:%d] BAM %s is not indexed?\n", __func__, __LINE__, res->bam_fn);
        fflush(stderr);
        exit(1);
    }
    bam_hdr_t *header = sam_hdr_read(in);

    uint32_t flank = 1000;
    refcache_t *rs = init_refcache(res->ref_fn, flank, flank);

    qc_cov_record_t rec;
    memset(&rec, 0, sizeof(qc_cov_record_t));
    qc_cov_window_t w;
    while (1) {
        wqueue_get(qc_cov_window, res->q, &w);
        if (w.tid == -1) break;

        // Tabulate coverages across window
        uint32_t *all_covgs = calloc(w.end - w.beg, sizeof(uint32_t));
        uint32_t *q40_covgs = calloc(w.end - w.beg, sizeof(uint32_t));

        // Prep read filter counts
        filter_counts_t *counts = init_filter_counts();

        // Iterate through reads that overlap window
        hts_itr_t *iter = sam_itr_queryi(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
        bam1_t *b = bam_init1();
        int ret;
        while ((ret = sam_itr_next(in, iter, b))>0) {
            bam1_core_t *c = &b->core;

            // Count all reads
            counts->n_reads++;

            // 0-based reference position
            uint32_t rpos = c->pos;

            // Read-based filtering
            if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->filt.min_read_len) {
                if (w.beg < rpos) counts->n_too_short++;
                continue;
            }
            if (c->flag > 0) { // only when any flag is set
                if (c->flag & BAM_FUNMAP) {
                    if (w.beg < rpos) counts->n_unmapped++;
                    continue;
                }
                if (conf->filt.filter_secondary && c->flag & BAM_FSECONDARY) {
                    if (w.beg < rpos) counts->n_secondary++;
                    continue;
                }
                if (conf->filt.filter_duplicate && c->flag & BAM_FDUP) {
                    if (w.beg < rpos) counts->n_duplicate++;
                    continue;
                }
                if (conf->filt.filter_ppair && c->flag & BAM_FPAIRED && !(c->flag & BAM_FPROPER_PAIR)) {
                    if (w.beg < rpos) counts->n_improper_pair++;
                    continue;
                }
                if (conf->filt.filter_qcfail && c->flag & BAM_FQCFAIL) {
                    if (w.beg < rpos) counts->n_qc_fail++;
                    continue;
                }
            }

            uint8_t *nm = bam_aux_get(b, "NM");
            if (nm && bam_aux2i(nm) > conf->filt.max_nm) {
                if (w.beg < rpos) counts->n_nm_fail++;
                continue;
            }

            uint8_t *as = bam_aux_get(b, "AS");
            if (as && bam_aux2i(as) < conf->filt.min_score) {
                if (w.beg < rpos) counts->n_bad_as++;
                continue;
            }

            // If the read is shorter than max_retention, then we can't ever have cnt_ret > max_retention
            // Therefore, only do the calculation if it might actually be true
            if (c->l_qseq < conf->filt.max_retention) {
                // Moved cache fetching into if-statement as this is expensive and only needs be done for
                // checking retention
                // If max_retention is low enough that this function is called for every read, we do have
                // extra calls being made to fetch, but refcache_fetch first checks to see if requested
                // region is already cached, so not too expensive to have it here
                char *chrm = header->target_name[w.tid];
                refcache_fetch(rs, chrm, w.beg>flank ? w.beg-flank : 1, w.end+flank);

                uint8_t bsstrand = get_bsstrand(rs, b, conf->filt.min_base_qual, 0);
                uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
                if (cnt_ret > conf->filt.max_retention) {
                    if (w.beg < rpos) counts->n_retention++;
                    continue;
                }
            }

            // Count reads passing all filters
            counts->n_passed++;
            if (c->qual >= 40) {
                if (w.beg < rpos) counts->n_q40++;
            }

            // Process CIGAR string to find if a base is covered or not
            int i;
            uint32_t j;
            for (i=0; i<c->n_cigar; ++i) {
                uint32_t op    = bam_cigar_op(bam_get_cigar(b)[i]);
                uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                switch(op) {
                    case BAM_CMATCH:
                    case BAM_CEQUAL:
                    case BAM_CDIFF:
                        for (j=0; j<oplen; ++j) {
                            if (w.beg <= rpos+j && rpos+j < w.end) {
                                uint32_t idx = rpos + j - w.beg;
                                all_covgs[idx] += 1;
                                if (c->qual >= 40) {
                                    q40_covgs[idx] += 1;
                                }
                            }
                        }
                        rpos += oplen;
                        break;
                    case BAM_CINS:
                    case BAM_CSOFT_CLIP:
                        // Normally, the read position would be incremented here, but we don't actually need
                        // the read position, so leaving these separate as a historical thing rather than a
                        // practical thing
                        break;
                    case BAM_CDEL:
                    case BAM_CREF_SKIP:
                        rpos += oplen;
                        break;
                    case BAM_CHARD_CLIP:
                    case BAM_CPAD:
                        break;
                    default:
                        fprintf(stderr, "Unknown cigar %u\n", op);
                        abort();
                }
            }
        }

        // Produce coverage output
        rec.maps = init_maps();
        format_coverage_data(rec.maps, all_covgs, q40_covgs, w.end-w.beg, w.cpg, w.top, w.bot);

        // Set record block id and put output maps and read counts into output queue
        rec.block_id = w.block_id;
        rec.counts = counts;
        wqueue_put2(qc_cov_record, res->rq, rec);

        // Clean up
        bam_destroy1(b);
        hts_itr_destroy(iter);

        free(q40_covgs);
        free(all_covgs);
        free(w.bot);
        free(w.top);
        free(w.cpg);
    }

    // Final clean up
    free_refcache(rs);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    hts_close(in);

    return 0;
}

void covg_conf_init(covg_conf_t *conf) {
    conf->bt = bisc_threads_init();
    conf->filt = meth_filter_init();
}

// Print usage for tool
static int usage() {
    covg_conf_t conf;
    covg_conf_init(&conf);

    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: coverage [options] <ref.fa> <cpgs.bed.gz> <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -P STR    Prefix for output file names\n");
    fprintf(stderr, "    -B STR    Bottom 10 percent GC content windows BED file\n");
    fprintf(stderr, "    -T STR    Top 10 percent GC content windows BED file\n");
    fprintf(stderr, "    -s INT    Step size of windows [%d]\n", conf.bt.step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf.bt.n_threads);
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -b INT    Minimum base quality [%u]\n", conf.filt.min_base_qual);
    fprintf(stderr, "    -a INT    Minimum alignment score (from AS-tag) [%u]\n", conf.filt.min_score);
    fprintf(stderr, "    -t INT    Max cytosine retention in a read [%u]\n", conf.filt.max_retention);
    fprintf(stderr, "    -l INT    Minimum read length [%u]\n", conf.filt.min_read_len);
    fprintf(stderr, "    -u        NO filtering of duplicate reads\n");
    fprintf(stderr, "    -p        NO filtering of improper pair\n");
    fprintf(stderr, "    -n INT    Maximum NM tag [%d]\n", conf.filt.max_nm);
    fprintf(stderr, "    -h        Print usage\n");
    fprintf(stderr, "\n");

    return 1;
}

// Main function
int main_qc_coverage(int argc, char *argv[]) {
    // Command line interface option values
    char *prefix = 0;
    char *top_fn = 0;
    char *bot_fn = 0;

    covg_conf_t conf;
    covg_conf_init(&conf);

    // Process command line arguments
    int c;
    if (argc < 2) { return usage(); }
    while ((c=getopt(argc, argv, ":@:B:P:T:a:b:l:n:s:t:chpu")) >= 0) {
        switch (c) {
            case '@': ensure_number(optarg); conf.bt.n_threads = atoi(optarg); break;
            case 'B': bot_fn = optarg; break;
            case 'P': prefix = optarg; break;
            case 's': ensure_number(optarg); conf.bt.step = atoi(optarg); break;
            case 'T': top_fn = optarg; break;
            case 't': conf.filt.max_retention = atoi(optarg); break;
            case 'l': conf.filt.min_read_len = atoi(optarg); break;
            case 'n': conf.filt.max_nm = atoi(optarg); break;
            case 'b': conf.filt.min_base_qual = atoi(optarg); break;
            case 'a': conf.filt.min_score = atoi(optarg); break;
            case 'c': conf.filt.filter_secondary = 0; break;
            case 'u': conf.filt.filter_duplicate = 0; break;
            case 'p': conf.filt.filter_ppair = 0; break;
            case 'h': usage(); return 0;
            case ':': usage(); fprintf(stderr, "Option needs an argument: -%c\n", optopt); return 1;
            case '?': usage(); fprintf(stderr, "Unrecognized option: -%c\n", optopt); return 1;
            default: return usage();
        }
    }

    // Missing required arguments
    if (optind + 3 > argc) {
        usage();
        fprintf(stderr, "CpG BED file, reference FASTA, or BAM input is missing\n");
        return 1;
    }

    // Set required argument values
    char *reffn = argv[optind++];
    char *cpg_bed_fn = argv[optind++];
    char *infn = argv[optind++];

    // Validate the number of threads and step size
    if (conf.bt.n_threads < 1) {
        fprintf(stderr, "[WARNING] Number of threads is < 1. Setting to default\n");
        conf.bt.n_threads = N_THREADS_DEFAULT;
    }
    if (conf.bt.step < 1) {
        fprintf(stderr, "[WARNING] Window step size is < 1. Setting to default\n");
        conf.bt.step = STEP_SIZE_DEFAULT;
    }

    // Read input BAM to get chromosome sizes for setting windows
    htsFile *in = hts_open(infn, "rb");
    if (in == NULL) {
        fprintf(stderr, "%s unable to be opened\n", infn);
        return 1;
    }
    bam_hdr_t *header = sam_hdr_read(in);

    // Setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .has_top = (top_fn) ? 1 : 0,
        .has_bot = (bot_fn) ? 1 : 0,
        .q = wqueue_init(qc_cov_record, conf.bt.step),
        .prefix = prefix ? prefix : NULL,
    };
    pthread_create(&writer, NULL, coverage_write_func, &writer_conf);

    // Setup multithreading
    wqueue_t(qc_cov_window) *wq = wqueue_init(qc_cov_window, conf.bt.step);
    pthread_t *processors = calloc(conf.bt.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.bt.n_threads, sizeof(result_t));
    int i;
    for (i=0; i<conf.bt.n_threads; ++i) {
        results[i].conf = &conf;
        results[i].bam_fn = infn;
        results[i].ref_fn = reffn;
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    qc_cov_window_t w;
    memset(&w, 0, sizeof(qc_cov_window_t));
    uint32_t wbeg;
    int64_t block_id = 0;

    // Setup windows and add to queue for processing
    htsFile *cpg_bed = hts_open(cpg_bed_fn, "r");
    tbx_t *cpg_tbx = tbx_index_load3(cpg_bed_fn, NULL, 0);
    if (cpg_tbx == NULL) {
        fprintf(stderr, "Make sure %s is tab-indexed via tabix\n", cpg_bed_fn);
        return 1;
    }

    htsFile *top_bed = top_fn ? hts_open(top_fn, "r") : NULL;
    tbx_t *top_tbx = top_fn ? tbx_index_load3(top_fn, NULL, 0) : NULL;
    if (top_bed != NULL && top_tbx == NULL) {
        fprintf(stderr, "Make sure %s is tab-indexed via tabix\n", top_fn);
        return 1;
    }

    htsFile *bot_bed = bot_fn ? hts_open(bot_fn, "r") : NULL;
    tbx_t *bot_tbx = bot_fn ? tbx_index_load3(bot_fn, NULL, 0) : NULL;
    if (bot_bed != NULL && bot_tbx == NULL) {
        fprintf(stderr, "Make sure %s is tab-indexed via tabix\n", bot_fn);
        return 1;
    }

    size_t j;
    for (j=0; j<header->n_targets; ++j) {
        uint32_t len = header->target_len[j];
        for (wbeg=0; wbeg<len; wbeg += conf.bt.step, block_id++) {
            w.tid = j;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.bt.step;
            if (w.end > len) w.end = len;
            w.cpg = set_bit_array(cpg_bed, cpg_tbx, j, wbeg, wbeg+conf.bt.step > len ? len : wbeg+conf.bt.step, 2);
            w.top = set_bit_array(top_bed, top_tbx, j, wbeg, wbeg+conf.bt.step > len ? len : wbeg+conf.bt.step, 3);
            w.bot = set_bit_array(bot_bed, bot_tbx, j, wbeg, wbeg+conf.bt.step > len ? len : wbeg+conf.bt.step, 3);
            wqueue_put(qc_cov_window, wq, &w);
        }
    }

    // "Windows" for telling thread pool we've reached the end of our queue
    for (i=0; i<conf.bt.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(qc_cov_window, wq, &w);
    }

    // Collect our threads back into serial processing
    for (i=0; i<conf.bt.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    // Tell the writer that we've reached the end of processing and can stop collecting results
    qc_cov_record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(qc_cov_record, writer_conf.q, rec);
    pthread_join(writer, NULL);

    // Clean up
    if (bot_fn) {
        tbx_destroy(bot_tbx);
        hts_close(bot_bed);
    }
    if (top_fn) {
        tbx_destroy(top_tbx);
        hts_close(top_bed);
    }
    tbx_destroy(cpg_tbx);
    hts_close(cpg_bed);
    wqueue_destroy(qc_cov_record, writer_conf.q);
    free(results);
    free(processors);
    wqueue_destroy(qc_cov_window, wq);
    hts_close(in);
    bam_hdr_destroy(header);

    return 0;
}
