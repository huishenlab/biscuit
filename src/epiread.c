/* convert bam to epiread format with supplied SNP bed file
 * 
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2021 Wanding.Zhou@vai.org
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
 **/
#include "pileup.h"
#include "wzmisc.h"

DEFINE_VECTOR(int_v, int)
DEFINE_VECTOR(char_v, char)

typedef struct episnp_chrom1_t {
    char *chrm;
    size_t n;
    uint32_t *locs;  // snp locations
} episnp_chrom1_t;

DEFINE_VECTOR(episnp_chrom1_v, episnp_chrom1_t)

void destroy_episnp(episnp_chrom1_v *episnp) {
    uint32_t i;
    for (i=0; i<episnp->size; ++i) {
        episnp_chrom1_t *e = ref_episnp_chrom1_v(episnp, i);
        free(e->locs);
    }
    free_episnp_chrom1_v(episnp);
}

// get all episnps from one chromosomes
static inline episnp_chrom1_t *get_episnp1(episnp_chrom1_v *episnp, char *chrm) {

    uint32_t i;
    episnp_chrom1_t *episnp1;
    for (i=0; i<episnp->size; ++i) {
        episnp1 = ref_episnp_chrom1_v(episnp, i);
        if (strcmp(episnp1->chrm, chrm) == 0) return episnp1;
    }
    return NULL;
}

static inline episnp_chrom1_t *get_n_insert_episnp1(episnp_chrom1_v *episnp, char *chrm) {

    episnp_chrom1_t *episnp1 = get_episnp1(episnp, chrm);
    if (!episnp1) {
        episnp1 = next_ref_episnp_chrom1_v(episnp);
        episnp1->chrm = strdup(chrm);
        episnp1->locs = NULL;
        episnp1->n = 0;
    }
    return episnp1;
}

typedef struct {
    char *bam_fn; // on stack
    char *ref_fn; // on stack
    episnp_chrom1_v *snp;
    wqueue_t(window) *q;
    wqueue_t(record) *rq;
    conf_t *conf;
} result_t;

static void *epiread_write_func(void *data) {

    writer_conf_t *c = (writer_conf_t*) data;

    FILE *out;
    if (c->outfn) out=fopen(c->outfn, "w");
    else out=stdout;

    int64_t next_block = 0;
    record_v *records = init_record_v(20);

    while (1) {
        record_t rec;
        wqueue_get(record, c->q, &rec);
        if(rec.block_id == RECORD_QUEUE_END) break;
        if (rec.block_id == next_block) {
            do {
                if (rec.s.s)
                    fputs(rec.s.s, out);
                free(rec.s.s);

                /* get next block from shelf if available else return OBSOLETE 
                 * and retrieve new block from queue */
                next_block++;
                pop_record_by_block_id(records, next_block, &rec);
            } while (rec.block_id != RECORD_SLOT_OBSOLETE);
        } else { // shelf the block if not next
            put_into_record_v(records, rec);
        }
    }

    free_record_v(records);
    if (c->outfn) { // for stdout, will close at the end of main
        fflush(out);
        fclose(out);
    }
    return 0;
}

#define episnp_test(snps, i) snps[(i)>>3]&(1<<((i)&0x7))
#define episnp_set(snps, i) snps[(i)>>3] |= 1<<((i)&0x7)

static void format_epiread(
        kstring_t *epi, bam1_t *b, uint8_t bsstrand, char *chrm, window_t *w, uint8_t *snps, conf_t *conf,
        int_v *snp_p, int_v *hcg_p, int_v *gch_p, int_v *cg_p,
        char_v *snp_c, char_v *hcg_c, char_v *gch_c, char_v *cg_c) {

    uint32_t k;

    if (conf->is_nome) { // nome-seq
        int first_epi = get_int_v(hcg_p, 0) < get_int_v(gch_p, 0) ? get_int_v(hcg_p, 0) : get_int_v(gch_p, 0);

        if (first_epi > 0 && (unsigned) first_epi >= w->beg && (unsigned) first_epi < w->end) {
            ksprintf(epi, "%s\t%s\t%c\t%c",
                    chrm,
                    bam_get_qname(b),
                    (b->core.flag&BAM_FREAD2) ? '2' : '1',
                    bsstrand ? '-' : '+');

            // HCG context (0-based)
            if (hcg_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(hcg_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<hcg_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(hcg_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(hcg_c, 0));
                for (k=1; k<hcg_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(hcg_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // GCH context (0-based)
            if (gch_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(gch_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<gch_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(gch_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(gch_c, 0));
                for (k=1; k<gch_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(gch_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // SNP (0-based)
            if (snp_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(snp_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<snp_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(snp_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(snp_c, 0));
                for (k=1; k<snp_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(snp_c, k));
            } else if (snps) {
                kputs("\t.\t.", epi);
            } else {
                kputs("\t\t", epi);
            }

            // end line
            kputc('\n', epi);
        }
    } else { // bs-seq
        if (get_int_v(cg_p, 0) > 0 && (unsigned) get_int_v(cg_p, 0) >= w->beg && (unsigned) get_int_v(cg_p, 0) < w->end) {
            ksprintf(epi, "%s\t%s\t%c\t%c",
                    chrm,
                    bam_get_qname(b),
                    (b->core.flag&BAM_FREAD2) ? '2' : '1',
                    bsstrand ? '-' : '+');

            // CpG context (0-based)
            if (cg_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(cg_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<cg_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(cg_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(cg_c, 0));
                for (k=1; k<cg_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(cg_c, k));
            } else {
                kputs("\t.\t.", epi);
            }

            // SNP (0-based)
            if (snp_p->size > 0) {
                ksprintf(epi, "\t%d", get_int_v(snp_p, 0)-1);
                if (conf->print_all_locations) {
                    for (k=1; k<snp_p->size; ++k)
                        ksprintf(epi, ",%d", get_int_v(snp_p, k)-1);
                }

                ksprintf(epi, "\t%c", get_char_v(snp_c, 0));
                for (k=1; k<snp_c->size; ++k)
                    ksprintf(epi, "%c", get_char_v(snp_c, k));
            } else if (snps) {
                kputs("\t.\t.", epi);
            } else {
                kputs("\t\t", epi);
            }

            // end line
            kputc('\n', epi);
        }
    }
}

// format one bam record to epi-read format
static void format_epiread_pairwise(
        kstring_t *epi, char *chrm, window_t *w, conf_t *conf,
        int_v *snp_p, int_v *hcg_p, int_v *gch_p, int_v *cg_p,
        char_v *snp_c, char_v *hcg_c, char_v *gch_c, char_v *cg_c) {

    uint32_t j, k;
    for (k=0; k<snp_p->size; ++k) {
        // avoid double counting between windows
        if (!((unsigned) get_int_v(snp_p, k) >= w->beg && (unsigned) get_int_v(snp_p, k) < w->end))
            continue;

        if (conf->is_nome) { // nome-seq
            for (j=0; j<hcg_p->size; ++j) { // SNP and HCG context
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(hcg_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(hcg_c, j));
            }
            for (j=0; j<gch_p->size; ++j) { // SNP and GCH context
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(gch_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(gch_c, j));
            }
        } else { // bs-seq
            for (j=0; j<cg_p->size; ++j) {
                // chrm, snp position, cpg position, snp calling, cytosine calling
                ksprintf(epi, "%s\t%d\t%d\t%c\t%c\n",
                        chrm,
                        get_int_v(snp_p, k),
                        get_int_v(cg_p, j),
                        get_char_v(snp_c, k),
                        get_char_v(cg_c, j));
            }
        }
    }
}


static void *process_func(void *data) {

    result_t *res = (result_t*) data;
    conf_t *conf = (conf_t*) res->conf;
    htsFile *in = hts_open(res->bam_fn, "rb");
    hts_idx_t *idx = sam_index_load(in, res->bam_fn);
    if (!idx) {
        fprintf(stderr, "[%s:%d] BAM %s is not indexed?\n",
                __func__, __LINE__, res->bam_fn);
        fflush(stderr);
        exit(1);
    }
   
    bam_hdr_t *header = sam_hdr_read(in);
    refcache_t *rs = init_refcache(res->ref_fn, 1000, 1000);
    uint32_t j;

    record_t rec;
    memset(&rec, 0, sizeof(record_t));
    window_t w;
    while (1) {
        wqueue_get(window, res->q, &w);
        if (w.tid == -1) break;

        rec.tid = w.tid;
        char *chrm = header->target_name[w.tid];

        uint32_t snp_beg = w.beg>1000 ? w.beg-1000 : 1; // start location of snps
        uint32_t snp_end = w.end+1000;
        // make snp lookup table
        uint8_t *snps = NULL;
        if (res->snp) { // if snp is supplied
            snps = calloc((snp_end-snp_beg)/8+1, sizeof(uint8_t));
            episnp_chrom1_t *episnp1 = get_episnp1(res->snp, chrm);
            if (episnp1) { // if chromosome is found in snp file
                for (j=0; j<episnp1->n; ++j) {
                    uint32_t l = episnp1->locs[j];
                    if (l>=snp_beg && l<snp_end) {
                        episnp_set(snps, l-snp_beg);
                    }
                }
            }
        }
    
        rec.s.l = rec.s.m = 0; rec.s.s = 0; // the epiread string

        refcache_fetch(rs, chrm, w.beg>100?w.beg-100:1, w.end+100);
        hts_itr_t *iter = sam_itr_queryi(idx, w.tid, w.beg>1?(w.beg-1):1, w.end);
        bam1_t *b = bam_init1();
        int ret;

        while ((ret = sam_itr_next(in, iter, b))>0) {
            uint8_t bsstrand = get_bsstrand(rs, b, conf->min_base_qual, 0);

            // read-based filtering
            bam1_core_t *c = &b->core;
            if (c->qual < conf->min_mapq) continue;
            if (c->l_qseq < 0 || (unsigned) c->l_qseq < conf->min_read_len) continue;
            if (c->flag > 0) { // only when any flag is set
                if (conf->filter_secondary && c->flag & BAM_FSECONDARY) continue;
                if (conf->filter_duplicate && c->flag & BAM_FDUP) continue;
                if (conf->filter_ppair && c->flag & BAM_FPAIRED && !(c->flag & BAM_FPROPER_PAIR)) continue;
                if (conf->filter_qcfail && c->flag & BAM_FQCFAIL) continue;
            }

            uint8_t *nm = bam_aux_get(b, "NM");
            if (nm && bam_aux2i(nm)>conf->max_nm) continue;

            uint8_t *as = bam_aux_get(b, "AS");
            if (as && bam_aux2i(as) < conf->min_score) continue;

            uint32_t cnt_ret = cnt_retention(rs, b, bsstrand);
            if (cnt_ret > conf->max_retention) continue;

            // pairwise epiread format variables
            int_v  *snp_p = init_int_v(10);  // snp position
            char_v *snp_c = init_char_v(10); // snp character
            int_v  *cg_p=0, *hcg_p=0, *gch_p=0;
            char_v *cg_c=0, *hcg_c=0, *gch_c=0;
            if (conf->is_nome) {
                hcg_p = init_int_v(10);  // hcg positions
                hcg_c = init_char_v(10); // hcg characters
                gch_p = init_int_v(10);  // gch positions
                gch_c = init_char_v(10); // gch characters
            } else {
                cg_p = init_int_v(10);  // cpg positions
                cg_c = init_char_v(10); // cpg characters
            }

            int i; uint32_t j;
            uint32_t rpos = c->pos+1, qpos = 0, rmpos = c->mpos + 1;

            char rb, qb;
            for (i=0; i<c->n_cigar; ++i) {
                uint32_t op    = bam_cigar_op(bam_get_cigar(b)[i]);
                uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
                switch(op) {
                    case BAM_CMATCH:
                        for (j=0; j<oplen; ++j) {
                            rb = refcache_getbase_upcase(rs, rpos+j);
                            qb = bscall(b, qpos+j);

                            // skip bases with low base quality
                            if (bam_get_qual(b)[qpos+j] < conf->min_base_qual)
                                continue;

                            // read-position-based filtering
                            if (qpos+j+1 < conf->min_dist_end_5p ||
                                c->l_qseq < (int32_t)(qpos+j+1 + conf->min_dist_end_3p))
                                continue;

                            /* If read 2 in a proper pair, skip counting overlapped cytosines
                             * Right now I assume read 1 and read 2 are the same length and there is no gap.
                             * Better solution would be to log mate end in the alignment (using bam_endpos).
                             * The filtering of double counting is only effective when reads are properly paired.
                             *  
                             * The filtering remove bases from read 2 (usually the synthesized read)
                             * falling into the overlapped region.
                             */
                            if (conf->filter_doublecnt &&
                                (c->flag & BAM_FPROPER_PAIR) &&
                                (c->flag & BAM_FREAD2) &&
                                rpos+j >= max(rpos, rmpos) &&
                                rpos+j <= min(rpos + c->l_qseq, rmpos + c->l_qseq))
                                continue;
                            

                            // reference is a G
                            if (bsstrand && rb == 'G' && rpos+j-1 >= rs->beg) {
                                if (conf->is_nome) { // nome-seq
                                    if (rpos+j+1 <= rs->end) { // prevent overflow
                                        char rb0 = refcache_getbase_upcase(rs, rpos+j-1); // previous base
                                        char rb1 = refcache_getbase_upcase(rs, rpos+j+1); // next base
                                        if (rb0 == 'C' && rb1 != 'C') { // HCG context
                                            // Note: measure G in CpG context, record location of C
                                            push_int_v(hcg_p, (int) rpos+j-1);
                                            if (qb == 'A') {
                                                push_char_v(hcg_c, 'T');
                                            } else if (qb == 'G') {
                                                push_char_v(hcg_c, 'C');
                                            } else {
                                                push_char_v(hcg_c, 'N');
                                            }
                                        } else if (rb0 != 'C' && rb1 == 'C') { // GCH context
                                            push_int_v(gch_p, (int) rpos+j);
                                            if (qb == 'A') {
                                                push_char_v(gch_c, 'T');
                                            } else if (qb == 'G') {
                                                push_char_v(gch_c, 'C');
                                            } else {
                                                push_char_v(gch_c, 'N');
                                            }
                                        }
                                    }
                                } else { // bs-seq
                                    char rb0 = refcache_getbase_upcase(rs, rpos+j-1); // previous base
                                    if (rb0 == 'C') { // CpG context
                                        // Note: measure G in CpG context, record location of C
                                        push_int_v(cg_p, (int) rpos+j-1);
                                        if (qb == 'A') {
                                            push_char_v(cg_c, 'T');
                                        } else if (qb == 'G') {
                                            push_char_v(cg_c, 'C');
                                        } else {
                                            push_char_v(cg_c, 'N');
                                        }
                                    }
                                }
                            }

                            // reference is a C
                            if (!bsstrand && rb == 'C' && rpos+j+1 <= rs->end) {
                                if (conf->is_nome) { // nome-seq
                                    if (rpos+j-1 >= rs->beg) { // to prevent underflow
                                        char rb0 = refcache_getbase_upcase(rs, rpos+j-1); // previous base
                                        char rb1 = refcache_getbase_upcase(rs, rpos+j+1); // next base
                                        if (rb0 != 'G' && rb1 == 'G') { // HCG context
                                            // measure C in CpG context
                                            push_int_v(hcg_p, (int) rpos+j);
                                            if (qb == 'T') {
                                                push_char_v(hcg_c, 'T');
                                            } else if (qb == 'C') {
                                                push_char_v(hcg_c, 'C');
                                            } else {
                                                push_char_v(hcg_c, 'N');
                                            }
                                        } else if (rb0 == 'G' && rb1 != 'G') { // GCH context
                                            push_int_v(gch_p, (int) rpos+j);
                                            if (qb == 'T') {
                                                push_char_v(gch_c, 'T');
                                            } else if (qb == 'C') {
                                                push_char_v(gch_c, 'C');
                                            } else {
                                                push_char_v(gch_c, 'N');
                                            }
                                        }
                                    }
                                } else { // bs-seq
                                    char rb1 = refcache_getbase_upcase(rs, rpos+j+1); // next base
                                    if (rb1 == 'G') { // CpG context
                                        push_int_v(cg_p, (int) rpos+j-1);
                                        if (qb == 'T') {
                                            push_char_v(cg_c, 'T');
                                        } else if (qb == 'C') {
                                            push_char_v(cg_c, 'C');
                                        } else {
                                            push_char_v(cg_c, 'N');
                                        }
                                    }
                                }
                            }

                            // append SNP info if present
                            uint32_t snp_ind = rpos+j-snp_beg;
                            if (snps && episnp_test(snps, snp_ind)) {
                                push_char_v(snp_c, qb);
                                push_int_v(snp_p, rpos+j);
                            }
                        }
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CINS:
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                        rpos += oplen;
                        break;
                    case BAM_CSOFT_CLIP:
                        qpos += oplen;
                        break;
                    default:
                        fprintf(stderr, "Unknown cigar, %u\n", op);
                        abort();
                }
            }

            // produce epiread output
            if (conf->epiread_pair) {
                format_epiread_pairwise(
                        &rec.s, chrm, &w, conf,
                        snp_p, hcg_p, gch_p, cg_p,
                        snp_c, hcg_c, gch_c, cg_c);
            } else {
                format_epiread(
                        &rec.s, b, bsstrand, chrm, &w, snps, conf,
                        snp_p, hcg_p, gch_p, cg_p,
                        snp_c, hcg_c, gch_c, cg_c);
            }

            // clean up
            free_int_v(snp_p); free_char_v(snp_c);
            if (conf->is_nome) {
                free_int_v(hcg_p); free_char_v(hcg_c);
                free_int_v(gch_p); free_char_v(gch_c);
            } else {
                free_int_v(cg_p); free_char_v(cg_c);
            }
        }

        // run through cytosines
        rec.block_id = w.block_id;
        // put output string to output queue
        wqueue_put2(record, res->rq, rec);

        bam_destroy1(b);
        hts_itr_destroy(iter);
        free(snps);
    }

    free_refcache(rs);
    hts_close(in);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);

    return 0;
}

episnp_chrom1_v *bed_init_episnp(char *snp_bed_fn) {

    episnp_chrom1_v *episnp = init_episnp_chrom1_v(2);
    kstring_t line;
    line.l = line.m = 0; line.s = 0;

    // read SNP bed file
    episnp_chrom1_t *episnp1 = 0;
    char *tok;
    FILE *fh = fopen(snp_bed_fn,"r");
    while (1) {
        int c=fgetc(fh);
        if (c=='\n' || c==EOF) {
            if (strcount_char(line.s, '\t')>=2) {
                tok = strtok(line.s, "\t");

                if (!episnp1 || strcmp(episnp1->chrm, tok) != 0)
                    episnp1 = get_n_insert_episnp1(episnp, tok);

                episnp1->locs = realloc(episnp1->locs, (episnp1->n+1)*sizeof(uint32_t));
                tok = strtok(NULL, "\t");
                ensure_number(tok);
                episnp1->locs[episnp1->n] = atoi(tok)+1;
                episnp1->n++;
            }

            line.l = 0;
            if (c==EOF) {
                break;
            }
        } else {
            kputc(c, &line);
        }
    }
    free(line.s);

    return episnp;
}

static int usage(conf_t *conf) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit epiread [options] <ref.fa> <in.bam>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -B STR    Bed input for SNP display in epiread output\n");
    fprintf(stderr, "    -g STR    Region (optional, will process the whole bam if not specified)\n");
    fprintf(stderr, "    -s STR    Step of window dispatching [%d]\n", conf->step);
    fprintf(stderr, "    -@ INT    Number of threads [%d]\n", conf->n_threads);
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o STR    Output file [stdout]\n");
    fprintf(stderr, "    -P        Pairwise mode [off]\n");
    fprintf(stderr, "    -N        NOMe-seq mode [off]\n");
    fprintf(stderr, "    -A        Print all CpG and SNP locations in location column [off]\n");
    fprintf(stderr, "    -v        Verbose (print additional info for diagnostics) [off]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -b INT    Minimum base quality [%u]\n", conf->min_base_qual);
    fprintf(stderr, "    -m INT    Minimum mapping quality [%u]\n", conf->min_mapq);
    fprintf(stderr, "    -a INT    Minimum alignment score (from AS-tag) [%u]\n", conf->min_score);
    fprintf(stderr, "    -t INT    Max cytosine retention in a read [%u]\n", conf->max_retention);
    fprintf(stderr, "    -l INT    Minimum read length [%u]\n", conf->min_read_len);
    fprintf(stderr, "    -5 INT    Minimum distance to 5' end of a read [%u]\n", conf->min_dist_end_5p);
    fprintf(stderr, "    -3 INT    Minimum distance to 3' end of a read [%u]\n", conf->min_dist_end_3p);
    fprintf(stderr, "    -c        NO filtering secondary mapping\n");
    fprintf(stderr, "    -d        Double count cytosines in overlapping mate reads (avoided\n");
    fprintf(stderr, "                  by default)\n");
    fprintf(stderr, "    -u        NO filtering of duplicate\n");
    fprintf(stderr, "    -p        NO filtering of improper pair\n");
    fprintf(stderr, "    -n INT    Maximum NM tag [%d]\n", conf->max_nm);
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "\n");

    return 1;
}

int main_epiread(int argc, char *argv[]) {

    int c;
    char *reg = 0;
    char *outfn = 0;
    char *statsfn = 0;
    char *snp_bed_fn = 0;

    conf_t conf;
    memset(&conf, 0, sizeof(conf_t));
    conf.step = 100000;
    conf.n_threads = 3;
    conf.min_base_qual = 20;
    conf.min_mapq = 40;
    conf.min_score = 40;
    conf.max_retention = 999999;
    conf.min_read_len = 10;
    conf.min_dist_end_5p = 3;
    conf.min_dist_end_3p = 3;
    conf.filter_qcfail = 1;
    conf.filter_secondary = 1;
    conf.filter_doublecnt = 1;
    conf.filter_duplicate = 1;
    conf.filter_ppair = 1;
    conf.max_nm = 999999;
    conf.print_all_locations = 0;
    conf.is_nome = 0;
    conf.verbose = 0;
    conf.epiread_pair = 0;

    if (argc<2) return usage(&conf);
    while ((c=getopt(argc, argv, ":@:B:o:g:s:t:l:5:3:n:b:m:a:ANcduPpvh"))>=0) {
        switch (c) {
            case 'B': snp_bed_fn = optarg; break;
            case 'o': outfn = optarg; break;
            case 'g': reg = optarg; break;
            case '@': conf.n_threads = atoi(optarg); break;
            case 's': conf.step = atoi(optarg); break;
            case 't': conf.max_retention = atoi(optarg); break;
            case 'l': conf.min_read_len = atoi(optarg); break;
            case '5': conf.min_dist_end_5p = atoi(optarg); break;
            case '3': conf.min_dist_end_3p = atoi(optarg); break;
            case 'n': conf.max_nm = atoi(optarg); break;
            case 'b': conf.min_base_qual = atoi(optarg); break;
            case 'm': conf.min_mapq = atoi(optarg); break;
            case 'a': conf.min_score = atoi(optarg); break;
            case 'A': conf.print_all_locations = 1; break;
            case 'N': conf.is_nome = 1; break;
            case 'c': conf.filter_secondary = 0; break;
            case 'd': conf.filter_doublecnt = 0; break;
            case 'u': conf.filter_duplicate = 0; break;
            case 'p': conf.filter_ppair = 0; break;
            case 'P': conf.epiread_pair = 1; break;
            case 'v': conf.verbose = 1; break;
            case 'h': return usage(&conf);
            case ':': usage(&conf); wzfatal("Option needs an argument: -%c\n", optopt);
            case '?': usage(&conf); wzfatal("Unrecognized option: -%c\n", optopt);
            default: return usage(&conf);
        }
    }

    if (optind + 2 > argc) {
        usage(&conf);
        wzfatal("Reference or bam input is missing\n");
    }
    char *reffn = argv[optind++];
    char *infn = argv[optind++];

    episnp_chrom1_v *episnp = snp_bed_fn ?
        bed_init_episnp(snp_bed_fn) : NULL;

    wqueue_t(window) *wq = wqueue_init(window, 100000);
    pthread_t *processors = calloc(conf.n_threads, sizeof(pthread_t));
    result_t *results = calloc(conf.n_threads, sizeof(result_t));
    int i; unsigned j;
    htsFile *in = hts_open(infn, "rb");
    bam_hdr_t *header = sam_hdr_read(in);

    // sort sequence name by alphabetic order, chr1, chr10, chr11 ...
    target_v *targets = init_target_v(50);
    target_t *t;
    for (i=0; i<header->n_targets; ++i) {
        t = next_ref_target_v(targets);
        t->tid = i;
        t->name = header->target_name[i];
        t->len = header->target_len[i];
    }

    qsort(targets->buffer, targets->size,
            sizeof(target_t), compare_targets);

    // setup writer
    pthread_t writer;
    writer_conf_t writer_conf = {
        .q = wqueue_init(record, 100000),
        .outfn = outfn,
        .statsfn = statsfn,
        .header = 0,
        .targets = targets,
        .conf = &conf,
    };
    pthread_create(&writer, NULL, epiread_write_func, &writer_conf);
    for (i=0; i<conf.n_threads; ++i) {
        results[i].q = wq;
        results[i].rq = writer_conf.q;
        results[i].snp = episnp;
        results[i].ref_fn = reffn;
        results[i].bam_fn = infn;
        results[i].conf = &conf;
        pthread_create(&processors[i], NULL, process_func, &results[i]);
    }

    window_t w; memset(&w, 0, sizeof(window_t));
    uint32_t wbeg;
    int64_t block_id=0;

    // process bam
    if (reg) { // regional
        int tid;
        uint32_t beg, end;
        pileup_parse_region(reg, header, &tid, (int*) &beg, (int*) &end);
        // chromosome are assumed to be less than 2**29
        beg++; end++;
        if (beg<=0) beg = 1;
        if (end>header->target_len[tid]) end = header->target_len[tid];
        for (wbeg = beg; wbeg < end; wbeg += conf.step, block_id++) {
            w.tid = tid;
            w.block_id = block_id;
            w.beg = wbeg;
            w.end = wbeg + conf.step;
            if (w.end > end) w.end = end;
            wqueue_put(window, wq, &w);
        }
    } else { // entire bam
        for (j=0; j<targets->size; ++j) {
            t = ref_target_v(targets, j);
            for (wbeg = 1; wbeg < t->len; wbeg += conf.step, block_id++) {
                w.tid = t->tid;
                w.block_id = block_id;
                w.beg = wbeg;
                w.end = wbeg+conf.step;
                if (w.end > t->len) w.end = t->len;
                wqueue_put(window, wq, &w);
            }
        }
    }
    for (i=0; i<conf.n_threads; ++i) {
        w.tid = -1;
        wqueue_put(window, wq, &w);
    }

    for (i=0; i<conf.n_threads; ++i) {
        pthread_join(processors[i], NULL);
    }

    record_t rec = { .block_id = RECORD_QUEUE_END };
    wqueue_put2(record, writer_conf.q, rec);
    pthread_join(writer, NULL);
    wqueue_destroy(record, writer_conf.q);

    free_target_v(targets);
    free(results);
    free(processors);
    wqueue_destroy(window, wq);
    hts_close(in);
    bam_hdr_destroy(header);

    return 0;
}
