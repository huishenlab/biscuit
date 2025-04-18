/* Correct bisulfite strand information if it is very inconsistent with C2T/G2A count 
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2020 Wanding.Zhou@vai.org
 *               2021-2025 Jacob.Morrison@vai.org
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

#include "bsstrand.h"

conversion_tag_t bam_tag_get_bsstrand(bam1_t *b) {
    char *s;

    s = (char*) bam_aux_get(b, "ZS"); /* bsmap flag */
    if (s) {
        s++;
        if (*s == '+') return TAG_BSW;
        else if (*s == '-') return TAG_BSC;
    }

    s = (char*) bam_aux_get(b, "YD");     /* bwa-meth flag */
    if (s) {
        s++;
        if (*s == 'f') return TAG_BSW;
        if (*s == 'r') return TAG_BSC;
        if (*s == 'c') return TAG_CONFLICT;
        if (*s == 'u') return TAG_UNKNOWN;
    }

    s = (char*) bam_aux_get(b, "XG");     /* bismark flag */
    if (s) {
        s++;
        if (strcmp((char*)s, "CT")==0) return TAG_BSW;
        else if (strcmp((char*)s, "GA")==0) return TAG_BSC;
    }

    /* otherwise, guess the bsstrand from nCT and nGA */
    return TAG_UNKNOWN;
}

int bsstrand_func(bam1_t *b, samFile *out, bam_hdr_t *header, void *data) {

    bsstrand_data_t *d = (bsstrand_data_t*)data;
    bsstrand_conf_t *conf = d->conf;
    const bam1_core_t *c = &b->core;

    if (c->flag & BAM_FUNMAP){
        if (out) 
            if (sam_write1(out, header, b) < 0)
                wzfatal("Cannot write bam.\n");
        d->n_unmapped++;
        return 0;
    }

    refcache_fetch(d->rs, header->target_name[c->tid], c->pos, bam_endpos(b)+1);
    uint32_t i, rpos=c->pos+1, qpos=0;
    int32_t nC2T = 0, nG2A = 0;
    uint32_t j;
    char rbase, qbase;

    for (i=0; i<c->n_cigar; ++i) {
        uint32_t op = bam_cigar_op(bam_get_cigar(b)[i]);
        uint32_t oplen = bam_cigar_oplen(bam_get_cigar(b)[i]);
        switch(op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for(j=0; j<oplen; ++j) {
                    rbase = toupper(refcache_getbase(d->rs, rpos+j));
                    qbase = bscall(b, qpos+j);
                    if (rbase == 'C' && qbase == 'T') nC2T += 1;
                    if (rbase == 'G' && qbase == 'A') nG2A += 1;
                    /* printf("%c vs %c\n", toupper(rbase), qbase); */
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
            case BAM_CHARD_CLIP:
                qpos += oplen;
                break;
            default:
                fprintf(stderr, "Unknown cigar, %u\n", op);
                abort();
        }
    }

    /* inference */
    conversion_tag_t bsstrand;
    if (nC2T == 0 && nG2A == 0) {
        bsstrand = TAG_UNKNOWN;
    } else {
        double s = min(nG2A, nC2T)/max(nG2A, nC2T);
        if (nC2T > nG2A) {
            if (nG2A == 0 || s <= 0.5) {
                bsstrand = TAG_BSW;
            } else {
                bsstrand = TAG_CONFLICT;
            }
        } else {
            if (nC2T == 0 || s <= 0.5) {
                bsstrand = TAG_BSC;
            } else {
                bsstrand = TAG_CONFLICT;
            }
        }
    }

    conversion_tag_t tag = bam_tag_get_bsstrand(b);
    d->confusion[tag*4+bsstrand]++;

    /* inference compared to tag */
    if (conf->correct_bsstrand) {
        uint8_t *ytag = bam_aux_get(b, "YD");
        if (ytag) { /* store original to OY and update YD */
            ytag++;
            /* bam_aux_append(b, "OY", 'A', 1, ytag); */
            if (bsstrand != tag) {
                *ytag = conversion_tags[bsstrand];
                d->n_corr++;
            }
        } else { /* just append YD */
            uint8_t data = conversion_tags[bsstrand];
            bam_aux_append(b, "YD", 'A', 1, &data);
        }
    }

    // R1_FOR, R1_REV, R2_FOR, R2_REV
    d->strandcnt[(c->flag & BAM_FREAD1 ? 0 : 1) * 8 +
        (c->flag & BAM_FREVERSE ? 1 : 0) * 4 + tag]++;

    if (conf->output_count) {
        bam_aux_append(b, "YC", 'i', 4, (uint8_t*) &nC2T);
        bam_aux_append(b, "YG", 'i', 4, (uint8_t*) &nG2A);
    }

    if (out) 
        if (sam_write1(out, header, b) < 0)
            wzfatal("Cannot write bam.\n");
    d->n_mapped++;

    return 0;
}

static void usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: biscuit bsstrand [options] <ref.fa> <in.bam> [out.bam]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -g STR    Region (optional, will process the whole bam if not specified)\n");
    fprintf(stderr, "    -y        Append count of C>T (YC tag) and G>A (YG tag) in out.bam\n");
    fprintf(stderr, "    -c        Correct bsstrand in out.bam, YD tag will be replaced if it exists\n");
    fprintf(stderr, "                  and created if not\n");
    fprintf(stderr, "    -h        This help\n");
    fprintf(stderr, "\n");
}

/* if nC2T > 1 and nG2A > 1: then fail and mark "?"
   if nC2T - nG2A > 2 and "-" => "+"
   if nG2A - nC2T > 2 and "+" => "-"
   output a summary of how many reads are inconsistent */
int main_bsstrand(int argc, char *argv[]) {
    int c;
    char *reg = 0; /* region */
    bsstrand_conf_t conf = {0};

    kstring_t call = generate_command_line_string(argc, argv);

    if (argc < 2) { usage(); return 1; }
    while ((c = getopt(argc, argv, ":g:cyh")) >= 0) {
        switch (c) {
            case 'g': reg = optarg; break;
            case 'y': conf.output_count = 1; break;
            case 'c': conf.correct_bsstrand = 1; break;
            case 'h': usage(); return 1;
            case ':': usage(); wzfatal("Option needs an argument: -%c\n", optopt); break;
            case '?': usage(); wzfatal("Unrecognized option: -%c\n", optopt); break;
            default: usage(); return 1;
        }
    }
    char *reffn = optind < argc ? argv[optind++] : NULL;
    char *infn  = optind < argc ? argv[optind++] : NULL;
    char *outfn = optind < argc ? argv[optind++] : NULL;
    if (!reffn || !infn) {
        usage();
        wzfatal("Please provide reference and input bam.\n");
    }

    bsstrand_data_t d = {0}; // all counts reset to 0
    d.rs = init_refcache(reffn, 100, 100000);
    d.conf = &conf;
    bam_filter(infn, outfn, reg, &d, call.s, bsstrand_func);

    /*** output stats ***/
    fprintf(stderr, "Mapped reads: %d\n", d.n_mapped);
    fprintf(stderr, "Unmapped reads: %d\n", d.n_unmapped);
    fprintf(stderr, "Corrected reads: %d (%1.2f%%)\n", d.n_corr, (double)d.n_corr/(double)d.n_mapped*100.);

    int i;
    /* Mapping strand vs conversion strand */
    fprintf(stderr, "\nStrand Distribution:\n");
    fprintf(stderr, "strand\\BS      BSW (f)      BSC (r)\n");
    fprintf(stderr, "     R1 (f):   ");
    for (i=0;i<2;++i) { fprintf(stderr, "%-13d", d.strandcnt[i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "     R1 (r):   ");
    for (i=0;i<2;++i) { fprintf(stderr, "%-13d", d.strandcnt[4+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "     R2 (f):   ");
    for (i=0;i<2;++i) { fprintf(stderr, "%-13d", d.strandcnt[8+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "     R2 (r):   ");
    for (i=0;i<2;++i) { fprintf(stderr, "%-13d", d.strandcnt[12+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    for (i=0; i<2; ++i) {
        fprintf(stderr, "\nR%d mapped to OT/OB:   %d", i+1, 
                d.strandcnt[i*8+0*4+TAG_BSW] + d.strandcnt[i*8+1*4+TAG_BSC]);
        fprintf(stderr, "\nR%d mapped to CTOT/CTOB: %d", i+1,
                d.strandcnt[i*8+1*4+TAG_BSW] + d.strandcnt[i*8+0*4+TAG_BSC]);
    }
    fprintf(stderr, "\n");

    /* confusion matrix of conversion counts */
    fprintf(stderr, "\nConfusion counts (single-end):\n");
    fprintf(stderr, "orig\\infer      BSW (f)      BSC (r)      Conflict (c) Unknown (u)\n");
    fprintf(stderr, "     BSW (f):   ");
    for (i=0;i<4;++i) { fprintf(stderr, "%-13d", d.confusion[i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "     BSC (r):   ");
    for (i=0;i<4;++i) { fprintf(stderr, "%-13d", d.confusion[4+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "Conflict (c):   ");
    for (i=0;i<4;++i) { fprintf(stderr, "%-13d", d.confusion[8+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, " Unknown (u):   ");
    for (i=0;i<4;++i) { fprintf(stderr, "%-13d", d.confusion[12+i]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");

    free_refcache(d.rs);
    free(call.s);

    return 0;
}
