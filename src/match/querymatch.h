/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef QUERYMATCH_H
#define QUERYMATCH_H

#include <inttypes.h>
#include <stdio.h>
#include "core/error_api.h"
#include "core/readmode.h"
#include "core/encseq.h"
#include "core/unused_api.h"
#include "core/arraydef.h"
#include "querymatch-align.h"
#include "karlin_altschul_stat.h"
#include "querymatch-display.h"
#include "seq_or_encseq.h"

typedef struct GtQuerymatch GtQuerymatch;

GT_DECLAREARRAYSTRUCT(GtQuerymatch);

GtQuerymatch *gt_querymatch_new(void);

void gt_querymatch_file_set(GtQuerymatch *querymatch, FILE *fp);

void gt_querymatch_db_keyvalues_set(GtQuerymatch *querymatch,
                                    GtUword db_totallength,
                                    GtUword db_numofsequences);

void gt_querymatch_init(GtQuerymatch *querymatch,
                        const GtKarlinAltschulStat *karlin_altschul_stat,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtUword dbseqlen,
                        GtWord score,
                        GtUword distance,
                        GtUword mismatches,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength,
                        const char *db_desc,
                        const char *query_desc);

bool gt_querymatch_read_line(GtQuerymatch *querymatchptr,
                             bool withseqlength,
                             const char *line_ptr,
                             bool selfmatch,
                             GtUword seedpos1,
                             GtUword seedpos2,
                             GtUword seedlen,
                             const GtEncseq *dbencseq,
                             const GtEncseq *queryencseq);

bool gt_querymatch_process(GtQuerymatch *querymatchptr,
                           const GtKarlinAltschulStat *karlin_altschul_stat,
                           const GtSeqorEncseq *dbes,
                           const GtSeqorEncseq *queryes,
                           bool greedyextension);

void gt_querymatch_delete(GtQuerymatch *querymatch);

bool gt_querymatch_complete(GtQuerymatch *querymatchptr,
                            const GtKarlinAltschulStat *karlin_altschul_stat,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtUword dbseqlen,
                            GtWord score,
                            GtUword distance,
                            GtUword mismatches,
                            bool selfmatch,
                            uint64_t queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtSeqorEncseq *dbes,
                            const GtSeqorEncseq *queryes,
                            GtUword query_totallength,
                            GtUword seedpos1,
                            GtUword seedpos2,
                            GtUword seedlen,
                            bool greedyextension);

void gt_querymatch_outoptions_set(GtQuerymatch *querymatch,
                GtQuerymatchoutoptions *querymatchoutoptions);

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch);

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch);

const GtUchar *gt_querymatch_querysequence(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dbseqnum(const GtQuerymatch *querymatch);

GtUword gt_querymatch_dblen(const GtQuerymatch *querymatch);

bool gt_querymatch_queryreverse(const GtQuerymatch *querymatch);

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen);

GtUword gt_querymatch_distance(const GtQuerymatch *querymatch);

GtWord gt_querymatch_score(const GtQuerymatch *querymatch);

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch);

void gt_querymatch_enhanced_prettyprint(double evalue,double bit_score,
                                        const GtQuerymatch *querymatch);

void gt_querymatch_coordinates_out(const GtQuerymatch *querymatch);

bool gt_querymatch_check_final(double *evalue_ptr,
                               double *bit_score_ptr,
                               const GtQuerymatch *querymatch,
                               GtUword userdefinedleastlength,
                               GtUword errorpercentage,
                               double evalue_threshold);

void gt_querymatch_query_readmode_set(GtQuerymatch *querymatch,
                                      GtReadmode query_readmode);

void gt_querymatch_verify_alignment_set(GtQuerymatch *querymatch);

GtReadmode gt_querymatch_query_readmode(const GtQuerymatch *querymatch);

GtWord gt_querymatch_distance2score(GtUword distance,GtUword alignedlen);

/* Returns true, iff the given seed start position is below the querymatch. */
bool gt_querymatch_overlap(const GtQuerymatch *querymatch,
                           GtUword nextseed_db_end_relative,
                           GtUword nextseed_query_end_relative,
                           bool use_db_pos);

bool gt_querymatch_previousmatches_overlap(GtQuerymatch *querymatch,
                            GtUword nextseed_db_end_relative,
                            GtUword nextseed_query_end_relative);

void gt_querymatch_previousmatches_add(GtQuerymatch *querymatch);

void gt_querymatch_previousmatches_clear(GtQuerymatch *querymatch);

bool gt_querymatch_has_seed(const GtQuerymatch *querymatch);

GtUword gt_querymatch_querystart_fwdstrand(const GtQuerymatch *querymatch);

void gt_querymatch_table_add(GtArrayGtQuerymatch *querymatch_table,
                             const GtQuerymatch *querymatch);

void gt_querymatch_table_sort(GtArrayGtQuerymatch *querymatch_table,
                              bool ascending);

GtQuerymatch *gt_querymatch_table_get(const GtArrayGtQuerymatch
                                        *querymatch_table,GtUword idx);

void gt_querymatch_display_set(GtQuerymatch *querymatch,
                               const GtSeedExtendDisplayFlag *display_flag);

#endif
