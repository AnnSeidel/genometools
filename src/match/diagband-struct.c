/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/str_api.h"
#include "core/minmax.h"
#include "core/intbits.h"
#include "match/diagband-struct.h"

/* called with real bpos */
GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  return 1 + ((amaxlen + bmaxlen) >> logdiagbandwidth);
}

typedef uint32_t GtDiagbandseedScore;

struct GtDiagbandStruct
{
  GtUword amaxlen, bmaxlen, logdiagbandwidth, num_diagbands, used_diagbands,
          reset_from_matches, reset_with_memset;
  GtDiagbandseedScore *score;
  GtDiagbandseedPosition *lastpos;
  bool bpos_sorted;
};

static GtUword gt_diagbandseed_diagonalband(const GtDiagbandStruct
                                              *diagband_struct,
                                            GtDiagbandseedPosition apos,
                                            GtDiagbandseedPosition bpos)
{
  if (diagband_struct->bpos_sorted)
  {
    return GT_DIAGBANDSEED_DIAGONAL(diagband_struct->amaxlen,apos,bpos)
           >> diagband_struct->logdiagbandwidth;
  }
  return GT_DIAGBANDSEED_DIAGONAL(diagband_struct->bmaxlen,bpos,apos)
         >> diagband_struct->logdiagbandwidth;
}

bool gt_diagband_struct_empty(const GtDiagbandStruct *diagband_struct)
{
  return diagband_struct->used_diagbands == 0 ? true : false;
}

GtDiagbandStruct *gt_diagband_struct_new(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  GtDiagbandStruct *diagband_struct = gt_malloc(sizeof *diagband_struct);

  diagband_struct->used_diagbands = 0;
  diagband_struct->num_diagbands
    = gt_diagband_struct_num_diagbands(amaxlen,bmaxlen,logdiagbandwidth);
  diagband_struct->bpos_sorted = true;
  diagband_struct->amaxlen = amaxlen;
  diagband_struct->bmaxlen = bmaxlen;
  diagband_struct->logdiagbandwidth = logdiagbandwidth;
  /* diagband_score[0] and diagband_score[num_diagbands+1] remain zero as
     boundaries */
  diagband_struct->score = gt_calloc(diagband_struct->num_diagbands + 2,
                                     sizeof *diagband_struct->score);
  diagband_struct->score++; /* so we need not increment the index when
                               accessing score */
  diagband_struct->lastpos
    = gt_calloc(diagband_struct->num_diagbands,
                sizeof *diagband_struct->lastpos);
  diagband_struct->reset_from_matches = 0;
  diagband_struct->reset_with_memset = 0;
  return diagband_struct;
}

void gt_diagband_struct_bpos_sorted_set(GtDiagbandStruct *diagband_struct,
                                        bool value)
{
  gt_assert(diagband_struct != NULL);
  diagband_struct->bpos_sorted = value;
}

void gt_diagbandseed_maxlen_update(GtDiagbandStruct *diagband_struct,
                                   GtUword amaxlen,GtUword bmaxlen)
{
  gt_assert(diagband_struct != NULL);
  diagband_struct->amaxlen = amaxlen;
  diagband_struct->bmaxlen = bmaxlen;
}

void gt_diagband_struct_reset_counts(const GtDiagbandStruct *diagband_struct,
                                     FILE *stream)
{
  fprintf(stream,"# number of resets of all used diagonal bands: " GT_WU,
          diagband_struct->reset_with_memset +
          diagband_struct->reset_from_matches);
  if (diagband_struct->reset_from_matches > 0)
  {
    fprintf(stream,"; resets from matches: " GT_WU,
            diagband_struct->reset_from_matches);
  }
  fprintf(stream,"\n");
}

void gt_diagband_struct_delete(GtDiagbandStruct *diagband_struct)
{
  if (diagband_struct != NULL)
  {
    diagband_struct->score--; /* need to recover original base adress */
    gt_free(diagband_struct->score);
    gt_free(diagband_struct->lastpos);
    gt_free(diagband_struct);
  }
}

/* for a given match of length <matchlength> ending a positions <apos> and
   <bpos> the sequence from A and from B, respectively, update the
   diagonal band score, which is the number of positions on B covered by the
   match. If previous matches have been added to the band before, then
   the positions overlapping with these on the B-sequence are not counted.*/

static void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                             GtDiagbandseedPosition apos,
                                             GtDiagbandseedPosition bpos,
                                             GtDiagbandseedPosition matchlength)
{
  GtUword diagband_idx, keypos;

  gt_assert(diagband_struct != NULL);
  diagband_idx = gt_diagbandseed_diagonalband(diagband_struct, apos, bpos);
  keypos = diagband_struct->bpos_sorted ? bpos : apos;
  gt_assert(diagband_idx < diagband_struct->num_diagbands);
  if (diagband_struct->lastpos[diagband_idx] == 0 /* first matches */||
      /* match with end position keypos begins strictly after previous match */
      diagband_struct->lastpos[diagband_idx] + matchlength <= keypos)
  {
    /* no overlap */
    diagband_struct->lastpos[diagband_idx] = keypos;
    if (diagband_struct->score[diagband_idx] == 0)
    {
      diagband_struct->used_diagbands++;
    }
    diagband_struct->score[diagband_idx] += matchlength;
  } else
  {
    /* overlap: add positions after last counted position */
    if (diagband_struct->lastpos[diagband_idx] < keypos)
    {
      const GtUword addlength = keypos - diagband_struct->lastpos[diagband_idx];

      diagband_struct->lastpos[diagband_idx] = keypos;
      if (diagband_struct->score[diagband_idx] == 0)
      {
        diagband_struct->used_diagbands++;
      }
      diagband_struct->score[diagband_idx] += addlength;
    }
  }
}

static GtUword gt_diagband_struct_dband_coverage(
                                    const GtDiagbandStruct *diagband_struct,
                                    GtUword diagband_idx)
{
  gt_assert(diagband_struct != NULL);
  return (GtUword) MAX(diagband_struct->score[diagband_idx + 1],
                       diagband_struct->score[diagband_idx - 1])
         + (GtUword) diagband_struct->score[diagband_idx];
}

GtUword gt_diagband_struct_coverage(const GtDiagbandStruct *diagband_struct,
                                    GtDiagbandseedPosition apos,
                                    GtDiagbandseedPosition bpos)
{
  GtUword diagband_idx;

  gt_assert(diagband_struct != NULL);
  diagband_idx = gt_diagbandseed_diagonalband(diagband_struct,apos,bpos);
  return gt_diagband_struct_dband_coverage(diagband_struct,diagband_idx);
}

void gt_diagband_struct_mem_multi_update(GtDiagbandStruct *diagband_struct,
                                         const GtDiagbandseedMaximalmatch
                                           *memstore,
                                         GtUword numofmatches)
{
  GtUword idx;

  gt_assert(memstore != NULL);
  for (idx = 0; idx < numofmatches; idx++)
  {
    gt_diagband_struct_single_update(diagband_struct,
                                     memstore[idx].apos,
                                     memstore[idx].bpos,
                                     memstore[idx].len);
  }
}

void gt_diagband_struct_seed_multi_update(GtDiagbandStruct *diagband_struct,
                                          const GtSeedpairPositions *seedstore,
                                          GtUword segment_length,
                                          GtUword seedlength)
{
  GtUword idx;

  gt_assert(seedstore != NULL);
  for (idx = 0; idx < segment_length; idx++)
  {
    gt_assert (idx == 0 ||
               (diagband_struct->bpos_sorted &&
                seedstore[idx-1].bpos <= seedstore[idx].bpos) ||
               (!diagband_struct->bpos_sorted &&
                seedstore[idx-1].apos <= seedstore[idx].apos));
    gt_diagband_struct_single_update(diagband_struct,
                                     seedstore[idx].apos,
                                     seedstore[idx].bpos,
                                     seedlength);
  }
}

void gt_diagband_struct_reset(GtDiagbandStruct *diagband_struct,
                              const GtSeedpairPositions *seedstore,
                              const GtDiagbandseedMaximalmatch *memstore,
                              GtUword segment_length)
{
  gt_assert(diagband_struct != NULL);
  if (diagband_struct->used_diagbands * 3 >= diagband_struct->num_diagbands)
  { /* >= 33% of diagbands are used */
    memset(diagband_struct->score,0,
           sizeof *diagband_struct->score * diagband_struct->num_diagbands);
    memset(diagband_struct->lastpos,0,
           sizeof *diagband_struct->lastpos * diagband_struct->num_diagbands);
    diagband_struct->reset_with_memset++;
  } else
  {
    GtUword idx;

    if (seedstore != NULL)
    {
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = gt_diagbandseed_diagonalband(diagband_struct,
                                         seedstore[idx].apos,
                                         seedstore[idx].bpos);
        diagband_struct->score[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    } else
    {
      gt_assert(memstore != NULL);
      for (idx = 0; idx < segment_length; idx++)
      {
        const GtUword diagband_idx
          = gt_diagbandseed_diagonalband(diagband_struct,
                                         memstore[idx].apos,
                                         memstore[idx].bpos);
        diagband_struct->score[diagband_idx] = 0;
        diagband_struct->lastpos[diagband_idx] = 0;
      }
    }
    diagband_struct->reset_from_matches++;
  }
  diagband_struct->used_diagbands = 0;
}

struct GtDiagbandStatistics
{
  bool compute_sum;
  bool forward;
  GtUword sumscore;
  GtBitsequence *track;
};

GtDiagbandStatistics *gt_diagband_statistics_new(const GtStr
                                                   *diagband_distance_arg,
                                                 bool forward)
{
  const char *arg = gt_str_get(diagband_distance_arg);
  GtDiagbandStatistics *diagband_statistics
    = gt_malloc(sizeof *diagband_statistics);

  diagband_statistics->forward = forward;
  diagband_statistics->compute_sum = false;
  if (strcmp(arg,"sum") == 0)
  {
    diagband_statistics->compute_sum = true;
  } else
  {
    gt_assert(false);
  }
  diagband_statistics->sumscore = 0;
  diagband_statistics->track = NULL;
  return diagband_statistics;
}

void gt_diagband_statistics_delete(GtDiagbandStatistics *diagband_statistics)
{
  if (diagband_statistics != NULL)
  {
    gt_free(diagband_statistics->track);
    gt_free(diagband_statistics);
  }
}

void gt_diagband_statistics_display(const GtDiagbandStatistics
                                           *diagband_statistics)
{
  gt_assert(diagband_statistics != NULL);
  if (diagband_statistics->compute_sum)
  {
    printf("# forward=%s, sum_diagband_score=" GT_WU "\n",
            diagband_statistics->forward ? "true" : "false",
            diagband_statistics->sumscore);
  } else
  {
    gt_assert(false);
  }
}

static void gt_diagband_statistics_score_add(
                                      GtDiagbandStatistics *diagband_statistics,
                                      const GtDiagbandStruct *diagband_struct,
                                      GtUword diagband_idx)
{
  if (!GT_ISIBITSET(diagband_statistics->track,diagband_idx))
  {
    diagband_statistics->sumscore += diagband_struct->score[diagband_idx];
    GT_SETIBIT(diagband_statistics->track,diagband_idx);
  }
}

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                GT_UNUSED bool bpos_sorted,
                                /* remove GT_UNUSED once arguments are used */
                                GT_UNUSED const GtEncseq *aencseq,
                                GT_UNUSED const GtEncseq *bencseq,
                                GT_UNUSED GtUword aseqnum,
                                GT_UNUSED GtUword bseqnum,
                                const GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                GT_UNUSED unsigned int seedlength,
                                const GtSeedpairPositions *seedstore,
                                GtUword segment_length)
{
  GtUword idx;
  GtDiagbandStatistics *diagband_statistics
    = (GtDiagbandStatistics *) v_diagband_statistics;

  if (diagband_statistics->track == NULL)
  {
    GT_INITBITTAB(diagband_statistics->track,diagband_struct->num_diagbands);
  } else
  {
    GT_CLEARBITTAB(diagband_statistics->track,diagband_struct->num_diagbands);
  }
  if (seedstore != NULL)
  {
    for (idx = 0; idx < segment_length; idx++)
    {
      const GtUword diagband_idx
        = gt_diagbandseed_diagonalband(diagband_struct,
                                       seedstore[idx].apos,
                                       seedstore[idx].bpos);
      gt_diagband_statistics_score_add(diagband_statistics,diagband_struct,
                                       diagband_idx);
    }
  } else
  {
    gt_assert(memstore != NULL);
    for (idx = 0; idx < segment_length; idx++)
    {
      const GtUword diagband_idx
        = gt_diagbandseed_diagonalband(diagband_struct,
                                       memstore[idx].apos,
                                       memstore[idx].bpos);
      gt_diagband_statistics_score_add(diagband_statistics,diagband_struct,
                                       diagband_idx);
    }
  }
}
