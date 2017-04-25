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
#define GT_DIAGBANDSEED_DIAGONALBAND(AMAXLEN,LOGDIAGBANDWIDTH,APOS,BPOS)\
        (GT_DIAGBANDSEED_DIAGONAL(AMAXLEN,APOS,BPOS) >> (LOGDIAGBANDWIDTH))

#define SIZE 100

GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth)
{
  return 1 + ((amaxlen + bmaxlen) >> logdiagbandwidth);
}

typedef uint32_t GtDiagbandseedScore;

struct GtDiagbandStruct
{
  GtUword amaxlen, logdiagbandwidth, num_diagbands, used_diagbands,
          reset_from_matches, reset_with_memset;
  GtDiagbandseedScore *score;
  GtDiagbandseedPosition *lastpos;
};

typedef struct SeedArea{
  GtUword abegin, aend, bbegin, bend, diagidx; 
} SeedArea;

GtUword amaxlen, logdiagbandwidth; //TODO: GLOBALE VARIABLE, noch aendern!!!

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
  diagband_struct->amaxlen = amaxlen;
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

void gt_diagband_struct_reset_counts(const GtDiagbandStruct *diagband_struct,
                                     FILE *stream)
{
  fprintf(stream,"# number of resets of all used diagonal bands: " GT_WU,
          diagband_struct->reset_with_memset +
          diagband_struct->reset_from_matches);
  if (diagband_struct->reset_with_memset > 0)
  {
    fprintf(stream,"; simple resets: " GT_WU,
            diagband_struct->reset_with_memset);
  }
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

void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                      GtDiagbandseedPosition apos,
                                      GtDiagbandseedPosition bpos,
                                      GtDiagbandseedPosition matchlength)
{
  GtUword diagband_idx;

  gt_assert(diagband_struct != NULL);
  diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                              diagband_struct->logdiagbandwidth,
                                              apos,
                                              bpos);
  gt_assert(diagband_idx < diagband_struct->num_diagbands);
  if (diagband_struct->lastpos[diagband_idx] == 0 /* first matches */||
      /* match with end position bpos begins strictly after previous match */
      diagband_struct->lastpos[diagband_idx] + matchlength <= bpos)
  {
    /* no overlap */
    diagband_struct->lastpos[diagband_idx] = bpos;
    if (diagband_struct->score[diagband_idx] == 0)
    {
      diagband_struct->used_diagbands++;
    }
    diagband_struct->score[diagband_idx] += matchlength;
  } else
  {
    /* overlap: add positions after last counted position */
    if (diagband_struct->lastpos[diagband_idx] < bpos)
    {
      const GtUword addlength = bpos - diagband_struct->lastpos[diagband_idx];

      diagband_struct->lastpos[diagband_idx] = bpos;
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
  diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                              diagband_struct->logdiagbandwidth,
                                              apos, bpos);
  return gt_diagband_struct_dband_coverage(diagband_struct,diagband_idx);
}

void gt_diagband_struct_multi_update(GtDiagbandStruct *diagband_struct,
                                     const GtDiagbandseedMaximalmatch *memstore,
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
          = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                         diagband_struct->logdiagbandwidth,
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
          = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                         diagband_struct->logdiagbandwidth,
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

static int compare_seeds_by_bstart(const void *p1, const void *p2)
{
  const SeedArea *area1 = (const SeedArea *) p1;
  const SeedArea *area2 = (const SeedArea *) p2;
  
  if(area1->bbegin < area2->bbegin)
    return -1;
    
  if(area1->bbegin > area2->bbegin)
    return 1;
  
  //TODO: im else fall vll area kleiner einstufen, die naeher an der hauptdiagonale ist
  //statt kleinerer index? oder laengeres?
  
  if(area1->diagidx < area2->diagidx)
    return -1;
    
  if(area1->diagidx > area2->diagidx)
    return 1;
  
  return 0;

}


static int compare_seeds_by_diags(const void *p1, const void *p2)
{

  const GtSeedpairPositions *seed1 = (const GtSeedpairPositions*)p1;
  const GtSeedpairPositions *seed2 = (const GtSeedpairPositions*)p2;

  const GtUword diagband_idx1
  = GT_DIAGBANDSEED_DIAGONALBAND(amaxlen,
                                 logdiagbandwidth,
                                 seed1->apos,
                                 seed1->bpos);

    const GtUword diagband_idx2
  = GT_DIAGBANDSEED_DIAGONALBAND(amaxlen,
                                 logdiagbandwidth,
                                 seed2->apos,
                                 seed2->bpos);
  
  if (diagband_idx1>diagband_idx2)
    return 1;
  
  if (diagband_idx1<diagband_idx2)
    return -1;
    
  if (seed1->bpos < seed2->bpos)
    return -1;
    
  if (seed1->bpos > seed2->bpos)
    return 1;
  
  return 0;
}

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                /* remove GT_UNUSED once arguments are used */
                                const GtEncseq *aencseq,
                                const GtEncseq *bencseq,
                                GT_UNUSED GtUword aseqnum,
                                GT_UNUSED GtUword bseqnum,
                                const GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                unsigned int seedlength,
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
    //TODO: ACHTUNG GLOBALE VARIABLE, noch aendern, benoetigt fuer qsort!!!
    amaxlen = diagband_struct->amaxlen;
    logdiagbandwidth = diagband_struct->logdiagbandwidth;

    //for (idx = 0; idx < segment_length; idx++)
    //{
      ////printf("apos: %d, bpos: %d\n", seedstore[idx].apos, seedstore[idx].bpos);
      //const GtUword diagband_idx
        //= GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                       //diagband_struct->logdiagbandwidth,
                                       //seedstore[idx].apos,
                                       //seedstore[idx].bpos);
      ////gt_diagband_statistics_score_add(diagband_statistics,diagband_struct,
      ////                                 diagband_idx);
    //}
    GtUword curr_diagband_idx, prev_diagband_idx, count=0, size=0;
    GtSeedpairPositions *seedstore2 = gt_malloc(sizeof(*seedstore2)*segment_length);
    seedstore2 = memcpy(seedstore2,seedstore, segment_length*sizeof(*seedstore));
    
    //TODO: struct von pointern auf GtSeedpairPositions statt memcpy?

    //TODO: use radixsort?
    qsort(seedstore2, segment_length, sizeof(GtSeedpairPositions), compare_seeds_by_diags);

    SeedArea *seed_area = gt_malloc(diagband_struct->used_diagbands*sizeof(*seed_area));
    size=diagband_struct->used_diagbands;

    
    
    prev_diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                       diagband_struct->logdiagbandwidth,
                                       seedstore2[0].apos,
                                       seedstore2[0].bpos);
    
    seed_area[0].abegin = seedstore2[0].apos-seedlength;
    seed_area[0].bbegin = seedstore2[0].bpos-seedlength;
    seed_area[0].aend = seedstore2[0].apos;
    seed_area[0].bend = seedstore2[0].bpos;
    seed_area[0].diagidx = prev_diagband_idx;
    count++;
    
    
      
    for (idx = 1; idx < segment_length; idx++)
    {
      
        curr_diagband_idx = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                       diagband_struct->logdiagbandwidth,
                                       seedstore2[idx].apos,
                                       seedstore2[idx].bpos);
        
        if (curr_diagband_idx == prev_diagband_idx)
        {
          
          
          GtUword bbegin = seedstore2[idx].bpos-seedlength;
          if (seed_area[idx-1].bend +1 <= bbegin)
          {
            if(seed_area[idx-1].bend < seedstore2[idx].bpos)
            {
              seed_area[count].aend = seedstore2[idx].apos;
              seed_area[count].bend = seedstore2[idx].bpos;
            }
          }
          else //TODO: else faelle verknuepfen
          {
            if (count == size)
            {
              size += SIZE;
              seed_area = gt_realloc(seed_area, sizeof(*seed_area)*size);
            }
            seed_area[count].abegin = seedstore2[idx].apos-seedlength+1;
            seed_area[count].bbegin = seedstore2[idx].bpos-seedlength+1;
            seed_area[count].aend = seedstore2[idx].apos;
            seed_area[count].bend = seedstore2[idx].bpos;
            seed_area[count].diagidx = curr_diagband_idx;
            count++;

          }
          
        }
        else
        {//TODO: entferne area wenn nur ein einzelnes kmer auf der diagonalen liegt
          
          if (count == size)
          {
            size += SIZE;
            seed_area = gt_realloc(seed_area, sizeof(*seed_area)*size);
          }
          
          seed_area[count].abegin = seedstore2[idx].apos-seedlength+1;
          seed_area[count].bbegin = seedstore2[idx].bpos-seedlength+1;
          seed_area[count].aend = seedstore2[idx].apos;
          seed_area[count].bend = seedstore2[idx].bpos;
          seed_area[count].diagidx = curr_diagband_idx;
          count++;
        }
  
      }
      gt_free(seedstore2);

      // sort seed_areas by start in bseq
      qsort(seed_area, count, sizeof(SeedArea), compare_seeds_by_bstart);

      SeedArea *pre_area = &seed_area[0];
      GtUword nextstart_min= pre_area->bend+1,
      mismatches=0,
      len =0, jdx;
      bool anchor = false;

      // traverse all seed_areas
      //TODO: min length?
      for(idx=1; idx<count; idx++)
      {
        if (seed_area[idx].bbegin < nextstart_min)
          continue;
        
        if (seed_area[idx].diagidx != pre_area->diagidx)
        {
          if(anchor)
            len += pre_area->bend-pre_area->bbegin+1;
          pre_area = &seed_area[idx];
          nextstart_min = pre_area->bend+1;
          anchor = false;
        }
        else /* find two seed_areas on the same diagonal, current area starts
                earliest at prev_area+1*/
        {
          GtUword range = seed_area[idx].bbegin-pre_area->bend+1;
          for (jdx = 0; jdx < range; jdx++)
          {
            len++;
            
            if (gt_encseq_get_encoded_char(aencseq, pre_area->aend+1+jdx,GT_READMODE_FORWARD) !=
            gt_encseq_get_encoded_char(bencseq, pre_area->bend+1+jdx,GT_READMODE_FORWARD))
            {
              mismatches++;
            }
          }
          pre_area = &seed_area[idx];
          nextstart_min = pre_area->bend+1;
          anchor=true;
        }
      }
      
      printf("substitutions per site: %f\n", mismatches/(float)len);
      gt_free(seed_area);
  
  } else
  {
    gt_assert(memstore != NULL);
    for (idx = 0; idx < segment_length; idx++)
    {
      const GtUword diagband_idx
        = GT_DIAGBANDSEED_DIAGONALBAND(diagband_struct->amaxlen,
                                       diagband_struct->logdiagbandwidth,
                                       memstore[idx].apos,
                                       memstore[idx].bpos);
      gt_diagband_statistics_score_add(diagband_statistics,diagband_struct,
                                       diagband_idx);
    }
  }
}
