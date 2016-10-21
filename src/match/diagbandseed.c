/*
  Copyright (c) 2015-2016 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2015-2016 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015-2016 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "core/arraydef.h"
#include "core/codetype.h"
#include "core/complement.h"
#include "core/cstr_api.h"
#include "core/encseq.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/timer_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "match/declare-readfunc.h"
#include "match/diagbandseed.h"
#include "match/kmercodes.h"
#include "match/querymatch.h"
#include "match/querymatch-align.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"

#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

#define GT_DIAGBANDSEED_SEQNUM_UNDEF UINT_MAX
#define GT_DIAGBANDSEED_FMT          "in " GT_WD ".%06ld seconds.\n"

typedef uint32_t GtDiagbandseedPosition;
typedef uint32_t GtDiagbandseedSeqnum;
typedef uint32_t GtDiagbandseedScore;

typedef struct { /* 8 + 4 + 4 bytes */
  GtCodetype code;              /* only sort criterion */
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
} GtDiagbandseedKmerPos;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedKmerPos);
DECLAREBufferedfiletype(GtDiagbandseedKmerPos);
DECLAREREADFUNCTION(GtDiagbandseedKmerPos);

typedef struct
{ /* 4 + 4 + 4 + 4 bytes */
  GtDiagbandseedSeqnum bseqnum, /*  2nd important sort criterion */
                       aseqnum; /* most important sort criterion */
  GtDiagbandseedPosition apos,
                         bpos;  /*  3rd important sort criterion */
} GtDiagbandseedSeedPair;

struct GtDiagbandseedInfo
{
  const GtEncseq *aencseq;
  const GtEncseq *bencseq;
  const GtDiagbandseedExtendParams *extp;
  GtRange *seedpairdistance;
  GtUword maxfreq,
          memlimit,
          anumseqranges,
          bnumseqranges;
  unsigned int seedlength;
  GtDiagbandseedPairlisttype splt;
  bool norev,
       nofwd,
       verify,
       verbose,
       debug_kmer,
       debug_seedpair,
       use_kmerfile,
       trimstat_on;
};

struct GtDiagbandseedExtendParams
{
  GtUword errorpercentage,
          userdefinedleastlength,
          logdiagbandwidth,
          mincoverage,
          maxalignedlendifference,
          history_size,
          perc_mat_history,
          sensitivity,
          alignmentwidth;
  GtXdropscore xdropbelowscore;
  GtExtendCharAccess extend_char_access;
  unsigned int display_flag;
  double matchscore_bias;
  bool use_apos,
       extendgreedy,
       extendxdrop,
       weakends,
       benchmark,
       always_polished_ends,
       verify_alignment;
};

typedef struct
{
  GtArrayGtDiagbandseedKmerPos *list;
  GtDiagbandseedSeqnum seqnum;
  GtDiagbandseedPosition endpos;
  const GtEncseq *encseq;
  GtSpecialrangeiterator *sri;
  GtRange *specialrange;
  GtUword last_specialpos,
          prev_separator,
          next_separator;
  unsigned int seedlength;
  GtReadmode readmode;
} GtDiagbandseedProcKmerInfo;

/* * * * * CONSTRUCTORS AND DESTRUCTORS * * * * */

GtDiagbandseedInfo *gt_diagbandseed_info_new(const GtEncseq *aencseq,
                                             const GtEncseq *bencseq,
                                             GtUword maxfreq,
                                             GtUword memlimit,
                                             unsigned int seedlength,
                                             bool norev,
                                             bool nofwd,
                                             GtRange *seedpairdistance,
                                             GtDiagbandseedPairlisttype splt,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             bool trimstat_on,
                                             const GtDiagbandseedExtendParams
                                               *extp,
                                             GtUword anumseqranges,
                                             GtUword bnumseqranges)
{
  GtDiagbandseedInfo *info = gt_malloc(sizeof *info);
  info->aencseq = aencseq;
  info->bencseq = bencseq;
  info->maxfreq = maxfreq;
  info->memlimit = memlimit;
  info->seedlength = seedlength;
  info->norev = norev;
  info->nofwd = nofwd;
  info->seedpairdistance = seedpairdistance;
  info->splt = splt;
  info->verify = verify;
  info->verbose = verbose;
  info->debug_kmer = debug_kmer;
  info->debug_seedpair = debug_seedpair;
  info->use_kmerfile = use_kmerfile;
  info->extp = extp;
  info->anumseqranges = anumseqranges;
  info->bnumseqranges = bnumseqranges;
  info->trimstat_on = trimstat_on;
  return info;
}

void gt_diagbandseed_info_delete(GtDiagbandseedInfo *info)
{
  if (info != NULL) {
    gt_free(info);
  }
}

GtDiagbandseedExtendParams *gt_diagbandseed_extend_params_new(
                                GtUword errorpercentage,
                                GtUword userdefinedleastlength,
                                GtUword logdiagbandwidth,
                                GtUword mincoverage,
                                unsigned int display_flag,
                                bool use_apos,
                                GtXdropscore xdropbelowscore,
                                bool extendgreedy,
                                bool extendxdrop,
                                GtUword maxalignedlendifference,
                                GtUword history_size,
                                GtUword perc_mat_history,
                                GtExtendCharAccess extend_char_access,
                                GtUword sensitivity,
                                double matchscore_bias,
                                bool weakends,
                                bool benchmark,
                                GtUword alignmentwidth,
                                bool always_polished_ends,
                                bool verify_alignment)
{
  GtDiagbandseedExtendParams *extp = gt_malloc(sizeof *extp);
  extp->errorpercentage = errorpercentage;
  extp->userdefinedleastlength = userdefinedleastlength;
  extp->logdiagbandwidth = logdiagbandwidth;
  extp->mincoverage = mincoverage;
  extp->display_flag = display_flag;
  extp->use_apos = use_apos;
  extp->xdropbelowscore = xdropbelowscore;
  extp->extendgreedy = extendgreedy;
  extp->extendxdrop = extendxdrop;
  extp->maxalignedlendifference = maxalignedlendifference;
  extp->history_size = history_size;
  extp->perc_mat_history = perc_mat_history;
  extp->extend_char_access = extend_char_access;
  extp->sensitivity = sensitivity;
  extp->matchscore_bias = matchscore_bias;
  extp->weakends = weakends;
  extp->benchmark = benchmark;
  extp->alignmentwidth = alignmentwidth;
  extp->always_polished_ends = always_polished_ends;
  extp->verify_alignment = verify_alignment;
  return extp;
}

void gt_diagbandseed_extend_params_delete(GtDiagbandseedExtendParams *extp)
{
  if (extp != NULL) {
    gt_free(extp);
  }
}

/* * * * * K-MER LIST CREATION * * * * */

/* Estimate the number of k-mers in the given encseq. */
static GtUword gt_seed_extend_numofkmers(const GtEncseq *encseq,
                                         unsigned int seedlength,
                                         const GtSequenceRangeWithMaxLength
                                           *seqrange)
{
  GtUword lastpos, numofpos, subtract, ratioofspecial;

  const GtUword totalnumofspecial = gt_encseq_specialcharacters(encseq),
                totalnumofpos = gt_encseq_total_length(encseq),
                firstpos = gt_encseq_seqstartpos(encseq, seqrange->start),
                numofseq = seqrange->end - seqrange->start + 1;
  lastpos = (seqrange->end + 1 == gt_encseq_num_of_sequences(encseq)
             ? totalnumofpos
             : gt_encseq_seqstartpos(encseq, seqrange->end + 1) - 1);
  gt_assert(lastpos >= firstpos);
  numofpos = lastpos - firstpos;

  subtract = MIN(seedlength - 1, gt_encseq_min_seq_length(encseq)) + 1;
  gt_assert(numofpos + 1 >= numofseq * subtract);
  ratioofspecial = MIN(totalnumofspecial * numofpos / totalnumofpos, numofpos);
  return numofpos - MAX(numofseq * subtract - 1, ratioofspecial);
}

/* Returns the position of the next separator following specialrange.start.
   If the end of the encseq is reached, the position behind is returned. */
static GtUword gt_diagbandseed_update_separatorpos(GtRange *specialrange,
                                                   GtSpecialrangeiterator *sri,
                                                   const GtEncseq *encseq,
                                                   GtUword totallength,
                                                   GtReadmode readmode)
{
  gt_assert(sri != NULL && specialrange != NULL && encseq != NULL);
  do {
    GtUword idx;
    for (idx = specialrange->start; idx < specialrange->end; idx++) {
      if (gt_encseq_position_is_separator(encseq, idx, readmode)) {
        specialrange->start = idx + 1;
        return idx;
      }
    }
  } while (gt_specialrangeiterator_next(sri, specialrange));
  return totallength;
}

/* Add given code and its seqnum and position to a kmer list. */
static void gt_diagbandseed_processkmercode(void *prockmerinfo,
                                            bool firstinrange,
                                            GtUword startpos,
                                            GtCodetype code)
{
  const GtUword array_incr = 256;
  GtDiagbandseedProcKmerInfo *pkinfo;
  GtDiagbandseedKmerPos *kmerposptr = NULL;

  gt_assert(prockmerinfo != NULL);
  pkinfo = (GtDiagbandseedProcKmerInfo *) prockmerinfo;
  GT_GETNEXTFREEINARRAY(kmerposptr,
                        pkinfo->list,
                        GtDiagbandseedKmerPos,
                        array_incr + 0.2 *
                        pkinfo->list->allocatedGtDiagbandseedKmerPos);

  /* check separator positions and determine next seqnum and endpos */
  if (firstinrange) {
    const GtUword endpos = startpos + pkinfo->seedlength - 1;
    while (endpos >= pkinfo->next_separator) {
      pkinfo->seqnum++;
      pkinfo->prev_separator = pkinfo->next_separator + 1;
      pkinfo->next_separator
        = gt_diagbandseed_update_separatorpos(pkinfo->specialrange,
                                              pkinfo->sri,
                                              pkinfo->encseq,
                                              pkinfo->last_specialpos,
                                              pkinfo->readmode);
      gt_assert(pkinfo->next_separator >= pkinfo->prev_separator);
    }
    gt_assert(endpos >= pkinfo->prev_separator);
    gt_assert(startpos < pkinfo->next_separator);
    if (pkinfo->readmode == GT_READMODE_FORWARD) {
      pkinfo->endpos = (GtDiagbandseedPosition) (endpos -
                                                 pkinfo->prev_separator);
    } else {
      pkinfo->endpos = (GtDiagbandseedPosition) (pkinfo->next_separator - 1 -
                                                 startpos);
    }
  }

  /* save k-mer code */
  kmerposptr->code = (pkinfo->readmode == GT_READMODE_FORWARD
                      ? code : gt_kmercode_reverse(code, pkinfo->seedlength));
  /* save endpos and seqnum */
  gt_assert(pkinfo->endpos != UINT_MAX);
  kmerposptr->endpos = pkinfo->endpos;
  pkinfo->endpos = (pkinfo->readmode == GT_READMODE_FORWARD
                    ? pkinfo->endpos + 1 : pkinfo->endpos - 1);
  kmerposptr->seqnum = pkinfo->seqnum;
}

/* Uses GtKmercodeiterator for fetching the kmers. */
static void gt_diagbandseed_get_kmers_kciter(GtDiagbandseedProcKmerInfo *pkinfo)
{
  GtKmercodeiterator *kc_iter = NULL;
  const GtKmercode *kmercode = NULL;
  bool firstinrange = true;
  GtUword maxpos = 0, position;

  /* initialise GtKmercodeiterator */
  gt_assert(pkinfo != NULL);
  position = gt_encseq_seqstartpos(pkinfo->encseq, pkinfo->seqnum);
  kc_iter = gt_kmercodeiterator_encseq_new(pkinfo->encseq,
                                           pkinfo->readmode,
                                           pkinfo->seedlength,
                                           position);
  if (pkinfo->seedlength <= pkinfo->last_specialpos) {
    maxpos = pkinfo->last_specialpos + 1 - pkinfo->seedlength;
  }

  /* iterate */
  while (position < maxpos) {
    kmercode = gt_kmercodeiterator_encseq_next(kc_iter);
    if (!kmercode->definedspecialposition) {
      gt_diagbandseed_processkmercode((void *) pkinfo,
                                      firstinrange,
                                      position,
                                      kmercode->code);
      firstinrange = false;
    } else {
      firstinrange = true;
    }
    position++;
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Return a sorted list of k-mers of given seedlength from specified encseq.
 * Only sequences in seqrange will be taken into account.
 * The caller is responsible for freeing the result. */
GtArrayGtDiagbandseedKmerPos gt_diagbandseed_get_kmers(
                                   const GtEncseq *encseq,
                                   unsigned int seedlength,
                                   GtReadmode readmode,
                                   const GtSequenceRangeWithMaxLength *seqrange,
                                   bool debug_kmer,
                                   bool verbose,
                                   GtUword known_size,
                                   FILE *stream)
{
  GtArrayGtDiagbandseedKmerPos list;
  GtDiagbandseedProcKmerInfo pkinfo;
  GtRange specialrange;
  GtTimer *timer = NULL;
  GtUword listlen = known_size;
  const GtUword totallength = gt_encseq_total_length(encseq);

  gt_assert(encseq != NULL);
  if (known_size > 0)
  {
    listlen = known_size;
  } else
  {
    listlen = gt_seed_extend_numofkmers(encseq, seedlength, seqrange);
  }
  if (verbose) {
    timer = gt_timer_new();
    fprintf(stream, "# Start fetching %u-mers (expect " GT_WU ")...\n",
            seedlength, listlen);
    gt_timer_start(timer);
  }

  GT_INITARRAY(&list, GtDiagbandseedKmerPos);
  GT_CHECKARRAYSPACEMULTI(&list, GtDiagbandseedKmerPos, listlen);

  pkinfo.list = &list;
  pkinfo.seqnum = seqrange->start;
  pkinfo.endpos = 0;
  pkinfo.encseq = encseq;
  pkinfo.seedlength = seedlength;
  pkinfo.readmode = readmode;
  if (seqrange->end + 1 == gt_encseq_num_of_sequences(encseq)) {
    pkinfo.last_specialpos = totallength;
  } else {
    /* start position of following sequence, minus separator position */
    pkinfo.last_specialpos
      = gt_encseq_seqstartpos(encseq, seqrange->end + 1) - 1;
  }
  pkinfo.prev_separator = gt_encseq_seqstartpos(encseq, seqrange->start);
  if (gt_encseq_has_specialranges(encseq)) {
    bool search = true;
    pkinfo.sri = gt_specialrangeiterator_new(encseq, true);
    while (search && gt_specialrangeiterator_next(pkinfo.sri, &specialrange)) {
      search = specialrange.end < pkinfo.prev_separator ? true : false;
    }
    specialrange.start = pkinfo.prev_separator;
    pkinfo.specialrange = &specialrange;
    pkinfo.next_separator
      = gt_diagbandseed_update_separatorpos(pkinfo.specialrange,
                                            pkinfo.sri,
                                            pkinfo.encseq,
                                            totallength,
                                            pkinfo.readmode);
  } else {
    pkinfo.sri = NULL;
    pkinfo.specialrange = NULL;
    pkinfo.next_separator = pkinfo.last_specialpos;
  }

  if (gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    /* Use fast access to encseq, requires 2bit-enc and absence of wildcards. */
    getencseqkmers_twobitencoding_slice(encseq,
                                        readmode,
                                        seedlength,
                                        seedlength,
                                        false,
                                        gt_diagbandseed_processkmercode,
                                        (void *) &pkinfo,
                                        NULL,
                                        NULL,
                                        pkinfo.prev_separator,
                                        pkinfo.last_specialpos);
  } else {
    /* Use GtKmercodeiterator for encseq access */
    gt_diagbandseed_get_kmers_kciter(&pkinfo);
  }
  if (gt_encseq_has_specialranges(encseq)) {
    gt_specialrangeiterator_delete(pkinfo.sri);
  }
  listlen = list.nextfreeGtDiagbandseedKmerPos;

  /* reduce size of array to number of entries */
  /* list.allocatedGtDiagbandseedKmerPos = listlen;
  gt_realloc(list.spaceGtDiagbandseedKmerPos,
             listlen * sizeof (GtDiagbandseedKmerPos)); */

  if (debug_kmer) {
    const GtDiagbandseedKmerPos *idx = list.spaceGtDiagbandseedKmerPos;
    const GtDiagbandseedKmerPos *end = idx + listlen;
    while (idx < end) {
      fprintf(stream, "# Kmer (" GT_LX ",%"PRIu32",%"PRIu32")\n",
              idx->code, idx->endpos, idx->seqnum);
      idx++;
    }
  }

  if (verbose) {
    fprintf(stream, "# ...found " GT_WU " %u-mers ", listlen, seedlength);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_start(timer);
  }

  /* sort list */
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)
                                   list.spaceGtDiagbandseedKmerPos,
                                   listlen);
  if (verbose) {
    fprintf(stream, "# ...sorted " GT_WU " %u-mers ", listlen, seedlength);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }

  return list;
}

/* * * * * SEEDPAIR LIST CREATION * * * * */

typedef struct {
  GtArrayGtDiagbandseedKmerPos segment;
  bool at_end;
  /* for list based iterator */
  const GtArrayGtDiagbandseedKmerPos *origin_list;
  const GtDiagbandseedKmerPos *listend;
  GtDiagbandseedKmerPos *listptr;
  /* for file based iterator */
  GtBufferedfile_GtDiagbandseedKmerPos kmerstream;
  GtDiagbandseedKmerPos buffer;
} GtDiagbandseedKmerIterator;

static void gt_diagbandseed_kmer_iter_reset(GtDiagbandseedKmerIterator *ki)
{
  gt_assert(ki != NULL);
  ki->at_end = false;
  if (ki->origin_list != NULL) { /* list based */
    ki->listptr = ki->origin_list->spaceGtDiagbandseedKmerPos;
    ki->segment.spaceGtDiagbandseedKmerPos = ki->listptr;
    if (ki->origin_list->nextfreeGtDiagbandseedKmerPos == 0) {
      ki->at_end = true;
    }
  } else { /* file based */
    int rval;
    ki->kmerstream.nextread = ki->kmerstream.nextfree = 0;
    rewind(ki->kmerstream.fp);
    rval = gt_readnextfromstream_GtDiagbandseedKmerPos(&ki->buffer,
                                                       &ki->kmerstream);
    if (rval != 1) {
      ki->at_end = true;
    }
  }
}

static GtDiagbandseedKmerIterator *gt_diagbandseed_kmer_iter_new_list(
                                      const GtArrayGtDiagbandseedKmerPos *list)
{
  GtDiagbandseedKmerIterator *ki = gt_malloc(sizeof *ki);
  GT_INITARRAY(&ki->segment, GtDiagbandseedKmerPos);
  gt_assert(list != NULL);
  ki->origin_list = list;
  ki->listend = list->spaceGtDiagbandseedKmerPos +
                list->nextfreeGtDiagbandseedKmerPos;
  gt_diagbandseed_kmer_iter_reset(ki);
  return ki;
}

static
GtDiagbandseedKmerIterator *gt_diagbandseed_kmer_iter_new_file(FILE *fp)
{
  GtDiagbandseedKmerIterator *ki = gt_malloc(sizeof *ki);
  GT_INITARRAY(&ki->segment, GtDiagbandseedKmerPos);
  ki->origin_list = NULL;
  ki->listend = ki->listptr = NULL;
  gt_assert(fp != NULL);
  ki->kmerstream.fp = fp;
  ki->kmerstream.bufferedfilespace = gt_malloc(FILEBUFFERSIZE *
                                               sizeof (GtDiagbandseedKmerPos));
  gt_diagbandseed_kmer_iter_reset(ki);
  return ki;
}

static void gt_diagbandseed_kmer_iter_delete(GtDiagbandseedKmerIterator *ki)
{
  if (ki != NULL) {
    if (ki->origin_list == NULL) { /* file based */
      gt_free(ki->kmerstream.bufferedfilespace);
      gt_fa_fclose(ki->kmerstream.fp);
      GT_FREEARRAY(&ki->segment, GtDiagbandseedKmerPos);
    }
    gt_free(ki);
  }
}

static const GtArrayGtDiagbandseedKmerPos *gt_diagbandseed_kmer_iter_next(
                                              GtDiagbandseedKmerIterator *ki)
{
  GtCodetype code;
  if (ki->at_end) {
    return NULL;
  }
  ki->segment.nextfreeGtDiagbandseedKmerPos = 0; /* reset segment list */

  if (ki->origin_list != NULL) { /* list based */
    code = ki->listptr->code;
    ki->segment.spaceGtDiagbandseedKmerPos = ki->listptr;
    /* add element to segment list until code differs */
    do {
      ki->listptr++;
    } while (ki->listptr < ki->listend && code == ki->listptr->code);
    ki->segment.nextfreeGtDiagbandseedKmerPos
      += (GtUword) (ki->listptr - ki->segment.spaceGtDiagbandseedKmerPos);
    if (ki->listptr >= ki->listend) {
      ki->at_end = true;
    }
  } else { /* file based */
    int rval;
    code = ki->buffer.code;
    /* fill segment list from file, stop when code changes */
    do {
      GT_STOREINARRAY(&ki->segment, GtDiagbandseedKmerPos,
                      ki->segment.allocatedGtDiagbandseedKmerPos * 0.2 + 128,
                      ki->buffer);
      rval = gt_readnextfromstream_GtDiagbandseedKmerPos(&ki->buffer,
                                                         &ki->kmerstream);
    } while (rval == 1 && code == ki->buffer.code);
    if (rval != 1) {
      ki->at_end = true;
    }
  }
  return &ki->segment;
}

/* Evaluate the results of the seed pair count histogram */
static GtUword gt_diagbandseed_processhistogram(GtUword *histogram,
                                                GtUword maxfreq,
                                                GtUword maxgram,
                                                GtUword memlimit,
                                                GtUword mem_used,
                                                bool alist_blist_id,
                                                size_t sizeofunit)
{
  /* calculate available memory, take 98% of memlimit */
  GtUword count = 0, frequency = 0, mem_avail = 0.98 * memlimit;

  if (mem_avail > mem_used) {
    mem_avail = (mem_avail - mem_used) / sizeofunit;
  } else {
    mem_avail = 0;
    maxfreq = 0;
  }

  /* there is enough free memory */
  if (mem_avail > 0) {
    /* count seed pairs until available memory reached */
    for (frequency = 1; frequency <= maxgram && count < mem_avail;
         frequency++) {
      count += histogram[frequency - 1];
    }
    if (count > mem_avail) {
      gt_assert(frequency >= 2);
      frequency -= 2;
      gt_assert(count >= histogram[frequency]);
      count -= histogram[frequency];
    } else if (frequency == maxgram + 1) {
      frequency = GT_UWORD_MAX;
    }
    maxfreq = MIN(maxfreq, frequency);
  }

  /* determine minimum required memory for error message */
  if (maxfreq <= 1 && alist_blist_id) {
    count = (histogram[0] + histogram[1]) * sizeofunit;
    count = (count + mem_used) / 0.98;
  } else if (maxfreq == 0) {
    count = histogram[0] * sizeofunit;
    count = (count + mem_used) / 0.98;
  }
  histogram[maxgram] = count;
  return maxfreq;
}

static void gt_diagbandseed_encode_seedpair(GtBitbuffer *bb,
                                            uint8_t *bytestring,
                                            GtUword bytestring_length,
                                            const uint32_t *seedpair_values,
                                            const GtBitcount_type *bits_tab,
                                            unsigned int components)
{
  int idx;
  GtUword bytestring_offset = 0;

  for (idx = 0; idx < components; idx++)
  {
    if (bits_tab[idx] > 0)
    {
      bytestring_offset
        = gt_bitbuffer_write_bytestring(bb,
                                        bytestring,
                                        bytestring_offset,
                                        bytestring_length,
                                        (GtUword) seedpair_values[idx],
                                        bits_tab[idx]);
    }
  }
  gt_bitbuffer_flush(false,bb,bytestring + bytestring_offset);
}

static void gt_diagbandseed_decode_seedpair(GtBitbuffer *bb,
                                            const uint32_t *seedpair_values,
                                            const GtBitcount_type *bits_tab,
                                            const uint8_t *bytestring,
                                            unsigned int components)
{
  int idx;
  GtUword bytestring_offset = 0;

  gt_bitbuffer_reset_for_read(bb);
  for (idx = 0; idx < components; idx++)
  {
    if (bits_tab[idx] > 0)
    {
      GtUword value;
      bytestring_offset = gt_bitbuffer_read_bytestring (bb,
                                                        &value,
                                                        bytestring,
                                                        bytestring_offset,
                                                        bits_tab[idx]);
      if (value != (GtUword) seedpair_values[idx])
      {
        char buffer[GT_INTWORDSIZE+1];
        gt_bitsequence_tostring(buffer,value);
        fprintf(stderr,"line %d: value = " GT_WU
                       "(%s) != %u = seedpair_value[%d]\n",
                       __LINE__,value,
                       buffer + (64 - bits_tab[idx]),
                       seedpair_values[idx],idx);
        exit(EXIT_FAILURE);
      }
    }
  }
}

static GtBitcount_type gt_diagbandseed_sum_bits(const GtBitcount_type *bits_tab,
                                                unsigned int components)
{
  unsigned int idx;
  GtBitcount_type sum = 0;

  for (idx = 0; idx < components; idx++)
  {
    sum += bits_tab[idx];
  }
  return sum;
}

const int idx_aseqnum = 0, idx_bseqnum = 1, idx_bpos = 2, idx_apos = 3;

static void showbytestring(FILE *fp,const uint8_t *bytestring,
                           GtUword bytestring_length)
{
  char buffer[GT_INTWORDSIZE+1];
  GtUword idx;

  for (idx = 0; idx < bytestring_length; idx++)
  {
    gt_bitsequence_tostring(buffer,bytestring[idx]);
    fprintf(fp,"%s ",buffer + 56);
  }
}

static int gt_diagbandseed_seeds_compare_bytestring(const uint8_t *previous,
                                                    const uint8_t *current,
                                                    GtUword bytestring_length)
{
  int ret = memcmp(previous,current,(size_t) bytestring_length);
  if (ret < 0)
  {
    return -1;
  }
  if (ret > 0)
  {
    return 1;
  }
  return 0;
}

static void gt_diagbandseed_seedpaircmp(const uint32_t *sp1,
                                        const uint32_t *sp2)
{
  int idx;

  for (idx = 0; idx < 4; idx++)
  {
    if (sp1[idx] != sp2[idx])
    {
      fprintf(stderr,"idx %d: sp1 = %u != %u = sp2\n",idx,sp1[idx],sp2[idx]);
      exit(EXIT_FAILURE);
    }
    gt_assert(sp1[idx] == sp2[idx]);
  }
}

static GtUword gt_diagbandseed_seed2GtUword(const uint32_t *seedpair_values,
                                            const GtBitcount_type *bits_tab)
{
  GtUword encoding = 0;
  int idx;
  GtBitcount_type remain = CHAR_BIT * (GtBitcount_type) sizeof encoding;

  for (idx = 0; idx < 4; idx++)
  {
    if (bits_tab[idx] > 0)
    {
      gt_assert(remain >= bits_tab[idx] &&
                (GtUword) seedpair_values[idx] <=
                (((GtUword) 1) << bits_tab[idx]) - 1);
      remain -= bits_tab[idx];
      encoding |= (((GtUword) seedpair_values[idx]) << remain);
    } else
    {
      gt_assert(seedpair_values[idx] == 0);
    }
  }
  return encoding;
}

static void gt_diagbandseed_GtUword2seed(uint32_t *seedpair_values,
                                         GtUword encoding,
                                         const GtBitcount_type *bits_tab)
{
  int idx;
  GtBitcount_type remain = CHAR_BIT * (GtBitcount_type) sizeof encoding;

  for (idx = 0; idx < 4; idx++)
  {
    if (bits_tab[idx] > 0)
    {
      gt_assert(remain >= bits_tab[idx]);
      remain -= bits_tab[idx];
      seedpair_values[idx] = (encoding >> remain) &
                              ((((GtUword) 1) << bits_tab[idx]) - 1);
    } else
    {
      seedpair_values[idx] = 0;
    }
  }
}

const char *gt_diagbandseed_splt_comment(void)
{
  return "specify type of pairlist, possible values are ulong and bytestring";
}

GtDiagbandseedPairlisttype gt_diagbandseed_splt_get(const char *splt_string,
                                                    GtError *err)
{
  if (strcmp(splt_string,"ulong") == 0)
  {
    return GT_DIAGBANDSEED_SPLT_ULONG;
  }
  if (strcmp(splt_string,"bytestring") == 0)
  {
    return GT_DIAGBANDSEED_SPLT_BYTESTRING;
  }
  if (strcmp(splt_string,"") == 0)
  {
    return GT_DIAGBANDSEED_SPLT_STRUCT;
  }
  gt_error_set(err,"illegal parameter for option -splt: %s",
                    gt_diagbandseed_splt_comment());
  return -1;
}

GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

typedef struct
{
  GtArrayGtDiagbandseedSeedPair *mlist_struct;
  GtArrayGtUword *mlist_ulong;
  GtDiagbandseedPairlisttype splt;
  GtUword mask_tab[4], a_bseqnum_mask, aseqrange_start, bseqrange_start;
  GtBitcount_type bits_tab[4];
  int shift_tab[4];
} GtSeedpairlist;

#define GT_DIAGBANDSEED_ENCODE_SEEDPAIR(ASEQNUM,BSEQNUM,BPOS,APOS)\
             (((GtUword) (ASEQNUM) << seedpairlist->shift_tab[idx_aseqnum]) | \
              ((GtUword) (BSEQNUM) << seedpairlist->shift_tab[idx_bseqnum]) | \
              ((GtUword) (BPOS) << seedpairlist->shift_tab[idx_bpos]) | \
              ((GtUword) (APOS) << seedpairlist->shift_tab[idx_apos]))

static GtSeedpairlist *gt_seedpairlist_new(GtDiagbandseedPairlisttype splt,
                        const GtSequenceRangeWithMaxLength *aseqrange,
                        const GtSequenceRangeWithMaxLength *bseqrange)
{
  GtSeedpairlist *seedpairlist = gt_malloc(sizeof *seedpairlist);
  GtBitcount_type bits_seedpair;
  size_t sum_bytes;
  int idx;
  const GtUword anumofseq = aseqrange->end - aseqrange->start + 1,
                amaxlen = aseqrange->max_length,
                bnumofseq = bseqrange->end - bseqrange->start + 1,
                bmaxlen = bseqrange->max_length;

  seedpairlist->bits_tab[idx_aseqnum]
    = (GtBitcount_type) gt_radixsort_bits(anumofseq);
  seedpairlist->bits_tab[idx_bseqnum]
    = (GtBitcount_type) gt_radixsort_bits(bnumofseq);
  seedpairlist->bits_tab[idx_bpos]
    = (GtBitcount_type) gt_radixsort_bits(bmaxlen);
  seedpairlist->bits_tab[idx_apos]
    = (GtBitcount_type) gt_radixsort_bits(amaxlen);
  gt_assert(seedpairlist->bits_tab[idx_apos] > 0);
  gt_assert(seedpairlist->bits_tab[idx_bpos] > 0);
  for (idx = 0, bits_seedpair = 0; idx < 4; idx++)
  {
    bits_seedpair += seedpairlist->bits_tab[idx];
  }
  sum_bytes = gt_radixsort_bits2bytes(bits_seedpair);
  if (sum_bytes <= sizeof (GtUword))
  {
    int shift = (int) (sizeof (GtUword) * CHAR_BIT);
    for (idx = 0; idx < 4; idx++)
    {
      shift -= seedpairlist->bits_tab[idx];
      seedpairlist->shift_tab[idx] = shift;
      seedpairlist->mask_tab[idx]
        = (((GtUword) 1) << seedpairlist->bits_tab[idx]) - 1;
    }
  }
  seedpairlist->a_bseqnum_mask
    = (((GtUword) 1) << (seedpairlist->bits_tab[idx_aseqnum] +
                         seedpairlist->bits_tab[idx_bseqnum])) - 1;
  seedpairlist->aseqrange_start = aseqrange->start;
  seedpairlist->bseqrange_start = bseqrange->start;
  seedpairlist->mlist_struct = NULL;
  seedpairlist->mlist_ulong = NULL;
  if (splt == GT_DIAGBANDSEED_SPLT_ULONG && sum_bytes <= sizeof (GtUword))
  {
    seedpairlist->mlist_ulong = gt_malloc(sizeof *seedpairlist->mlist_ulong);
    GT_INITARRAY(seedpairlist->mlist_ulong, GtUword);
    seedpairlist->splt = splt;
  } else
  {
    if (splt == GT_DIAGBANDSEED_SPLT_BYTESTRING)
    {
      gt_assert(false);
    } else
    {
      seedpairlist->splt = GT_DIAGBANDSEED_SPLT_STRUCT;
      seedpairlist->mlist_struct
        = gt_malloc(sizeof *seedpairlist->mlist_struct);
      GT_INITARRAY(seedpairlist->mlist_struct, GtDiagbandseedSeedPair);
    }
  }
  return seedpairlist;
}

static void gt_seedpairlist_show_bits(FILE *stream,
                                      const GtSeedpairlist *seedpairlist)
{
  GtBitcount_type bits_seedpair
    = gt_diagbandseed_sum_bits(seedpairlist->bits_tab,4);
  fprintf(stream,"# bits_seedpair=%d, bytes_seedpair=%u with ",
                    (int) bits_seedpair,
                    (unsigned int) gt_radixsort_bits2bytes(bits_seedpair));
  fprintf(stream,"aseqnum=%hu bits,",seedpairlist->bits_tab[idx_aseqnum]);
  fprintf(stream,"apos=%hu bits,",seedpairlist->bits_tab[idx_apos]);
  fprintf(stream,"bseqnum=%hu bits,",seedpairlist->bits_tab[idx_bseqnum]);
  fprintf(stream,"bpos=%hu bits\n",seedpairlist->bits_tab[idx_bpos]);
}

static size_t gt_seedpairlist_sizeofunit(GtDiagbandseedPairlisttype splt)
{
  if (splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return sizeof (GtDiagbandseedSeedPair);
  } else
  {
    if (splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      return sizeof (GtUword);
    } else
    {
      gt_assert(false);
    }
  }
}

static void gt_seedpairlist_init(GtSeedpairlist *seedpairlist,
                                 GtUword known_size)
{
  if (known_size > 0) {
    gt_assert(seedpairlist != NULL);
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
    {
      gt_assert(seedpairlist->mlist_struct != NULL);
      GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_struct,
                              GtDiagbandseedSeedPair, known_size);
    } else
    {
      if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
      {
        gt_assert(seedpairlist->mlist_ulong != NULL);
        GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_ulong,GtUword, known_size);
      } else
      {
        gt_assert(false);
      }
    }
  }
}

static void gt_seedpairlist_delete(GtSeedpairlist *seedpairlist)
{
  if (seedpairlist != NULL)
  {
    if (seedpairlist->mlist_struct != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_struct, GtDiagbandseedSeedPair);
      gt_free(seedpairlist->mlist_struct);
    }
    if (seedpairlist->mlist_ulong != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_ulong, GtUword);
      gt_free(seedpairlist->mlist_ulong);
    }
    gt_free(seedpairlist);
  }
}

static GtUword gt_seedpairlist_length(const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL);
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    gt_assert(seedpairlist->mlist_struct != NULL);
    return seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair;
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      gt_assert(seedpairlist->mlist_ulong != NULL);
      return seedpairlist->mlist_ulong->nextfreeGtUword;
    } else
    {
      gt_assert(false);
    }
  }
  return 0;
}

static const GtDiagbandseedSeedPair *gt_seedpairlist_mlist_struct(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL &&
            seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT &&
            seedpairlist->mlist_struct != NULL);
  return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair;
}

static const GtUword *gt_seedpairlist_mlist_ulong(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL &&
            seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG &&
            seedpairlist->mlist_ulong != NULL);
  return seedpairlist->mlist_ulong->spaceGtUword;
}

static void gt_seedpairlist_add(GtSeedpairlist *seedpairlist,
                                bool knownsize,
                                const GtSequenceRangeWithMaxLength *aseqrange,
                                const GtSequenceRangeWithMaxLength *bseqrange,
                                GtDiagbandseedSeqnum aseqnum,
                                GtDiagbandseedSeqnum bseqnum,
                                GtDiagbandseedPosition bpos,
                                GtDiagbandseedPosition apos)
{
  gt_assert(seedpairlist != NULL);
  gt_assert(aseqnum >= aseqrange->start && aseqnum <= aseqrange->end &&
            bseqnum >= bseqrange->start && bseqnum <= bseqrange->end &&
            bpos < bseqrange->max_length && apos < aseqrange->max_length);
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    GtDiagbandseedSeedPair *seedpair = NULL;
    if (knownsize)
    {
      gt_assert(seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair <
                seedpairlist->mlist_struct->allocatedGtDiagbandseedSeedPair);
      seedpair = seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair +
                 seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair++;
    } else
    {
      GT_GETNEXTFREEINARRAY(seedpair,
                            seedpairlist->mlist_struct,
                            GtDiagbandseedSeedPair,
                            256 + 0.2 *
                            seedpairlist->
                               mlist_struct->allocatedGtDiagbandseedSeedPair);
    }
    seedpair->aseqnum = aseqnum;
    seedpair->bseqnum = bseqnum;
    seedpair->bpos = bpos;
    seedpair->apos = apos;
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      const GtUword encoding
        = GT_DIAGBANDSEED_ENCODE_SEEDPAIR(aseqnum - aseqrange->start,
                                          bseqnum - bseqrange->start,
                                          bpos,apos);

      if (knownsize)
      {
        gt_assert(seedpairlist->mlist_ulong->nextfreeGtUword <
                  seedpairlist->mlist_ulong->allocatedGtUword);
        seedpairlist->mlist_ulong->spaceGtUword[
                   seedpairlist->mlist_ulong->nextfreeGtUword++] = encoding;
      } else
      {
        GT_STOREINARRAY(seedpairlist->mlist_ulong,
                        GtUword,
                        256 + 0.2 * seedpairlist->mlist_ulong->allocatedGtUword,
                        encoding);
      }
    } else
    {
      gt_assert(false);
    }
  }
}

static void gt_diagbandseed_seedpairlist_sort(GtSeedpairlist *seedpairlist)
{
  GtUword mlistlen = gt_seedpairlist_length(seedpairlist);

  if (mlistlen > 0)
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
    {
      GtDiagbandseedSeedPair *mlist;
      gt_assert(seedpairlist->mlist_struct != NULL);
      mlist = seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair;
      gt_radixsort_inplace_Gtuint64keyPair((Gtuint64keyPair *) mlist, mlistlen);
    } else
    {
      if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
      {
        gt_assert(seedpairlist->mlist_ulong != NULL);
        gt_radixsort_inplace_ulong(seedpairlist->mlist_ulong->spaceGtUword,
                                   mlistlen);
      } else
      {
        gt_assert(false);
      }
    }
  }
}

static GtUword gt_seedpairlist_a_bseqnum (const GtSeedpairlist *seedpairlist,
                                          GtUword encoding)
{
  if (seedpairlist->a_bseqnum_mask == 0)
  {
    return 0;
  }
  return (encoding >> seedpairlist->shift_tab[idx_bseqnum]) &
         seedpairlist->a_bseqnum_mask;
}

static GtUword gt_seedpairlist_extract(const GtSeedpairlist *seedpairlist,
                                       GtUword encoding,
                                       int compidx)
{
  return (encoding >> seedpairlist->shift_tab[compidx]) &
         seedpairlist->mask_tab[compidx];
}

static GtUword gt_seedpairlist_extract_idx(const GtSeedpairlist *seedpairlist,
                                           GtUword spidx,
                                           int compidx)
{
  gt_assert(seedpairlist != NULL &&
            seedpairlist->mlist_ulong != NULL &&
            seedpairlist->mlist_ulong->spaceGtUword != NULL);
  return gt_seedpairlist_extract(seedpairlist,
                                 seedpairlist->mlist_ulong->spaceGtUword[spidx],
                                 compidx);
}

static GtDiagbandseedSeqnum gt_seedpairlist_aseqnum(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx]
           .aseqnum;
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    GtDiagbandseedSeqnum aseqnum = (GtDiagbandseedSeqnum)
           (gt_seedpairlist_extract_idx(seedpairlist,spidx,idx_aseqnum) +
           seedpairlist->aseqrange_start);
    return aseqnum;
  }
  gt_assert(false);
  return 0;
}

static GtDiagbandseedSeqnum gt_seedpairlist_bseqnum(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx]
           .bseqnum;
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    return (GtDiagbandseedSeqnum)
           (gt_seedpairlist_extract_idx(seedpairlist,spidx,idx_bseqnum) +
            seedpairlist->bseqrange_start);
  }
  gt_assert(false);
  return 0;
}

static GtDiagbandseedPosition gt_seedpairlist_apos(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx].apos;
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    return (GtDiagbandseedPosition)
           gt_seedpairlist_extract_idx(seedpairlist,spidx,idx_apos);
  }
  gt_assert(false);
  return 0;
}

static GtDiagbandseedPosition gt_seedpairlist_bpos(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx].bpos;
  }
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    return (GtDiagbandseedPosition)
           gt_seedpairlist_extract_idx(seedpairlist,spidx,idx_bpos);
  }
  gt_assert(false);
  return 0;
}

static void gt_diagbandseed_show_seed(FILE *stream,
                                      const GtSeedpairlist *seedpairlist,
                                      GtUword spidx)
{
  fprintf(stream, "(%"PRIu32 ",%"PRIu32 ",%"PRIu32 ",%"PRIu32")",
                  gt_seedpairlist_aseqnum(seedpairlist,spidx),
                  gt_seedpairlist_bseqnum(seedpairlist,spidx),
                  gt_seedpairlist_apos(seedpairlist,spidx),
                  gt_seedpairlist_bpos(seedpairlist,spidx));
}

static int gt_diagbandseed_seeds_compare(const GtSeedpairlist *seedpairlist,
                                         const GtUword current)
{
  GtDiagbandseedSeqnum p_aseqnum = gt_seedpairlist_aseqnum(seedpairlist,
                                                           current - 1),
                       c_aseqnum = gt_seedpairlist_aseqnum(seedpairlist,
                                                           current),
                       p_bseqnum, c_bseqnum;
  GtDiagbandseedPosition p_bpos, c_bpos, p_apos, c_apos;
  if (p_aseqnum < c_aseqnum)
  {
    return -1;
  }
  if (p_aseqnum > c_aseqnum)
  {
    return 1;
  }
  p_bseqnum = gt_seedpairlist_bseqnum(seedpairlist,current-1);
  c_bseqnum = gt_seedpairlist_bseqnum(seedpairlist,current);
  if (p_bseqnum < c_bseqnum)
  {
    return -1;
  }
  if (p_bseqnum > c_bseqnum)
  {
    return 1;
  }
  p_bpos = gt_seedpairlist_bpos(seedpairlist,current-1);
  c_bpos = gt_seedpairlist_bpos(seedpairlist,current);
  if (p_bpos < c_bpos)
  {
    return -1;
  }
  if (p_bpos > c_bpos)
  {
    return 1;
  }
  p_apos = gt_seedpairlist_apos(seedpairlist,current-1);
  c_apos = gt_seedpairlist_apos(seedpairlist,current);
  if (p_apos < c_apos)
  {
    return -1;
  }
  if (p_apos > c_apos)
  {
    return 1;
  }
  return 0;
}

static void gt_diagbandseed_seedpairlist_out(FILE *stream,
                                             const GtSeedpairlist *seedpairlist)
{
  GtUword spidx, mlistlen = gt_seedpairlist_length(seedpairlist);

  for (spidx = 0; spidx < mlistlen; spidx++)
  {
    gt_assert(spidx == 0 ||
              gt_diagbandseed_seeds_compare(seedpairlist,spidx) < 0);
    fprintf(stream, "# SeedPair ");
    gt_diagbandseed_show_seed(stream,seedpairlist,spidx);
    fprintf(stream, "\n");
  }
}

void gt_diagbandseed_seedpairlist_encode_decode(
                     const GtSeedpairlist *seedpairlist,
                     const GtBitcount_type *bits_tab,
                     bool debug_seedpair)
{
  GtUword spidx, mlistlen = gt_seedpairlist_length(seedpairlist);
  /* in order of their priority for sorting */
  uint32_t seedpair_values[4], seedpair_values2[4];
  GtBitbuffer *bb_read = gt_bitbuffer_new();
  GtBitbuffer *bb_write = gt_bitbuffer_new();
  GtBitcount_type bits_seedpair = gt_diagbandseed_sum_bits(bits_tab,4);
  GtUword bytestring_length = gt_radixsort_bits2bytes(bits_seedpair);
  uint8_t *bytestring = gt_malloc(sizeof *bytestring * bytestring_length);
  uint8_t *bytestring2 = gt_malloc(sizeof *bytestring2 * bytestring_length);
  uint8_t *previous = bytestring2, *current = bytestring;

  for (spidx = 0; spidx < mlistlen; spidx++)
  {
    uint8_t *tmp;
    uint32_t aseqnum = gt_seedpairlist_aseqnum(seedpairlist,spidx),
             bseqnum = gt_seedpairlist_bseqnum(seedpairlist,spidx),
             apos, bpos;
    seedpair_values[idx_aseqnum] = aseqnum;
    seedpair_values[idx_bseqnum] = bseqnum;
    seedpair_values[idx_bpos] = bpos = gt_seedpairlist_bpos(seedpairlist,spidx);
    seedpair_values[idx_apos] = apos = gt_seedpairlist_apos(seedpairlist,spidx);
    if (bytestring_length <= sizeof (GtUword))
    {
      GtUword encoding = gt_diagbandseed_seed2GtUword(seedpair_values,bits_tab);
      GtUword encoding2 = GT_DIAGBANDSEED_ENCODE_SEEDPAIR(aseqnum,bseqnum,
                                                          bpos,apos);
      gt_assert(encoding == encoding2);
      gt_diagbandseed_GtUword2seed(seedpair_values2,encoding,bits_tab);
      gt_diagbandseed_seedpaircmp(seedpair_values,seedpair_values2);
    }
    gt_diagbandseed_encode_seedpair(bb_read,
                                    current,
                                    bytestring_length,
                                    seedpair_values,
                                    bits_tab,
                                    4);
    gt_diagbandseed_decode_seedpair(bb_write,
                                    seedpair_values,
                                    bits_tab,
                                    current,
                                    4);
    if (debug_seedpair && spidx > 0)
    {
      int ret = gt_diagbandseed_seeds_compare_bytestring(previous,
                                                         current,
                                                         bytestring_length);
      if (ret >= 0)
      {
        printf("ret=%d\n",ret);
        fprintf(stderr,"prev=");
        gt_diagbandseed_show_seed(stderr,seedpairlist,spidx-1);
        showbytestring(stderr,previous,bytestring_length);
        fprintf(stderr,"\ncurr=");
        gt_diagbandseed_show_seed(stderr,seedpairlist,spidx);
        showbytestring(stderr,bytestring,bytestring_length);
        fprintf(stderr,"\n");
        exit(EXIT_FAILURE);
      }
    }
    tmp = previous;
    previous = current;
    current = tmp;
  }
  gt_bitbuffer_delete(bb_read);
  gt_bitbuffer_delete(bb_write);
  gt_free(bytestring);
  gt_free(bytestring2);
}

/* Fill a GtDiagbandseedSeedPair list of equal kmers from the iterators. */
static void gt_diagbandseed_merge(GtSeedpairlist *seedpairlist,
                                  GtUword *histogram,
                                  bool knowthesize,
                                  GtDiagbandseedKmerIterator *aiter,
                                  const GtSequenceRangeWithMaxLength *aseqrange,
                                  GtDiagbandseedKmerIterator *biter,
                                  const GtSequenceRangeWithMaxLength *bseqrange,
                                  GtUword maxfreq,
                                  GtUword maxgram,
                                  const GtRange *seedpairdistance,
                                  bool selfcomp)
{
  const GtArrayGtDiagbandseedKmerPos *alist, *blist;

  gt_assert(aiter != NULL && biter != NULL &&
            ((histogram == NULL && seedpairlist != NULL) ||
            (histogram != NULL && seedpairlist == NULL)));
  alist = gt_diagbandseed_kmer_iter_next(aiter);
  blist = gt_diagbandseed_kmer_iter_next(biter);
  while (alist != NULL && blist != NULL) {
    const GtDiagbandseedKmerPos *asegment = alist->spaceGtDiagbandseedKmerPos,
                                *bsegment = blist->spaceGtDiagbandseedKmerPos;
    GtUword alen = alist->nextfreeGtDiagbandseedKmerPos,
            blen = blist->nextfreeGtDiagbandseedKmerPos;
    if (asegment->code < bsegment->code) {
      alist = gt_diagbandseed_kmer_iter_next(aiter);
    } else
    {
      if (asegment->code > bsegment->code)
      {
        blist = gt_diagbandseed_kmer_iter_next(biter);
      } else
      {
        GtUword frequency = MAX(alen, blen);
        if (frequency <= maxfreq)
        {
          /* add all equal k-mers */
          frequency = MIN(maxgram, frequency);
          gt_assert(frequency > 0);
          if (histogram != NULL && !selfcomp)
          {
            histogram[frequency - 1] += alen * blen;
          } else
          {
            const GtDiagbandseedKmerPos *aptr, *bptr;
            for (aptr = asegment; aptr < asegment + alen; aptr++)
            {
              for (bptr = bsegment; bptr < bsegment + blen; bptr++)
              {
                if (!selfcomp || aptr->seqnum < bptr->seqnum ||
                    (aptr->seqnum == bptr->seqnum &&
                     aptr->endpos + seedpairdistance->start <= bptr->endpos &&
                     aptr->endpos + seedpairdistance->end >= bptr->endpos))
                {
                  /* no duplicates from the same dataset */
                  if (histogram == NULL)
                  {
                    /* save SeedPair in seedpairlist */
                    gt_seedpairlist_add(seedpairlist,
                                        knowthesize,
                                        aseqrange,
                                        bseqrange,
                                        aptr->seqnum,
                                        bptr->seqnum,
                                        bptr->endpos,
                                        aptr->endpos);
                  } else
                  {
                    /* count seed pair frequency in histogram */
                    histogram[frequency - 1]++;
                  }
                }
              }
            }
          }
        } /* else: ignore all equal elements */
        alist = gt_diagbandseed_kmer_iter_next(aiter);
        blist = gt_diagbandseed_kmer_iter_next(biter);
      }
    }
  }
}

/* Verify seed pairs in the original sequences */
static int gt_diagbandseed_verify(const GtSeedpairlist *seedpairlist,
                                  const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  unsigned int seedlength,
                                  bool reverse,
                                  bool verbose,
                                  FILE *stream,
                                  GtError *err) {
  const GtUword mlistlen = gt_seedpairlist_length(seedpairlist);
  GtTimer *timer = gt_timer_new();
  GtUword idx;
  char *buf1 = gt_malloc(3 * (seedlength + 1) * sizeof *buf1);
  char *buf2 = buf1 + 1 + seedlength;
  char *buf3 = buf2 + 1 + seedlength;
  buf1[seedlength] = buf2[seedlength] = buf3[seedlength] = '\0';

  if (verbose) {
    fprintf(stream, "# Start verifying seed pairs...\n");
    gt_timer_start(timer);
  }

  gt_assert(aencseq != NULL && bencseq != NULL);
  for (idx = 0; idx < mlistlen; idx++)
  {
    GtDiagbandseedSeqnum aseqnum = gt_seedpairlist_aseqnum(seedpairlist,idx),
                         bseqnum = gt_seedpairlist_bseqnum(seedpairlist,idx);
    GtDiagbandseedPosition bpos = gt_seedpairlist_bpos(seedpairlist,idx),
                           apos = gt_seedpairlist_apos(seedpairlist,idx);
    GtDiagbandseedPosition abs_apos, abs_bpos;

    /* extract decoded k-mers at seed pair positions */
    abs_apos = apos + gt_encseq_seqstartpos(aencseq, aseqnum);
    gt_encseq_extract_decoded(aencseq, buf1, abs_apos + 1 - seedlength,
                              abs_apos);

    if (!reverse) {
      abs_bpos = bpos + gt_encseq_seqstartpos(bencseq, bseqnum);
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos + 1 - seedlength,
                                abs_bpos);
      if (strcmp(buf1, buf2) != 0) {
        gt_error_set(err, "Wrong SeedPair (" "%"PRIu32
                          ",%"PRIu32",%"PRIu32",%"PRIu32"): %s != %s\n",
                     aseqnum, bseqnum, apos, bpos, buf1, buf2);
        gt_free(buf1);
        gt_timer_delete(timer);
        return -1;
      }
    } else {
      /* get reverse k-mer */
      char *bufptr;
      abs_bpos = gt_encseq_seqstartpos(bencseq, bseqnum) +
                 gt_encseq_seqlength(bencseq, bseqnum) - bpos - 1;
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos,
                                abs_bpos + seedlength - 1);

      for (bufptr = buf3; bufptr < buf3 + seedlength; bufptr++) {
        gt_complement(bufptr, buf2[seedlength + buf3 - bufptr - 1], NULL);
      }
      if (strcmp(buf1, buf3) != 0) {
        gt_error_set(err, "Wrong SeedPair (" "%"PRIu32
                          ",%"PRIu32",%"PRIu32",%"PRIu32"): %s != %s\n",
                     aseqnum, bseqnum, apos, bpos, buf1, buf3);
        gt_free(buf1);
        gt_timer_delete(timer);
        return -1;
      }
    }
  }
  if (verbose) {
    fprintf(stream, "# ...successfully verified each seed pair ");
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
  }
  gt_free(buf1);
  gt_timer_delete(timer);
  return 0;
}

/* Return estimated length of mlist, and maxfreq w.r.t. the given memlimit */
static int gt_diagbandseed_get_mlistlen_maxfreq(GtUword *mlistlen,
                                            GtUword *maxfreq,
                                            GtDiagbandseedKmerIterator *aiter,
                                            GtDiagbandseedKmerIterator *biter,
                                            GtUword memlimit,
                                            GtDiagbandseedPairlisttype splt,
                                            const GtRange *seedpairdistance,
                                            GtUword len_used,
                                            bool selfcomp,
                                            bool alist_blist_id,
                                            bool verbose,
                                            FILE *stream,
                                            GtError *err)
{
  const GtUword maxgram = MIN(*maxfreq, 8190) + 1; /* Cap on k-mer count */
  GtUword *histogram = NULL;
  GtTimer *timer = NULL;
  int had_err = 0;

  if (memlimit == GT_UWORD_MAX) {
    return 0; /* no histogram calculation */
  }

  if (verbose) {
    timer = gt_timer_new();
    fprintf(stream, "# Start calculating k-mer frequency histogram...\n");
    gt_timer_start(timer);
  }

  /* build histogram; histogram[maxgram] := estimation for mlistlen */
  histogram = gt_calloc(maxgram + 1, sizeof *histogram);
  gt_diagbandseed_merge(NULL, /* mlist not needed: just count */
                        histogram,
                        false,
                        aiter,
                        NULL,
                        biter,
                        NULL,
                        *maxfreq,
                        maxgram,
                        seedpairdistance,
                        selfcomp);
  *maxfreq = gt_diagbandseed_processhistogram(histogram,
                                              *maxfreq,
                                              maxgram,
                                              memlimit,
                                              len_used *
                                                sizeof (GtDiagbandseedKmerPos),
                                              alist_blist_id,
                                              gt_seedpairlist_sizeofunit(splt));
  *mlistlen = histogram[maxgram];
  gt_free(histogram);

  if (verbose) {
    gt_timer_show_formatted(timer,
                            "# ...finished histogram " GT_DIAGBANDSEED_FMT,
                            stream);
    gt_timer_delete(timer);
  }

  /* check maxfreq value */
  if (*maxfreq == 0 || (*maxfreq == 1 && alist_blist_id)) {
    gt_error_set(err,
                 "option -memlimit too strict: need at least " GT_WU "MB",
                 (*mlistlen >> 20) + 1);
    *mlistlen = 0;
    had_err = -1;
  } else if (verbose) {
    if (*maxfreq == GT_UWORD_MAX) {
      fprintf(stream, "# Disable k-mer maximum frequency, ");
    } else {
      fprintf(stream, "# Set k-mer maximum frequency to " GT_WU ", ", *maxfreq);
    }
    fprintf(stream, "expect " GT_WU " seed pairs.\n", *mlistlen);
  } else if (*maxfreq <= 5) {
    gt_warning("Only k-mers occurring <= " GT_WU " times will be considered, "
               "due to small memlimit.", *maxfreq);
  }

  return had_err;
}

/* Return a sorted list of SeedPairs from given Kmer-Iterators.
 * Parameter known_size > 0 can be given to allocate the memory beforehand.
 * The caller is responsible for freeing the result. */
static void gt_diagbandseed_get_seedpairs(GtSeedpairlist *seedpairlist,
                                          GtDiagbandseedKmerIterator *aiter,
                                          GtDiagbandseedKmerIterator *biter,
                                          GtUword maxfreq,
                                          GtUword known_size,
                                          const GtRange *seedpairdistance,
                                          bool selfcomp,
                                          const GtSequenceRangeWithMaxLength
                                            *aseqrange,
                                          const GtSequenceRangeWithMaxLength
                                            *bseqrange,
                                          bool debug_seedpair,
                                          bool verbose,
                                          FILE *stream)
{
  GtTimer *timer = NULL;
  GtUword mlistlen;
  const GtUword anumofseq = aseqrange->end - aseqrange->start + 1,
                amaxlen = aseqrange->max_length,
                bnumofseq = bseqrange->end - bseqrange->start + 1,
                bmaxlen = bseqrange->max_length;
  GtBitcount_type bits_tab[4];

  bits_tab[idx_aseqnum] = (GtBitcount_type) gt_radixsort_bits(anumofseq);
  bits_tab[idx_apos] = (GtBitcount_type) gt_radixsort_bits(amaxlen);
  bits_tab[idx_bseqnum] = (GtBitcount_type) gt_radixsort_bits(bnumofseq);
  bits_tab[idx_bpos] = (GtBitcount_type) gt_radixsort_bits(bmaxlen);
  if (verbose) {
    timer = gt_timer_new();
    if (known_size > 0) {
      fprintf(stream, "# Start building " GT_WU " seed pairs...\n",
                      known_size);
    } else {
      fprintf(stream, "# Start building seed pairs...\n");
    }
    gt_timer_start(timer);
  }

  /* allocate mlist space according to seed pair count */
  gt_seedpairlist_init(seedpairlist,known_size);

  /* create mlist */
  gt_diagbandseed_merge(seedpairlist,
                        NULL, /* histogram not needed: save seed pairs */
                        known_size > 0 ? true : false,
                        aiter,
                        aseqrange,
                        biter,
                        bseqrange,
                        maxfreq,
                        GT_UWORD_MAX, /* maxgram */
                        seedpairdistance,
                        selfcomp);
  mlistlen = gt_seedpairlist_length(seedpairlist);
  if (verbose) {
    fprintf(stream, "# ...collected " GT_WU " seed pairs ", mlistlen);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
  }

  /* sort mlist */
  if (mlistlen > 0) {
    if (verbose) {
      gt_timer_start(timer);
    }
     gt_diagbandseed_seedpairlist_sort(seedpairlist);
    if (verbose) {
      fprintf(stream, "# ...sorted " GT_WU " seed pairs ", mlistlen);
      gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    }

    if (debug_seedpair) {
      gt_diagbandseed_seedpairlist_out(stream,seedpairlist);
    }
  }
  if (verbose) {
    gt_timer_start(timer);
  }
  if (seedpairlist == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    gt_diagbandseed_seedpairlist_encode_decode(seedpairlist,bits_tab,
                                               debug_seedpair);
  }
  if (verbose) {
    fprintf(stream, "# ...encoding/decoding of seedpairs seed pairs ");
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }
}

static int gt_diagbandseed_update_dband(GtDiagbandseedScore *diagband_score,
                                       GtDiagbandseedPosition *diagband_lastpos,
                                       GtUword diag,
                                       GtUword bpos,
                                       unsigned int seedlength)
{
  GtUword addlength;

  if (bpos >= diagband_lastpos[diag] + seedlength)
  {
    /* no overlap */
    addlength = (GtUword) seedlength;
  } else
  {
    /* overlap: add positons after last counted position */
    gt_assert(diagband_lastpos[diag] <= bpos);
    addlength = bpos - diagband_lastpos[diag];
  }
  diagband_lastpos[diag] = bpos;
  if (addlength > 0)
  {
    diagband_score[diag] += addlength;
    if (diagband_score[diag] == addlength)
    {
      return 1;
    }
  }
  return 0;
}

static int gt_diagbandseed_possibly_extend(const GtQuerymatch *previousmatch,
                                        GtUword aseqnum,
                                        GtUword apos,
                                        GtUword bseqnum,
                                        GtUword bpos,
                                        bool use_apos,
                                        unsigned int seedlength,
                                        GtUword errorpercentage,
                                        GtUword userdefinedleastlength,
                                        const GtEncseq *aencseq,
                                        const GtEncseq *bencseq,
                                        GtProcessinfo_and_querymatchspaceptr
                                          *info_querymatch,
                                        GtReadmode query_readmode,
                                        GtExtendSelfmatchRelativeFunc
                                          extend_selfmatch_relative_function,
                                        GtExtendQuerymatchRelativeFunc
                                          extend_querymatch_relative_function)
{
  int ret = 0;
  if (previousmatch == NULL ||
      !gt_querymatch_overlap(previousmatch,apos,bpos,use_apos))
  {
    /* extend seed */
    const GtQuerymatch *querymatch;
    /* relative seed start position in A and B */
    const GtUword bstart = (GtUword) (bpos + 1 - seedlength);
    const GtUword astart = (GtUword) (apos + 1 - seedlength);

    ret = 1; /* perform extension */
    if (aencseq == bencseq)
    {
      querymatch = extend_selfmatch_relative_function(info_querymatch,
                                                      aencseq,
                                                      aseqnum,
                                                      astart,
                                                      bseqnum,
                                                      bstart,
                                                      seedlength,
                                                      query_readmode);
    } else
    {
      querymatch = extend_querymatch_relative_function(info_querymatch,
                                                       aencseq,
                                                       aseqnum,
                                                       astart,
                                                       bencseq,
                                                       bseqnum,
                                                       bstart,
                                                       seedlength,
                                                       query_readmode);
    }
    if (querymatch != NULL)
    {
      /* show extension results */
      if (gt_querymatch_check_final(querymatch, errorpercentage,
                                    userdefinedleastlength))
      {
        gt_querymatch_prettyprint(querymatch);
      }
      ret = 2; /* output match */
    }
  }
  return ret;
}

/*#define GT_DIAGBANDSEED_SEEDHISTOGRAM 100*/
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
static void gt_diagbandseed_seedhistogram_out(FILE *stream,
                                              const GtUword *seedhistogram)
{
  int seedcount;
  GtUword sum = 0;

  fprintf(stream, "# seed histogram:");
  for (seedcount = 0; seedcount < GT_DIAGBANDSEED_SEEDHISTOGRAM; seedcount++)
  {
    if (seedcount % 10 == 0)
    {
      fprintf(stream, "\n#\t");
    }
    sum += seedcount * seedhistogram[seedcount];
    fprintf(stream, GT_WU "\t", seedhistogram[seedcount]);
  }
  fprintf(stream, "\n# sum = " GT_WU "\n",sum);
}
#endif

/* * * * * SEED EXTENSION * * * * */

#define GT_DIAGBANDSEED_DIAG(AMAXLEN,APOS,BPOS)\
        (((AMAXLEN) + (GtUword) (BPOS) - (GtUword) (APOS)) \
         >> arg->logdiagbandwidth)

/* start seed extension for seed pairs in mlist */
static void gt_diagbandseed_process_seeds(const GtSeedpairlist *seedpairlist,
                                          const GtDiagbandseedExtendParams *arg,
                                          void *processinfo,
                                          GtQuerymatchoutoptions *querymoutopt,
                                          const GtEncseq *aencseq,
                                          GT_UNUSED
                                          const GtSequenceRangeWithMaxLength
                                            *aseqrange,
                                          const GtEncseq *bencseq,
                                          GT_UNUSED
                                          const GtSequenceRangeWithMaxLength
                                            *bseqrange,
                                          unsigned int seedlength,
                                          bool reverse,
                                          bool verbose,
                                          FILE *stream)
{
  GtDiagbandseedScore *diagband_score;
  GtDiagbandseedPosition *diagband_lastpos;
  GtExtendSelfmatchRelativeFunc extend_selfmatch_relative_function = NULL;
  GtExtendQuerymatchRelativeFunc extend_querymatch_relative_function = NULL;
  GtProcessinfo_and_querymatchspaceptr info_querymatch;
  /* Although the sequences of the parts processed are shorter, we need to
     set amaxlen and bmaxlen to the maximum size of all sequences
     to get the same division into diagonal bands for all parts and thus
     obtain results independent of the number of parts chosen. */
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq),
                bmaxlen = gt_encseq_max_seq_length(bencseq),
                mlistlen = gt_seedpairlist_length(seedpairlist),
                ndiags = 1 + ((amaxlen + bmaxlen) >> arg->logdiagbandwidth),
                minsegmentlen = (arg->mincoverage - 1) / seedlength + 1;
  GtUword count_extensions = 0;
  GtReadmode query_readmode;
  GtTimer *timer = NULL;
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
  GtUword *seedhistogram = gt_calloc(GT_DIAGBANDSEED_SEEDHISTOGRAM,
                                     sizeof *seedhistogram);
  GtUword seedcount = 0;
#endif

  /* select extension method */
  if (mlistlen == 0 || mlistlen < minsegmentlen) {
    return;
  }
  if (arg->extendgreedy) {
    extend_selfmatch_relative_function = gt_greedy_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_greedy_extend_querymatch_relative;
  } else if (arg->extendxdrop) {
    extend_selfmatch_relative_function = gt_xdrop_extend_selfmatch_relative;
    extend_querymatch_relative_function = gt_xdrop_extend_querymatch_relative;
  } else { /* no seed extension */
    return;
  }

  if (verbose) {
    GtStr *add_column_header;
    timer = gt_timer_new();
    fprintf(stream, "# Start %s seed pair extension ...\n",
                     arg->extendgreedy ? "greedy" : "xdrop");
    fprintf(stream,"# Parameters for selecting seeds: diagonal bands=" GT_WU
                    ", minimal segmentsize=" GT_WU ", minimal coverage=" GT_WU
                    "\n",
                    ndiags,minsegmentlen,arg->mincoverage);
    fprintf(stream, "# Columns: alen aseq astartpos strand blen bseq bstartpos "
            "score editdist identity");
    add_column_header = gt_querymatch_column_header(arg->display_flag);
    if (gt_str_length(add_column_header) > 0)
    {
      fputs(gt_str_get(add_column_header),stream);
    }
    fputc('\n',stream);
    gt_str_delete(add_column_header);
    gt_timer_start(timer);
  }
  info_querymatch.processinfo = processinfo;
  info_querymatch.querymatchspaceptr = gt_querymatch_new();
  gt_querymatch_display_set(info_querymatch.querymatchspaceptr,
                            arg->display_flag);
  if (gt_querymatch_evalue_display(arg->display_flag) ||
      gt_querymatch_bit_score_display(arg->display_flag))
  {
    info_querymatch.karlin_altschul_stat = gt_karlin_altschul_stat_new_gapped();
    gt_karlin_altschul_stat_add_keyvalues(info_querymatch.karlin_altschul_stat,
                                          gt_encseq_total_length(aencseq),
                                          gt_encseq_num_of_sequences(aencseq));
  } else
  {
    info_querymatch.karlin_altschul_stat = NULL;
  }
  if (arg->verify_alignment)
  {
    gt_querymatch_verify_alignment_set(info_querymatch.querymatchspaceptr);
  }
  if (querymoutopt != NULL) {
    gt_querymatch_outoptions_set(info_querymatch.querymatchspaceptr,
                                 querymoutopt);
  }
  query_readmode = reverse ? GT_READMODE_REVCOMPL
                           : GT_READMODE_FORWARD;
  gt_querymatch_query_readmode_set(info_querymatch.querymatchspaceptr,
                                   query_readmode);
  gt_querymatch_file_set(info_querymatch.querymatchspaceptr, stream);

  /* diagband_score[0] and diagband_score[ndiags+1] remain zero as boundaries */
  diagband_score = gt_calloc(ndiags + 2, sizeof *diagband_score);
  diagband_score++; /* so we need not increment the index */
  diagband_lastpos = gt_calloc(ndiags, sizeof *diagband_lastpos);
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    const GtDiagbandseedSeedPair
      *mlist = gt_seedpairlist_mlist_struct(seedpairlist),
      *mlistend = mlist + mlistlen,
      *maxsegm_start = mlistend - minsegmentlen,
      *nextsegm = mlist;

    /* iterate through segments of equal k-mers */
    while (nextsegm <= maxsegm_start)
    {
      const GtDiagbandseedSeedPair *currsegm = nextsegm, *seedpair;
      const GtDiagbandseedSeqnum currsegm_aseqnum = currsegm->aseqnum;
      const GtDiagbandseedSeqnum currsegm_bseqnum = currsegm->bseqnum;
      bool haspreviousmatch;
      GtUword diagbands_used;

      /* if insuffienct number of kmers in segment: skip whole segment */
      if (currsegm_aseqnum != currsegm[minsegmentlen - 1].aseqnum ||
          currsegm_bseqnum != currsegm[minsegmentlen - 1].bseqnum)
      {
        do
        {
          nextsegm++;
        } while (nextsegm < mlistend &&
                 currsegm_aseqnum == nextsegm->aseqnum &&
                 currsegm_bseqnum == nextsegm->bseqnum);
        continue; /* process next segment */
      }

      /* this segment begining with nextsegm pssibly has enough seeds */
      diagbands_used = 0;
      do
      {
        const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,nextsegm->apos,
                                                  nextsegm->bpos);

        gt_assert(diag < ndiags);
        diagbands_used += gt_diagbandseed_update_dband(diagband_score,
                                                       diagband_lastpos,
                                                       diag,
                                                       nextsegm->bpos,
                                                       seedlength);
        nextsegm++;
      } while (nextsegm < mlistend &&
               nextsegm->aseqnum == currsegm_aseqnum &&
               nextsegm->bseqnum == currsegm_bseqnum);

      /* from here on we only need the apos and bpos values of the segment, as
         the segment boundaries have been identified.
         second scan: test for mincoverage and overlap to previous extension,
         based on apos and bpos values. */
      for (seedpair = currsegm, haspreviousmatch = false; seedpair < nextsegm;
           seedpair++)
      {
        const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,seedpair->apos,
                                                  seedpair->bpos);
        if ((GtUword) MAX(diagband_score[diag + 1], diagband_score[diag - 1]) +
            (GtUword) diagband_score[diag]
            >= arg->mincoverage)
        {
          int ret = gt_diagbandseed_possibly_extend(
                          haspreviousmatch ? info_querymatch.querymatchspaceptr
                                           : NULL,
                          currsegm_aseqnum,
                          seedpair->apos,
                          currsegm_bseqnum,
                          seedpair->bpos,
                          arg->use_apos,
                          seedlength,
                          arg->errorpercentage,
                          arg->userdefinedleastlength,
                          aencseq,
                          bencseq,
                          &info_querymatch,
                          query_readmode,
                          extend_selfmatch_relative_function,
                          extend_querymatch_relative_function);
          if (ret >= 1)
          {
            count_extensions++;
          }
          if (ret == 2)
          {
            haspreviousmatch = true;
          }
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
          seedcount++;
#endif
        }
      }

      /* reset diagonal band scores */
      if (diagbands_used * 3 >= ndiags) /* >= 33% of diagbands are used */
      {
        memset(diagband_score,0,sizeof *diagband_score * ndiags);
        memset(diagband_lastpos,0,sizeof *diagband_lastpos * ndiags);
      } else
      {
        /* third scan: */
        for (seedpair = currsegm; seedpair < nextsegm; seedpair++)
        {
          const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,seedpair->apos,
                                                    seedpair->bpos);
          diagband_score[diag] = 0;
          diagband_lastpos[diag] = 0;
        }
      }
#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
      seedhistogram[MIN(GT_DIAGBANDSEED_SEEDHISTOGRAM - 1, seedcount)]++;
      seedcount = 0;
#endif
    }
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      const GtUword *mlist = gt_seedpairlist_mlist_ulong(seedpairlist),
                    *mlistend = mlist + mlistlen,
                    *maxsegm_start = mlistend - minsegmentlen,
                    *nextsegm = mlist;

      /* iterate through segments of equal k-mers */
      while (nextsegm <= maxsegm_start)
      {
        GtUword currsegm_aseqnum, currsegm_bseqnum;
        const GtUword *currsegm = nextsegm, *seedpair;
        const GtUword currsegm_a_bseqnum
          = gt_seedpairlist_a_bseqnum (seedpairlist,*currsegm);
        bool haspreviousmatch;
        GtUword diagbands_used;

        /* if insuffienct number of kmers in segment: skip whole segment */
        if (currsegm_a_bseqnum !=
            gt_seedpairlist_a_bseqnum (seedpairlist,currsegm[minsegmentlen-1]))
        {
          do
          {
            nextsegm++;
          } while (nextsegm < mlistend &&
                   currsegm_a_bseqnum ==
                   gt_seedpairlist_a_bseqnum (seedpairlist,*nextsegm));
          continue; /* process next segment */
        }

        /* this segment begining with nextsegm pssibly has enough seeds */
        diagbands_used = 0;
        do
        {
          const GtDiagbandseedPosition apos
            = gt_seedpairlist_extract(seedpairlist,*nextsegm,idx_apos),
          bpos = gt_seedpairlist_extract(seedpairlist,*nextsegm,idx_bpos);
          const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen, apos, bpos);

          gt_assert(diag < ndiags);
          diagbands_used += gt_diagbandseed_update_dband(diagband_score,
                                                         diagband_lastpos,
                                                         diag,
                                                         bpos,
                                                         seedlength);
          nextsegm++;
        } while (nextsegm < mlistend &&
                 currsegm_a_bseqnum ==
                 gt_seedpairlist_a_bseqnum (seedpairlist,*nextsegm));

        /* from here on we only need the apos and bpos values of the segment, as
           the segment boundaries have been identified.
           second scan: test for mincoverage and overlap to previous extension,
           based on apos and bpos values. */
        currsegm_aseqnum = gt_seedpairlist_extract(seedpairlist,*currsegm,
                                                   idx_aseqnum) +
                           aseqrange->start;
        currsegm_bseqnum = gt_seedpairlist_extract(seedpairlist,*currsegm,
                                                   idx_bseqnum) +
                           bseqrange->start;
        for (seedpair = currsegm, haspreviousmatch = false; seedpair < nextsegm;
             seedpair++)
        {
          const GtDiagbandseedPosition apos
            = gt_seedpairlist_extract(seedpairlist,*seedpair,idx_apos),
          bpos = gt_seedpairlist_extract(seedpairlist,*seedpair,idx_bpos);
          const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,apos,bpos);
          if ((GtUword) MAX(diagband_score[diag + 1], diagband_score[diag - 1])
              + (GtUword) diagband_score[diag]
              >= arg->mincoverage)
          {
            int ret = gt_diagbandseed_possibly_extend(
                           haspreviousmatch ? info_querymatch.querymatchspaceptr
                                            : NULL,
                           currsegm_aseqnum,
                           apos,
                           currsegm_bseqnum,
                           bpos,
                           arg->use_apos,
                           seedlength,
                           arg->errorpercentage,
                           arg->userdefinedleastlength,
                           aencseq,
                           bencseq,
                           &info_querymatch,
                           query_readmode,
                           extend_selfmatch_relative_function,
                           extend_querymatch_relative_function);
            if (ret >= 1)
            {
              count_extensions++;
            }
            if (ret == 2)
            {
              haspreviousmatch = true;
            }
  #ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
            seedcount++;
  #endif
          }
        }

        /* reset diagonal band scores */
        if (diagbands_used * 3 >= ndiags) /* >= 33% of diagbands are used */
        {
          memset(diagband_score,0,sizeof *diagband_score * ndiags);
          memset(diagband_lastpos,0,sizeof *diagband_lastpos * ndiags);
        } else
        {
          /* third scan: */
          for (seedpair = currsegm; seedpair < nextsegm; seedpair++)
          {
            const GtDiagbandseedPosition apos
              = gt_seedpairlist_extract(seedpairlist,*seedpair,idx_apos),
            bpos = gt_seedpairlist_extract(seedpairlist,*seedpair,idx_bpos);
            const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,apos,bpos);
            diagband_score[diag] = 0;
            diagband_lastpos[diag] = 0;
          }
        }
  #ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
        seedhistogram[MIN(GT_DIAGBANDSEED_SEEDHISTOGRAM - 1, seedcount)]++;
        seedcount = 0;
  #endif
      }
    } else
    {
      gt_assert(false);
    }
  }
  diagband_score--; /* need to recover original base adress */
  gt_free(diagband_score);
  gt_free(diagband_lastpos);
  gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  gt_karlin_altschul_stat_delete(info_querymatch.karlin_altschul_stat);

#ifdef GT_DIAGBANDSEED_SEEDHISTOGRAM
  gt_diagbandseed_seedhistogram_out(stream,seedhistogram);
  gt_free(seedhistogram);
#endif
  if (verbose)
  {
    fprintf(stream, "# ...finished " GT_WU " seed pair extension%s ",
            count_extensions, count_extensions > 1 ? "s" : "");
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }
}

/* * * * * ALGORITHM STEPS * * * * */

static char *gt_diagbandseed_kmer_filename(const GtEncseq *encseq,
                                           unsigned int seedlength,
                                           bool forward,
                                           unsigned int numparts,
                                           unsigned int partindex)
{
  char *filename;
  GtStr *str = gt_str_new_cstr(gt_encseq_indexname(encseq));
  gt_str_append_char(str, '.');
  gt_str_append_uint(str, seedlength);
  gt_str_append_char(str, forward ? 'f' : 'r');
  gt_str_append_uint(str, numparts);
  gt_str_append_char(str, '-');
  gt_str_append_uint(str, partindex + 1);
  gt_str_append_cstr(str, ".kmer");
  filename = gt_cstr_dup(gt_str_get(str));
  gt_str_delete(str);
  return filename;
}

/* Go through the different steps of the seed and extend algorithm. */
static int gt_diagbandseed_algorithm(const GtDiagbandseedInfo *arg,
                                     const GtArrayGtDiagbandseedKmerPos *alist,
                                     FILE *stream,
                                     const GtSequenceRangeWithMaxLength
                                       *aseqrange,
                                     GtUword partindex_a,
                                     const GtSequenceRangeWithMaxLength
                                       *bseqrange,
                                     GtUword partindex_b,
                                     GtFtTrimstat *trimstat,
                                     GtError *err)
{
  GtArrayGtDiagbandseedKmerPos blist;
  GtSeedpairlist *seedpairlist = NULL, *revseedpairlist = NULL;;
  GtDiagbandseedKmerIterator *aiter = NULL, *biter = NULL;
  GtUword alen = 0, blen = 0, mlistlen = 0, maxfreq, len_used;
  GtRange seedpairdistance = *arg->seedpairdistance;
  char *blist_file = NULL;
  int had_err = 0;
  bool alist_blist_id, both_strands, selfcomp, equalranges, use_blist = false;
  const GtDiagbandseedExtendParams *extp = NULL;
  GtFtPolishing_info *pol_info = NULL;
  void *processinfo = NULL;
  GtQuerymatchoutoptions *querymoutopt = NULL;

  gt_assert(arg != NULL && aseqrange != NULL && bseqrange != NULL);
  maxfreq = arg->maxfreq;
  selfcomp = (arg->bencseq == arg->aencseq &&
              aseqrange->start <= bseqrange->end &&
              aseqrange->end >= bseqrange->start)
              ? true : false;
  equalranges = (aseqrange->start == bseqrange->start &&
                 aseqrange->end == bseqrange->end) ? true : false;
  alist_blist_id = (selfcomp && !arg->nofwd && equalranges) ? true : false;
  both_strands = (arg->norev || arg->nofwd) ? false : true;
  if (!alist_blist_id) {
    seedpairdistance.start = 0UL;
  }

  if (arg->verbose && (arg->anumseqranges > 1 || arg->bnumseqranges > 1)) {
    fprintf(stream, "# Process part " GT_WU " (sequences " GT_WU "..." GT_WU
                    ") vs part " GT_WU " (sequences " GT_WU "..." GT_WU ")\n",
            partindex_a + 1, aseqrange->start,aseqrange->end,
            partindex_b + 1, bseqrange->start,bseqrange->end);
  }

  /* Create k-mer iterator for alist */
  if (alist == NULL) {
    char *alist_file;
    alist_file = gt_diagbandseed_kmer_filename(arg->aencseq, arg->seedlength,
                                               true, arg->anumseqranges,
                                               partindex_a);
    FILE *alist_fp = gt_fa_fopen(alist_file, "rb", err);
    if (alist_fp == NULL) {
      return -1;
    }
    alen = (GtUword)(gt_file_size(alist_file) / sizeof (GtDiagbandseedKmerPos));
    aiter = gt_diagbandseed_kmer_iter_new_file(alist_fp);
    gt_free(alist_file);
    alist_file = NULL;
  } else {
    gt_assert(alist != NULL);
    alen = alist->nextfreeGtDiagbandseedKmerPos;
    aiter = gt_diagbandseed_kmer_iter_new_list(alist);
  }

  /* Second k-mer list */
  if (alist_blist_id && alist != NULL) {
    biter = gt_diagbandseed_kmer_iter_new_list(alist);
    blen = alen;
  } else if (arg->use_kmerfile) {
    blist_file = gt_diagbandseed_kmer_filename(arg->bencseq, arg->seedlength,
                                               !arg->nofwd, arg->bnumseqranges,
                                               partindex_b);
    if (!gt_file_exists(blist_file)) {
      gt_free(blist_file);
      blist_file = NULL;
    }
  }
  if (blist_file != NULL) {
    FILE *blist_fp = gt_fa_fopen(blist_file, "rb", err);
    if (blist_fp == NULL) {
      gt_free(blist_file);
      gt_diagbandseed_kmer_iter_delete(aiter);
      return -1;
    }
    blen = (GtUword)(gt_file_size(blist_file) / sizeof (GtDiagbandseedKmerPos));
    gt_assert(biter == NULL);
    biter = gt_diagbandseed_kmer_iter_new_file(blist_fp);
    gt_free(blist_file);
    blist_file = NULL;
  } else if (!alist_blist_id) {
    const GtReadmode readmode = arg->nofwd ? GT_READMODE_COMPL
                                           : GT_READMODE_FORWARD;
    const GtUword known_size = (selfcomp && equalranges) ? alen : 0;
    blist = gt_diagbandseed_get_kmers(arg->bencseq,
                                      arg->seedlength,
                                      readmode,
                                      bseqrange,
                                      arg->debug_kmer,
                                      arg->verbose,
                                      known_size,
                                      stream);
    blen = blist.nextfreeGtDiagbandseedKmerPos;
    biter = gt_diagbandseed_kmer_iter_new_list(&blist);
    use_blist = true;
  }

  len_used = alen;
  if (!selfcomp || !arg->norev) {
    len_used += blen;
  }
  had_err = gt_diagbandseed_get_mlistlen_maxfreq(&mlistlen,
                                                 &maxfreq,
                                                 aiter,
                                                 biter,
                                                 arg->memlimit,
                                                 arg->splt,
                                                 &seedpairdistance,
                                                 len_used,
                                                 selfcomp,
                                                 alist_blist_id,
                                                 arg->verbose,
                                                 stream,
                                                 err);

  if (!had_err) {
    gt_diagbandseed_kmer_iter_reset(aiter);
    gt_diagbandseed_kmer_iter_reset(biter);
    seedpairlist = gt_seedpairlist_new(arg->splt,aseqrange,bseqrange);
    if (arg->verbose)
    {
      gt_seedpairlist_show_bits(stream,seedpairlist);
    }
    gt_diagbandseed_get_seedpairs(seedpairlist,
                                  aiter,
                                  biter,
                                  maxfreq,
                                  mlistlen,
                                  &seedpairdistance,
                                  selfcomp,
                                  aseqrange,
                                  bseqrange,
                                  arg->debug_seedpair,
                                  arg->verbose,
                                  stream);
    mlistlen = gt_seedpairlist_length(seedpairlist);
    if (arg->verify && mlistlen > 0) {
      had_err = gt_diagbandseed_verify(seedpairlist,
                                       arg->aencseq,
                                       arg->bencseq,
                                       arg->seedlength,
                                       arg->nofwd,
                                       arg->verbose,
                                       stream,
                                       err);
      if (had_err) {
        gt_seedpairlist_delete(seedpairlist);
        seedpairlist = NULL;
      }
    }
  }

  if (use_blist) {
    GT_FREEARRAY(&blist, GtDiagbandseedKmerPos);
  }
  use_blist = false;
  gt_diagbandseed_kmer_iter_delete(biter);
  if (had_err) {
    gt_diagbandseed_kmer_iter_delete(aiter);
    return had_err;
  }

  /* Create extension info objects */
  extp = arg->extp;
  if (extp->extendgreedy) {
    GtGreedyextendmatchinfo *grextinfo = NULL;
    const double weak_errorperc = (double)(extp->weakends
                                           ? MAX(extp->errorpercentage, 20)
                                           : extp->errorpercentage);

    pol_info = polishing_info_new_with_bias(weak_errorperc,
                                            extp->matchscore_bias,
                                            extp->history_size);
    grextinfo = gt_greedy_extend_matchinfo_new(extp->errorpercentage,
                                               extp->maxalignedlendifference,
                                               extp->history_size,
                                               extp->perc_mat_history,
                                               extp->userdefinedleastlength,
                                               extp->extend_char_access,
                                               extp->sensitivity,
                                               pol_info);
    if (extp->benchmark) {
      gt_greedy_extend_matchinfo_silent_set(grextinfo);
    }
    if (trimstat != NULL)
    {
      gt_greedy_extend_matchinfo_trimstat_set(grextinfo,trimstat);
    }
    processinfo = (void *)grextinfo;
  } else if (extp->extendxdrop) {
    GtXdropmatchinfo *xdropinfo = NULL;
    gt_assert(extp->extendgreedy == false);
    xdropinfo = gt_xdrop_matchinfo_new(extp->userdefinedleastlength,
                                       extp->errorpercentage,
                                       extp->xdropbelowscore,
                                       extp->sensitivity);
    if (extp->benchmark) {
      gt_xdrop_matchinfo_silent_set(xdropinfo);
    }
    processinfo = (void *)xdropinfo;
  }
  if (extp->extendxdrop || extp->alignmentwidth > 0 ||
      extp->verify_alignment) {
    querymoutopt = gt_querymatchoutoptions_new(true,
                                               false,
                                               extp->alignmentwidth);
    if (extp->extendxdrop || extp->extendgreedy) {
      const GtUword sensitivity = extp->extendxdrop ? 100UL : extp->sensitivity;
      gt_querymatchoutoptions_extend(querymoutopt,
                                     extp->errorpercentage,
                                     extp->maxalignedlendifference,
                                     extp->history_size,
                                     extp->perc_mat_history,
                                     extp->extend_char_access,
                                     extp->weakends,
                                     sensitivity,
                                     extp->matchscore_bias,
                                     extp->always_polished_ends,
                                     extp->display_flag);
    }
  }

  /* process first mlist */
  gt_diagbandseed_process_seeds(seedpairlist,
                                arg->extp,
                                processinfo,
                                querymoutopt,
                                arg->aencseq,
                                aseqrange,
                                arg->bencseq,
                                bseqrange,
                                arg->seedlength,
                                arg->nofwd,
                                arg->verbose,
                                stream);
  gt_seedpairlist_delete(seedpairlist);

  /* Third (reverse) k-mer list */
  if (both_strands) {
    GtUword mrevlen = 0;
    GtArrayGtDiagbandseedKmerPos clist;

    gt_assert(blist_file == NULL && !use_blist);
    seedpairdistance.start = 0UL;
    if (arg->use_kmerfile) {
      blist_file = gt_diagbandseed_kmer_filename(arg->bencseq, arg->seedlength,
                                                 false, arg->bnumseqranges,
                                                 partindex_b);
      if (!gt_file_exists(blist_file)) {
        gt_free(blist_file);
        blist_file = NULL;
      }
    }
    if (blist_file != NULL) {
      FILE *blist_fp = gt_fa_fopen(blist_file, "rb", err);
      if (blist_fp == NULL) {
        had_err = -1;
      } else {
        biter = gt_diagbandseed_kmer_iter_new_file(blist_fp);
      }
      gt_free(blist_file);
    } else {
      clist = gt_diagbandseed_get_kmers(arg->bencseq,
                                        arg->seedlength,
                                        GT_READMODE_COMPL,
                                        bseqrange,
                                        arg->debug_kmer,
                                        arg->verbose,
                                        blen,
                                        stream);
      biter = gt_diagbandseed_kmer_iter_new_list(&clist);
      use_blist = true;
    }

    if (!had_err) {
      gt_diagbandseed_kmer_iter_reset(aiter);
      had_err = gt_diagbandseed_get_mlistlen_maxfreq(&mrevlen,
                                                     &maxfreq,
                                                     aiter,
                                                     biter,
                                                     arg->memlimit,
                                                     arg->splt,
                                                     &seedpairdistance,
                                                     len_used,
                                                     selfcomp,
                                                     alist_blist_id,
                                                     arg->verbose,
                                                     stream,
                                                     err);
    }

    if (!had_err) {
      gt_diagbandseed_kmer_iter_reset(aiter);
      gt_diagbandseed_kmer_iter_reset(biter);
      revseedpairlist = gt_seedpairlist_new(arg->splt,aseqrange,bseqrange);
      gt_diagbandseed_get_seedpairs(revseedpairlist,
                                    aiter,
                                    biter,
                                    maxfreq,
                                    mrevlen,
                                    &seedpairdistance,
                                    selfcomp,
                                    aseqrange,
                                    bseqrange,
                                    arg->debug_seedpair,
                                    arg->verbose,
                                    stream);
      mrevlen = gt_seedpairlist_length(revseedpairlist);
      if (arg->verify && mrevlen > 0) {
        had_err = gt_diagbandseed_verify(revseedpairlist,
                                         arg->aencseq,
                                         arg->bencseq,
                                         arg->seedlength,
                                         true,
                                         arg->verbose,
                                         stream,
                                         err);
        if (had_err) {
          gt_seedpairlist_delete(revseedpairlist);
        }
      }
    }
    if (use_blist) {
      GT_FREEARRAY(&clist, GtDiagbandseedKmerPos);
    }
    gt_diagbandseed_kmer_iter_delete(biter);
  }
  gt_diagbandseed_kmer_iter_delete(aiter);

  /* Process second (reverse) mlist */
  if (!had_err && both_strands) {
    gt_diagbandseed_process_seeds(revseedpairlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  arg->aencseq,
                                  aseqrange,
                                  arg->bencseq,
                                  bseqrange,
                                  arg->seedlength,
                                  true,
                                  arg->verbose,
                                  stream);
  }
  gt_seedpairlist_delete(revseedpairlist);

  /* Clean up */
  if (extp->extendgreedy) {
    polishing_info_delete(pol_info);
    gt_greedy_extend_matchinfo_delete((GtGreedyextendmatchinfo *)processinfo);
  } else if (extp->extendxdrop) {
    gt_xdrop_matchinfo_delete((GtXdropmatchinfo *)processinfo);
  }
  gt_querymatchoutoptions_delete(querymoutopt);
  return had_err;
}

#ifdef GT_THREADS_ENABLED
typedef struct{
  const GtDiagbandseedInfo *arg;
  const GtArrayGtDiagbandseedKmerPos *alist;
  FILE *stream;
  const GtSequenceRangeWithMaxLength *aseqranges;
  const GtSequenceRangeWithMaxLength *bseqranges;
  GtArray *combinations;
  int had_err;
  GtError *err;
}GtDiagbandseedThreadInfo;

static
void gt_diagbandseed_thread_info_set(GtDiagbandseedThreadInfo *ti,
                                     const GtDiagbandseedInfo *arg,
                                     const GtArrayGtDiagbandseedKmerPos *alist,
                                     FILE *stream,
                                     const GtSequenceRangeWithMaxLength
                                       *aseqranges,
                                     const GtSequenceRangeWithMaxLength
                                       *bseqranges,
                                     GtArray *combinations,
                                     GtError *err)
{
  gt_assert(ti != NULL);
  ti->arg = arg;
  ti->alist = alist;
  ti->stream = stream;
  ti->aseqranges = aseqranges;
  ti->bseqranges = bseqranges;
  ti->combinations = gt_array_clone(combinations);
  ti->had_err = 0;
  ti->err = err;
}

static void *gt_diagbandseed_thread_algorithm(void *thread_info)
{
  GtDiagbandseedThreadInfo *info = (GtDiagbandseedThreadInfo *)thread_info;
  if (gt_array_size(info->combinations) > 0) {
    const GtUwordPair *last = gt_array_get_last(info->combinations),
                      *comb;

    for (comb = gt_array_get_first(info->combinations); comb <= last; comb++) {
      info->had_err = gt_diagbandseed_algorithm(info->arg,
                                                info->alist,
                                                info->stream,
                                                info->aseqranges + comb->a,
                                                comb->a,
                                                info->bseqranges + comb->b,
                                                comb->b,
                                                NULL,
                                                info->err);
      if (info->had_err) break;
    }
  }
  gt_array_delete(info->combinations);
  return NULL;
}
#endif

static int gt_diagbandseed_write_kmers(const GtArrayGtDiagbandseedKmerPos *list,
                                       const char *path,
                                       unsigned int seedlength,
                                       bool verbose,
                                       GtError *err)
{
  FILE *stream;

  if (verbose) {
    printf("# Write " GT_WU " %u-mers to file %s\n",
           list->nextfreeGtDiagbandseedKmerPos, seedlength, path);
  }

  stream = gt_fa_fopen(path, "wb", err);
  if (stream != NULL) {
    gt_xfwrite(list->spaceGtDiagbandseedKmerPos,
               sizeof *list->spaceGtDiagbandseedKmerPos,
               list->nextfreeGtDiagbandseedKmerPos, stream);
    gt_fa_fclose(stream);
    return 0;
  } else {
    return -1;
  }
}

/* Run the algorithm by iterating over all combinations of sequence ranges. */
int gt_diagbandseed_run(const GtDiagbandseedInfo *arg,
                        const GtSequenceRangeWithMaxLength *aseqranges,
                        const GtSequenceRangeWithMaxLength *bseqranges,
                        const GtUwordPair *pick,
                        GtError *err)
{
  const bool self = arg->aencseq == arg->bencseq ? true : false;
  const bool apick = pick->a != GT_UWORD_MAX ? true : false;
  const bool bpick = pick->b != GT_UWORD_MAX ? true : false;
  GtArrayGtDiagbandseedKmerPos alist;
  GtUword aidx, bidx;
  int had_err = 0;
#ifdef GT_THREADS_ENABLED
  GtDiagbandseedThreadInfo *tinfo = gt_malloc(gt_jobs * sizeof *tinfo);
  FILE **stream;
  unsigned int tidx;
  GtFtTrimstat *trimstat = NULL;

  /* create output streams */
  stream = gt_malloc(gt_jobs * sizeof *stream);
  stream[0] = stdout;
  for (tidx = 1; !had_err && tidx < gt_jobs; tidx++) {
    stream[tidx] = gt_xtmpfp_generic(NULL, TMPFP_OPENBINARY | TMPFP_AUTOREMOVE);
  }
#endif

  /* create all missing k-mer lists for bencseq */
  if (arg->use_kmerfile) {
    unsigned int count;
    for (count = 0; count < 2; count++) {
      const bool fwd = count == 0 ? true : false;
      if (fwd && (self || arg->nofwd)) continue;
      if (!fwd && arg->norev) continue;

      for (bidx = 0; !had_err && bidx < arg->bnumseqranges; bidx++) {
        char *path;
        if (bpick && pick->b != bidx) continue;

        path = gt_diagbandseed_kmer_filename(arg->bencseq, arg->seedlength, fwd,
                                             arg->bnumseqranges, bidx);
        if (!gt_file_exists(path)) {
          GtArrayGtDiagbandseedKmerPos blist;
          GtReadmode readmode = fwd ? GT_READMODE_FORWARD : GT_READMODE_COMPL;

          blist = gt_diagbandseed_get_kmers(arg->bencseq, arg->seedlength,
                                            readmode, bseqranges + bidx,
                                            arg->debug_kmer, arg->verbose, 0,
                                            stdout);
          had_err = gt_diagbandseed_write_kmers(&blist, path, arg->seedlength,
                                                arg->verbose, err);
          GT_FREEARRAY(&blist, GtDiagbandseedKmerPos);
        }
        gt_free(path);
      }
    }
  }
  if (!had_err && arg->trimstat_on)
  {
    trimstat = gt_ft_trimstat_new();
  }
  for (aidx = 0; !had_err && aidx < arg->anumseqranges; aidx++) {
    /* create alist here to prevent redundant calculations */
    char *path = NULL;
    bool use_alist = false;
    if (apick && pick->a != aidx) continue;

    if (arg->use_kmerfile) {
      path = gt_diagbandseed_kmer_filename(arg->aencseq, arg->seedlength, true,
                                           arg->anumseqranges, aidx);
    }

    if (!arg->use_kmerfile || !gt_file_exists(path)) {
      use_alist = true;
      alist = gt_diagbandseed_get_kmers(arg->aencseq,
                                        arg->seedlength,
                                        GT_READMODE_FORWARD,
                                        aseqranges + aidx,
                                        arg->debug_kmer,
                                        arg->verbose,
                                        0,
                                        stdout);
      if (arg->use_kmerfile) {
        had_err = gt_diagbandseed_write_kmers(&alist, path, arg->seedlength,
                                              arg->verbose, err);
      }
    }
    if (arg->use_kmerfile) {
      gt_free(path);
    }
    bidx = self ? aidx : 0;

#ifdef GT_THREADS_ENABLED
    if (gt_jobs <= 1) {
#endif
      while (!had_err && bidx < arg->bnumseqranges) {
        if (!bpick || pick->b == bidx) {
          /* start algorithm with chosen sequence ranges */
          had_err = gt_diagbandseed_algorithm(arg,
                                              use_alist ? &alist : NULL,
                                              stdout,
                                              aseqranges + aidx,
                                              aidx,
                                              bseqranges + bidx,
                                              bidx,
                                              trimstat,
                                              err);
        }
        bidx++;
      }
#ifdef GT_THREADS_ENABLED
    } else if (!arg->use_kmerfile) {
      const GtUword num_runs = bpick ? 1 : arg->bnumseqranges - bidx;
      const GtUword num_runs_per_thread = (num_runs - 1) / gt_jobs + 1;
      const GtUword num_threads = (num_runs - 1) / num_runs_per_thread + 1;
      GtArray *combinations = gt_array_new(sizeof (GtUwordPair));
      GtArray *threads = gt_array_new(sizeof (GtThread *));

      gt_assert(bidx < arg->bnumseqranges);
      gt_assert(num_threads <= gt_jobs);
      gt_assert(!bpick || num_threads == 1);

      /* start additional threads */
      for (tidx = 1; !had_err && tidx < num_threads; tidx++) {
        GtThread *thread;
        GtUword idx;
        bidx += num_runs_per_thread;
        const GtUword end = MIN(bidx + num_runs_per_thread, arg->bnumseqranges);

        for (idx = bidx; idx < end; idx++) {
          GtUwordPair comb = {aidx, idx};
          gt_array_add(combinations, comb);
        }
        gt_diagbandseed_thread_info_set(tinfo + tidx,
                                        arg,
                                        use_alist ? &alist : NULL,
                                        stream[tidx],
                                        aseqranges,
                                        bseqranges,
                                        combinations,
                                        err);
        gt_array_reset(combinations);
        if ((thread = gt_thread_new(gt_diagbandseed_thread_algorithm,
                                    tinfo + tidx, err)) != NULL) {
          gt_array_add(threads, thread);
        } else {
          had_err = -1;
        }
      }

      /* start main thread */
      if (!had_err) {
        GtUword idx;
        bidx = self ? aidx : 0;
        for (idx = bidx;
             idx < MIN(bidx + num_runs_per_thread, arg->bnumseqranges);
             ++idx) {
          if (!bpick || pick->b == idx) {
            GtUwordPair comb = {aidx, idx};
            gt_array_add(combinations, comb);
          }
        }
        gt_diagbandseed_thread_info_set(tinfo,
                                        arg,
                                        use_alist ? &alist : NULL,
                                        stream[0],
                                        aseqranges,
                                        bseqranges,
                                        combinations,
                                        err);
        gt_diagbandseed_thread_algorithm(tinfo);
      }

      /* clean up */
      for (tidx = 0; tidx < gt_array_size(threads); tidx++) {
        GtThread *thread = *(GtThread**) gt_array_get(threads, tidx);
        if (!had_err) {
          gt_thread_join(thread);
        }
        gt_thread_delete(thread);
      }
      gt_array_delete(threads);
      for (tidx = 0; tidx < num_threads && !had_err; tidx++) {
        had_err = tinfo[tidx].had_err;
      }
      gt_array_delete(combinations);
    }
#endif
    if (use_alist) {
      GT_FREEARRAY(&alist, GtDiagbandseedKmerPos);
    }
  }
#ifdef GT_THREADS_ENABLED
  if (gt_jobs > 1 && arg->use_kmerfile) {
    GtArray *combinations[gt_jobs];
    GtArray *threads = gt_array_new(sizeof (GtThread *));
    GtUword counter = 0;
    for (tidx = 0; tidx < gt_jobs; tidx++) {
      combinations[tidx] = gt_array_new(sizeof (GtUwordPair));
    }
    for (aidx = 0; aidx < arg->anumseqranges; aidx++) {
      if (apick && pick->a != aidx) continue;
      for (bidx = self ? aidx : 0; bidx < arg->bnumseqranges; bidx++) {
        if (!bpick || pick->b == bidx) {
          GtUwordPair comb = {aidx, bidx};
          gt_array_add(combinations[counter++ % gt_jobs], comb);
        }
      }
    }

    for (tidx = 1; !had_err && tidx < gt_jobs; tidx++) {
      GtThread *thread;
      gt_diagbandseed_thread_info_set(tinfo + tidx,
                                      arg,
                                      NULL,
                                      stream[tidx],
                                      aseqranges,
                                      bseqranges,
                                      combinations[tidx],
                                      err);
      if ((thread = gt_thread_new(gt_diagbandseed_thread_algorithm,
                                  tinfo + tidx, err)) != NULL) {
        gt_array_add(threads, thread);
      } else {
        had_err = -1;
      }
    }
    /* start main thread */
    if (!had_err) {
      gt_diagbandseed_thread_info_set(tinfo,
                                      arg,
                                      NULL,
                                      stream[0],
                                      aseqranges,
                                      bseqranges,
                                      combinations[0],
                                      err);
      gt_diagbandseed_thread_algorithm(tinfo);
    }

    /* clean up */
    for (tidx = 0; tidx < gt_array_size(threads); tidx++) {
      GtThread *thread = *(GtThread**) gt_array_get(threads, tidx);
      if (!had_err) {
        gt_thread_join(thread);
      }
      gt_thread_delete(thread);
    }
    for (tidx = 0; tidx < gt_array_size(threads) && !had_err; tidx++) {
      had_err = tinfo[tidx].had_err;
    }
    gt_array_delete(threads);
    for (tidx = 0; tidx < gt_jobs; tidx++) {
      gt_array_delete(combinations[tidx]);
    }
  }
  gt_free(tinfo);

  /* print the threads' output to stdout */
  for (tidx = 1; tidx < gt_jobs; tidx++) {
    int cc;
    rewind(stream[tidx]);
    while ((cc = fgetc(stream[tidx])) != EOF) {
      putchar(cc);
    }
    gt_fa_xfclose(stream[tidx]);
  }
  gt_free(stream);

#endif
  gt_ft_trimstat_delete(trimstat,arg->verbose);
  return had_err;
}
