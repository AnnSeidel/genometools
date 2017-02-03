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
#include <sys/time.h>
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
#include "core/spacecalc.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "match/declare-readfunc.h"
#include "match/kmercodes.h"
#include "match/querymatch.h"
#include "match/querymatch-align.h"
#include "match/seed-extend.h"
#include "match/sfx-mappedstr.h"
#include "match/sfx-suffixer.h"
#include "match/diagbandseed.h"

#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

/* We need to use 6 digits for the micro seconds */
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
  const GtEncseq *aencseq,
                 *bencseq;
  const GtDiagbandseedExtendParams *extp;
  const GtRange *seedpairdistance;
  GtUword maxfreq,
          memlimit;
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
          sensitivity;
  GtXdropscore xdropbelowscore;
  GtExtendCharAccess a_extend_char_access,
                     b_extend_char_access;
  const GtSeedExtendDisplayFlag *display_flag;
  double matchscore_bias;
  bool use_apos,
       extendgreedy,
       extendxdrop,
       weakends,
       benchmark,
       always_polished_ends,
       verify_alignment,
       only_selected_seqpairs,
       cam_generic;
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
                                             const GtRange *seedpairdistance,
                                             GtDiagbandseedPairlisttype splt,
                                             bool verify,
                                             bool verbose,
                                             bool debug_kmer,
                                             bool debug_seedpair,
                                             bool use_kmerfile,
                                             bool trimstat_on,
                                             const GtDiagbandseedExtendParams
                                               *extp)
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
                                const GtSeedExtendDisplayFlag *display_flag,
                                bool use_apos,
                                GtXdropscore xdropbelowscore,
                                bool extendgreedy,
                                bool extendxdrop,
                                GtUword maxalignedlendifference,
                                GtUword history_size,
                                GtUword perc_mat_history,
                                GtExtendCharAccess a_extend_char_access,
                                GtExtendCharAccess b_extend_char_access,
                                bool cam_generic,
                                GtUword sensitivity,
                                double matchscore_bias,
                                bool weakends,
                                bool benchmark,
                                bool always_polished_ends,
                                bool verify_alignment,
                                bool only_selected_seqpairs)
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
  extp->a_extend_char_access = a_extend_char_access;
  extp->b_extend_char_access = b_extend_char_access;
  extp->cam_generic = cam_generic;
  extp->sensitivity = sensitivity;
  extp->matchscore_bias = matchscore_bias;
  extp->weakends = weakends;
  extp->benchmark = benchmark;
  extp->always_polished_ends = always_polished_ends;
  extp->verify_alignment = verify_alignment;
  extp->only_selected_seqpairs = only_selected_seqpairs;
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
                                         GtUword seqrange_start,
                                         GtUword seqrange_end)
{
  GtUword lastpos, numofpos, subtract, ratioofspecial;

  const GtUword totalnumofspecial = gt_encseq_specialcharacters(encseq),
                totalnumofpos = gt_encseq_total_length(encseq),
                firstpos = gt_encseq_seqstartpos(encseq, seqrange_start),
                numofseq = seqrange_end - seqrange_start + 1;
  lastpos = (seqrange_end + 1 == gt_encseq_num_of_sequences(encseq)
             ? totalnumofpos
             : gt_encseq_seqstartpos(encseq, seqrange_end + 1) - 1);
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
static GtArrayGtDiagbandseedKmerPos gt_diagbandseed_get_kmers(
                                   const GtEncseq *encseq,
                                   unsigned int seedlength,
                                   GtReadmode readmode,
                                   GtUword seqrange_start,
                                   GtUword seqrange_end,
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
    listlen = gt_seed_extend_numofkmers(encseq, seedlength, seqrange_start,
                                        seqrange_end);
  }
  if (verbose) {
    timer = gt_timer_new();
    fprintf(stream, "# start fetching %u-mers (expect " GT_WU
                    ", allocate %.0f MB) ...\n",
            seedlength, listlen,
            GT_MEGABYTES(listlen * sizeof (*list.spaceGtDiagbandseedKmerPos)));
    gt_timer_start(timer);
  }

  GT_INITARRAY(&list, GtDiagbandseedKmerPos);
  GT_CHECKARRAYSPACEMULTI(&list, GtDiagbandseedKmerPos, listlen);

  pkinfo.list = &list;
  pkinfo.seqnum = seqrange_start;
  pkinfo.endpos = 0;
  pkinfo.encseq = encseq;
  pkinfo.seedlength = seedlength;
  pkinfo.readmode = readmode;
  if (seqrange_end + 1 == gt_encseq_num_of_sequences(encseq)) {
    pkinfo.last_specialpos = totallength;
  } else {
    /* start position of following sequence, minus separator position */
    pkinfo.last_specialpos
      = gt_encseq_seqstartpos(encseq, seqrange_end + 1) - 1;
  }
  pkinfo.prev_separator = gt_encseq_seqstartpos(encseq, seqrange_start);
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
  list.allocatedGtDiagbandseedKmerPos = listlen;
  list.spaceGtDiagbandseedKmerPos
    = gt_realloc(list.spaceGtDiagbandseedKmerPos,
                 listlen * sizeof (*list.spaceGtDiagbandseedKmerPos));

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
    fprintf(stream, "# ... collected " GT_WU " %u-mers ", listlen, seedlength);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_start(timer);
  }

  /* sort list */
  gt_radixsort_inplace_GtUwordPair((GtUwordPair *)
                                   list.spaceGtDiagbandseedKmerPos,
                                   listlen);
  if (verbose) {
    fprintf(stream, "# ... sorted " GT_WU " %u-mers ", listlen, seedlength);
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

static GtDiagbandseedKmerIterator *gt_diagbandseed_kmer_iter_new_file(FILE *fp)
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

/* Evaluate the results of the seed count histogram */
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
    /* count seeds until available memory reached */
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

const char *gt_diagbandseed_splt_comment(void)
{
  return "specify type of pairlist, possible values are struct, bytestring, "
         "and ulong";
}

static const char *gt_splt_arguments[] = {"struct","ulong","bytestring",""};

GtDiagbandseedPairlisttype gt_diagbandseed_splt_get(const char *splt_string,
                                                    GtError *err)
{
  size_t idx;
  for (idx = 0; idx < sizeof gt_splt_arguments/sizeof gt_splt_arguments[0];
       idx++)
  {
    if (strcmp(splt_string,gt_splt_arguments[idx]) == 0)
    {
      return (GtDiagbandseedPairlisttype) idx;
    }
  }
  gt_error_set(err,"illegal parameter for option -splt: %s",
                    gt_diagbandseed_splt_comment());
  return -1;
}

const int idx_aseqnum = 0, idx_bseqnum = 1, idx_bpos = 2, idx_apos = 3;

typedef struct
{
  GtDiagbandseedPosition apos, bpos;
} GtSeedpairPositions;

GT_DECLAREARRAYSTRUCT(GtDiagbandseedSeedPair);

typedef struct
{
  GtArrayGtDiagbandseedSeedPair *mlist_struct;
  GtArrayGtUword *mlist_ulong;
  GtArrayuint8_t *mlist_bytestring;
  GtDiagbandseedPairlisttype splt;
  GtUword mask_tab[4], transfer_mask,
          aseqrange_start, bseqrange_start,
          aseqrange_end, bseqrange_end,
          aseqrange_max_length, bseqrange_max_length;
  GtBitcount_type bits_seedpair,
                  bits_values[4],
                  bits_units[2],
                  bits_left_adjust[4];
  size_t bytes_seedpair;
  int shift_tab[4],
      transfer_shift,
      bits_unused_in2GtUwords;
} GtSeedpairlist;

#define GT_DIAGBANDSEED_ENCODE_SEQNUMS(ASEQNUM,BSEQNUM)\
             (((GtUword) (ASEQNUM) << seedpairlist->shift_tab[idx_aseqnum]) | \
              ((GtUword) (BSEQNUM) << seedpairlist->shift_tab[idx_bseqnum]))

#define GT_DIAGBANDSEED_ENCODE_POSITIONS(BPOS,APOS)\
             (((GtUword) (BPOS) << seedpairlist->shift_tab[idx_bpos]) | \
              ((GtUword) (APOS) << seedpairlist->shift_tab[idx_apos]))

#define GT_DIAGBANDSEED_ENCODE_SEEDPAIR(ASEQNUM,BSEQNUM,BPOS,APOS)\
        (GT_DIAGBANDSEED_ENCODE_SEQNUMS(ASEQNUM,BSEQNUM) | \
         GT_DIAGBANDSEED_ENCODE_POSITIONS(BPOS,APOS))

static GtSeedpairlist *gt_seedpairlist_new(GtDiagbandseedPairlisttype splt,
                        const GtSequencePartsInfo *aseqranges,
                        GtUword aidx,
                        const GtSequencePartsInfo *bseqranges,
                        GtUword bidx)
{
  GtSeedpairlist *seedpairlist = gt_malloc(sizeof *seedpairlist);
  int idx;
  const int allbits = sizeof (GtWord) * CHAR_BIT;
  const GtUword anumofseq
    = gt_sequence_parts_info_numofsequences_get(aseqranges,aidx),
    bnumofseq = gt_sequence_parts_info_numofsequences_get(bseqranges,bidx);

  seedpairlist->bits_values[idx_aseqnum]
    = (GtBitcount_type) gt_radixsort_bits(anumofseq);
  seedpairlist->bits_values[idx_bseqnum]
    = (GtBitcount_type) gt_radixsort_bits(bnumofseq);
  seedpairlist->bits_values[idx_bpos]
    = (GtBitcount_type)
      gt_radixsort_bits(gt_sequence_parts_info_max_length_get(bseqranges,bidx));
  seedpairlist->bits_values[idx_apos]
    = (GtBitcount_type)
      gt_radixsort_bits(gt_sequence_parts_info_max_length_get(aseqranges,aidx));
  seedpairlist->bits_units[0] = seedpairlist->bits_values[idx_aseqnum] +
                                seedpairlist->bits_values[idx_bseqnum];
  seedpairlist->bits_units[1] = seedpairlist->bits_values[idx_apos] +
                                seedpairlist->bits_values[idx_bpos];
  seedpairlist->bits_left_adjust[idx_aseqnum]
    = allbits - seedpairlist->bits_values[idx_aseqnum];
  seedpairlist->bits_left_adjust[idx_bseqnum]
    = allbits - seedpairlist->bits_units[0];
  seedpairlist->bits_left_adjust[idx_bpos]
    = allbits - seedpairlist->bits_values[idx_bpos];
  seedpairlist->bits_left_adjust[idx_apos]
    = allbits - seedpairlist->bits_units[1];
  seedpairlist->transfer_mask
    = (((GtUword) 1) << seedpairlist->bits_left_adjust[idx_bseqnum]) - 1;
  gt_assert(seedpairlist->bits_values[idx_apos] > 0);
  gt_assert(seedpairlist->bits_values[idx_bpos] > 0);
  for (idx = 0, seedpairlist->bits_seedpair = 0; idx < 4; idx++)
  {
    seedpairlist->bits_seedpair += seedpairlist->bits_values[idx];
    seedpairlist->mask_tab[idx]
      = (((GtUword) 1) << seedpairlist->bits_values[idx]) - 1;
  }
  /* the following is only used for bits_seedpair > allbits */
  seedpairlist->transfer_shift = seedpairlist->bits_seedpair - allbits;
  seedpairlist->bits_unused_in2GtUwords = 2 * allbits -
                                          seedpairlist->bits_seedpair;
  seedpairlist->bytes_seedpair
    = gt_radixsort_bits2bytes(seedpairlist->bits_seedpair);
  if (seedpairlist->bytes_seedpair <= sizeof (GtUword))
  {
    int shift = seedpairlist->bits_seedpair;
    for (idx = 0; idx < 4; idx++)
    {
      shift -= seedpairlist->bits_values[idx];
      seedpairlist->shift_tab[idx] = shift;
    }
  }
  seedpairlist->aseqrange_start
    = gt_sequence_parts_info_start_get(aseqranges,aidx);
  seedpairlist->bseqrange_start
    = gt_sequence_parts_info_start_get(bseqranges,bidx);
  seedpairlist->aseqrange_end
    = gt_sequence_parts_info_end_get(aseqranges,aidx);
  seedpairlist->bseqrange_end
    = gt_sequence_parts_info_end_get(bseqranges,bidx);
  seedpairlist->aseqrange_max_length
    = gt_sequence_parts_info_max_length_get(aseqranges,aidx);
  seedpairlist->bseqrange_max_length
    = gt_sequence_parts_info_max_length_get(bseqranges,bidx);
  seedpairlist->mlist_struct = NULL;
  seedpairlist->mlist_ulong = NULL;
  seedpairlist->mlist_bytestring = NULL;
  if (splt == GT_DIAGBANDSEED_SPLT_UNDEFINED)
  {
    if (seedpairlist->bytes_seedpair <= sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_SPLT_ULONG;
    } else
    {
      splt = GT_DIAGBANDSEED_SPLT_BYTESTRING;
    }
  }
  if (splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    if (seedpairlist->bytes_seedpair > sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_SPLT_BYTESTRING;
    }
  } else
  {
    if (splt == GT_DIAGBANDSEED_SPLT_BYTESTRING &&
        seedpairlist->bytes_seedpair <= sizeof (GtUword))
    {
      splt = GT_DIAGBANDSEED_SPLT_ULONG;
    }
  }
  if (splt == GT_DIAGBANDSEED_SPLT_ULONG)
  {
    gt_assert(seedpairlist->bytes_seedpair <= sizeof (GtUword));
    seedpairlist->mlist_ulong = gt_malloc(sizeof *seedpairlist->mlist_ulong);
    GT_INITARRAY(seedpairlist->mlist_ulong, GtUword);
  } else
  {
    if (splt == GT_DIAGBANDSEED_SPLT_BYTESTRING)
    {
      gt_assert(seedpairlist->bytes_seedpair > sizeof (GtUword));
      seedpairlist->mlist_bytestring
        = gt_malloc(sizeof *seedpairlist->mlist_bytestring);
      GT_INITARRAY(seedpairlist->mlist_bytestring, uint8_t);
    } else
    {
      seedpairlist->mlist_struct
        = gt_malloc(sizeof *seedpairlist->mlist_struct);
      GT_INITARRAY(seedpairlist->mlist_struct, GtDiagbandseedSeedPair);
    }
  }
  seedpairlist->splt = splt;
  return seedpairlist;
}

static void gt_seedpairlist_reset(GtSeedpairlist *seedpairlist)
{
  if (seedpairlist->mlist_ulong != NULL)
  {
    seedpairlist->mlist_ulong->nextfreeGtUword = 0;
  }
  if (seedpairlist->mlist_bytestring != NULL)
  {
    seedpairlist->mlist_bytestring->nextfreeuint8_t = 0;
  }
  if (seedpairlist->mlist_struct != NULL)
  {
    seedpairlist->mlist_struct->nextfreeGtDiagbandseedSeedPair = 0;
  }
}

static void gt_seedpairlist_show_bits(FILE *stream,
                                      const GtSeedpairlist *seedpairlist)
{
  fprintf(stream,"# splt=%s, bits_seedpair=%d, bytes_seedpair=%u with ",
                    gt_splt_arguments[seedpairlist->splt],
                    (int) seedpairlist->bits_seedpair,
                    (int) seedpairlist->bytes_seedpair);
  fprintf(stream,"aseqnum=%hu bits, ",seedpairlist->bits_values[idx_aseqnum]);
  fprintf(stream,"bseqnum=%hu bits, ",seedpairlist->bits_values[idx_bseqnum]);
  fprintf(stream,"bpos=%hu bits, ",seedpairlist->bits_values[idx_bpos]);
  fprintf(stream,"apos=%hu bits\n",seedpairlist->bits_values[idx_apos]);
}

static size_t gt_seedpairlist_sizeofunit(const GtSeedpairlist *seedpairlist)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    return sizeof (GtDiagbandseedSeedPair);
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      return sizeof (GtUword);
    } else
    {
      return seedpairlist->bytes_seedpair;
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
        gt_assert(seedpairlist->mlist_bytestring != NULL);
        GT_CHECKARRAYSPACEMULTI(seedpairlist->mlist_bytestring,uint8_t,
                                seedpairlist->bytes_seedpair * known_size);
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
    if (seedpairlist->mlist_bytestring != NULL)
    {
      GT_FREEARRAY(seedpairlist->mlist_bytestring, uint8_t);
      gt_free(seedpairlist->mlist_bytestring);
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
      gt_assert(seedpairlist->mlist_bytestring != NULL &&
                (seedpairlist->mlist_bytestring->nextfreeuint8_t %
                 seedpairlist->bytes_seedpair == 0));
      return seedpairlist->mlist_bytestring->nextfreeuint8_t/
             seedpairlist->bytes_seedpair;
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
  gt_assert(seedpairlist != NULL && seedpairlist->mlist_ulong != NULL);
  return seedpairlist->mlist_ulong->spaceGtUword;
}

const uint8_t *gt_seedpairlist_mlist_bytestring(
              const GtSeedpairlist *seedpairlist)
{
  gt_assert(seedpairlist != NULL && seedpairlist->mlist_bytestring != NULL);
  return seedpairlist->mlist_bytestring->spaceuint8_t;
}

static void gt_diagbandseed_one_GtUword2bytestring(uint8_t *bytestring,
                                                   GtUword bytestring_length,
                                                   GtUword value)
{
  int idx, shift;

  for (idx = 0, shift = (sizeof (GtUword) - 1) * CHAR_BIT;
       idx < bytestring_length; idx++, shift -= CHAR_BIT)
  {
    bytestring[idx] = (uint8_t) (value >> shift);
  }
}

static void gt_diagbandseed_encode_seedpair(uint8_t *bytestring,
                                            const GtSeedpairlist *seedpairlist,
                                            GtDiagbandseedSeqnum aseqnum,
                                            GtDiagbandseedSeqnum bseqnum,
                                            GtDiagbandseedPosition bpos,
                                            GtDiagbandseedPosition apos)
{
  GtUword value_seqnums, value_positions;

  value_seqnums
    = (((GtUword) aseqnum) << seedpairlist->bits_left_adjust[idx_aseqnum]) |
      (((GtUword) bseqnum) << seedpairlist->bits_left_adjust[idx_bseqnum]);
  value_positions
    = (((GtUword) bpos) << seedpairlist->bits_left_adjust[idx_bpos]) |
      (((GtUword) apos) << seedpairlist->bits_left_adjust[idx_apos]);
  value_seqnums |= (value_positions >> seedpairlist->bits_units[0]);
  gt_diagbandseed_one_GtUword2bytestring(bytestring,sizeof (GtUword),
                                         value_seqnums);
  value_positions <<= seedpairlist->bits_left_adjust[idx_bseqnum];
  gt_diagbandseed_one_GtUword2bytestring(bytestring + sizeof (GtUword),
                                         seedpairlist->bytes_seedpair
                                         - sizeof (GtUword),
                                         value_positions);
}

static void gt_seedpairlist_add(GtSeedpairlist *seedpairlist,
                                bool knownsize,
                                GtDiagbandseedSeqnum aseqnum,
                                GtDiagbandseedSeqnum bseqnum,
                                GtDiagbandseedPosition bpos,
                                GtDiagbandseedPosition apos)
{
  gt_assert(seedpairlist != NULL);
  gt_assert(aseqnum >= seedpairlist->aseqrange_start &&
            aseqnum <= seedpairlist->aseqrange_end &&
            bseqnum >= seedpairlist->bseqrange_start &&
            bseqnum <= seedpairlist->bseqrange_end &&
            apos < seedpairlist->aseqrange_max_length &&
            bpos < seedpairlist->bseqrange_max_length);
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
        = GT_DIAGBANDSEED_ENCODE_SEEDPAIR(aseqnum -
                                          seedpairlist->aseqrange_start,
                                          bseqnum -
                                          seedpairlist->bseqrange_start,
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
       uint8_t *bytestring;

       if (knownsize)
       {
         gt_assert(seedpairlist->mlist_bytestring->nextfreeuint8_t +
                   seedpairlist->bytes_seedpair <=
                   seedpairlist->mlist_bytestring->allocateduint8_t);
       } else
       {
         GT_CHECKARRAYSPACE_GENERIC(seedpairlist->mlist_bytestring,
                                    uint8_t,
                                    seedpairlist->bytes_seedpair,
                                    256 + 0.2 *
                                    seedpairlist->
                                    mlist_bytestring->allocateduint8_t);
       }
       bytestring = seedpairlist->mlist_bytestring->spaceuint8_t +
                    seedpairlist->mlist_bytestring->nextfreeuint8_t;
       seedpairlist->mlist_bytestring->nextfreeuint8_t
         += seedpairlist->bytes_seedpair;
       gt_diagbandseed_encode_seedpair(bytestring,
                                       seedpairlist,
                                       aseqnum - seedpairlist->aseqrange_start,
                                       bseqnum - seedpairlist->bseqrange_start,
                                       bpos,
                                       apos);
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
      gt_assert(seedpairlist->mlist_struct != NULL);
      gt_radixsort_inplace_Gtuint64keyPair(
            (Gtuint64keyPair *) seedpairlist->mlist_struct
                                            ->spaceGtDiagbandseedSeedPair,
            mlistlen);
    } else
    {
      if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
      {
        gt_assert(seedpairlist->mlist_ulong != NULL);
        gt_radixsort_inplace_ulong(seedpairlist->mlist_ulong->spaceGtUword,
                                   mlistlen);
      } else
      {
        gt_assert(seedpairlist->mlist_bytestring != NULL);
        gt_radixsort_inplace_flba(seedpairlist->mlist_bytestring->spaceuint8_t,
                                  gt_seedpairlist_length(seedpairlist),
                                  seedpairlist->bytes_seedpair);
      }
    }
  }
}

static GtUword gt_seedpairlist_a_bseqnum_ulong(
                                          const GtSeedpairlist *seedpairlist,
                                          GtUword encoding)
{
  return encoding >> seedpairlist->shift_tab[idx_bseqnum];
}

static GtUword gt_seedpairlist_extract_ulong(const GtSeedpairlist *seedpairlist,
                                             GtUword encoding,
                                             int compidx)
{
  return (encoding >> seedpairlist->shift_tab[compidx]) &
         seedpairlist->mask_tab[compidx];
}

static GtUword gt_seedpairlist_extract_ulong_at(
                                   const GtSeedpairlist *seedpairlist,
                                   GtUword spidx,
                                   int compidx)
{
  return gt_seedpairlist_extract_ulong(seedpairlist,
                                       seedpairlist->mlist_ulong->
                                                     spaceGtUword[spidx],
                                       compidx);
}

static GtUword gt_diagbandseed_bytestring2GtUword(const uint8_t *bytestring,
                                                  GtUword bytestring_length)
{
  int idx, shift;
  GtUword value = 0;

  for (idx = 0, shift = (sizeof (GtUword) - 1) * CHAR_BIT;
       idx < bytestring_length; idx++, shift -= CHAR_BIT)
  {
    value |= (((GtUword) bytestring[idx]) << shift);
  }
  return value;
}

static void gt_diagbandseed_decode_seedpair(GtDiagbandseedSeedPair *seedpair,
                                            const GtSeedpairlist *seedpairlist,
                                            GtUword offset)
{
  GtUword transfer, value_seqnums, value_positions;
  const uint8_t *bytestring = seedpairlist->mlist_bytestring->spaceuint8_t +
                              offset;

  value_seqnums
    = gt_diagbandseed_bytestring2GtUword(bytestring,sizeof (GtUword));
  seedpair->aseqnum = (GtDiagbandseedSeqnum)
                      (value_seqnums >>
                       seedpairlist->bits_left_adjust[idx_aseqnum]);
  seedpair->bseqnum = (GtDiagbandseedSeqnum)
                       ((value_seqnums >>
                         seedpairlist->bits_left_adjust[idx_bseqnum]) &
                         seedpairlist->mask_tab[idx_bseqnum]);
  transfer = (value_seqnums & seedpairlist->transfer_mask);
  value_positions
    = gt_diagbandseed_bytestring2GtUword(bytestring + sizeof (GtUword),
                                         seedpairlist->bytes_seedpair -
                                         sizeof (GtUword));
  value_positions >>= seedpairlist->bits_unused_in2GtUwords;
  value_positions |= (transfer << seedpairlist->transfer_shift);
  seedpair->bpos = (GtDiagbandseedPosition)
                    (value_positions >> seedpairlist->bits_values[idx_apos]);
  seedpair->apos = (GtDiagbandseedPosition)
                    (value_positions & seedpairlist->mask_tab[idx_apos]);
}

static void gt_seedpairlist_at(GtDiagbandseedSeedPair *seedpair,
                               const GtSeedpairlist *seedpairlist,
                               GtUword spidx)
{
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    *seedpair = seedpairlist->mlist_struct->spaceGtDiagbandseedSeedPair[spidx];
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      seedpair->aseqnum
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_aseqnum);
      seedpair->bseqnum
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_bseqnum);
      seedpair->bpos
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_bpos);
      seedpair->apos
        = gt_seedpairlist_extract_ulong_at(seedpairlist,spidx,idx_apos);
    } else
    {
      gt_diagbandseed_decode_seedpair(seedpair,seedpairlist,
                                      spidx * seedpairlist->bytes_seedpair);
    }
  }
}

static void gt_diagbandseed_show_seed(FILE *stream,
                                      const GtSeedpairlist *seedpairlist,
                                      GtUword spidx)
{
  GtDiagbandseedSeedPair seedpair;

  gt_seedpairlist_at(&seedpair,seedpairlist,spidx);
  fprintf(stream, "(%"PRIu32 ",%"PRIu32 ",%"PRIu32 ",%"PRIu32")",
                  seedpair.aseqnum,
                  seedpair.bseqnum,
                  seedpair.apos,
                  seedpair.bpos);
}

#ifndef NDEBUG
static int gt_diagbandseed_seeds_compare(const GtSeedpairlist *seedpairlist,
                                         const GtUword current)
{
  GtDiagbandseedSeedPair p_seedpair, c_seedpair;

  gt_assert(current > 0);
  gt_seedpairlist_at(&p_seedpair,seedpairlist,current - 1);
  gt_seedpairlist_at(&c_seedpair,seedpairlist,current);
  if (p_seedpair.aseqnum < c_seedpair.aseqnum)
  {
    return -1;
  }
  if (p_seedpair.aseqnum > c_seedpair.aseqnum)
  {
    return 1;
  }
  if (p_seedpair.bseqnum < c_seedpair.bseqnum)
  {
    return -1;
  }
  if (p_seedpair.bseqnum > c_seedpair.bseqnum)
  {
    return 1;
  }
  if (p_seedpair.bpos < c_seedpair.bpos)
  {
    return -1;
  }
  if (p_seedpair.bpos > c_seedpair.bpos)
  {
    return 1;
  }
  if (p_seedpair.apos < c_seedpair.apos)
  {
    return -1;
  }
  if (p_seedpair.apos > c_seedpair.apos)
  {
    return 1;
  }
  return 0;
}
#endif

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

/* Fill a GtDiagbandseedSeedPair list of equal kmers from the iterators. */
static void gt_diagbandseed_merge(GtSeedpairlist *seedpairlist,
                                  GtUword *histogram,
                                  bool knowthesize,
                                  GtDiagbandseedKmerIterator *aiter,
                                  GtDiagbandseedKmerIterator *biter,
                                  GtUword maxfreq,
                                  GtUword maxgram,
                                  const GtRange *seedpairdistance,
                                  bool selfcomp)
{
  const GtArrayGtDiagbandseedKmerPos *alist, *blist;
  const bool count_cartesian = (histogram != NULL && !selfcomp) ? true : false;

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
          if (count_cartesian)
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
                                        aptr->seqnum,
                                        bptr->seqnum,
                                        bptr->endpos,
                                        aptr->endpos);
                  } else
                  {
                    /* count seed frequency in histogram */
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

/* Verify seeds in the original sequences */
static int gt_diagbandseed_verify(const GtSeedpairlist *seedpairlist,
                                  const GtEncseq *aencseq,
                                  const GtEncseq *bencseq,
                                  unsigned int seedlength,
                                  bool reverse,
                                  bool verbose,
                                  FILE *stream,
                                  GtError *err) {
  const GtUword mlistlen = gt_seedpairlist_length(seedpairlist);
  GtTimer *timer = NULL;
  GtUword idx;
  char *buf1 = gt_malloc(3 * (seedlength + 1) * sizeof *buf1);
  char *buf2 = buf1 + 1 + seedlength;
  char *buf3 = buf2 + 1 + seedlength;
  buf1[seedlength] = buf2[seedlength] = buf3[seedlength] = '\0';

  if (verbose) {
    fprintf(stream, "# start verifying seeds ...\n");
    timer = gt_timer_new();
    gt_timer_start(timer);
  }

  gt_assert(aencseq != NULL && bencseq != NULL);
  for (idx = 0; idx < mlistlen; idx++)
  {
    GtDiagbandseedSeedPair seedpair;
    gt_seedpairlist_at(&seedpair,seedpairlist,idx);
    GtDiagbandseedPosition abs_apos, abs_bpos;

    /* extract decoded k-mers at seed positions */
    abs_apos = seedpair.apos + gt_encseq_seqstartpos(aencseq, seedpair.aseqnum);
    gt_encseq_extract_decoded(aencseq, buf1, abs_apos + 1 - seedlength,
                              abs_apos);

    if (!reverse) {
      abs_bpos = seedpair.bpos + gt_encseq_seqstartpos(bencseq,
                                                       seedpair.bseqnum);
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos + 1 - seedlength,
                                abs_bpos);
      if (strcmp(buf1, buf2) != 0) {
        gt_error_set(err, "Wrong SeedPair (" "%"PRIu32
                          ",%"PRIu32",%"PRIu32",%"PRIu32"): %s != %s\n",
                     seedpair.aseqnum, seedpair.bseqnum,
                     seedpair.apos, seedpair.bpos, buf1, buf2);
        gt_free(buf1);
        if (verbose)
        {
          gt_timer_delete(timer);
        }
        return -1;
      }
    } else {
      /* get reverse k-mer */
      char *bufptr;
      abs_bpos = gt_encseq_seqstartpos(bencseq, seedpair.bseqnum) +
                 gt_encseq_seqlength(bencseq, seedpair.bseqnum)
                 - seedpair.bpos - 1;
      gt_encseq_extract_decoded(bencseq, buf2, abs_bpos,
                                abs_bpos + seedlength - 1);

      for (bufptr = buf3; bufptr < buf3 + seedlength; bufptr++) {
        gt_complement(bufptr, buf2[seedlength + buf3 - bufptr - 1], NULL);
      }
      if (strcmp(buf1, buf3) != 0) {
        gt_error_set(err, "Wrong SeedPair (" "%"PRIu32
                          ",%"PRIu32",%"PRIu32",%"PRIu32"): %s != %s\n",
                     seedpair.aseqnum, seedpair.bseqnum,
                     seedpair.apos, seedpair.bpos, buf1, buf3);
        gt_free(buf1);
        if (verbose)
        {
          gt_timer_delete(timer);
        }
        return -1;
      }
    }
  }
  if (verbose) {
    fprintf(stream, "# ...successfully verified all seeds ");
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
    gt_timer_delete(timer);
  }
  gt_free(buf1);
  return 0;
}

/* Return estimated length of mlist, and maxfreq w.r.t. the given memlimit */
static int gt_diagbandseed_get_mlistlen_maxfreq(GtUword *mlistlen,
                                            GtUword *maxfreq,
                                            GtDiagbandseedKmerIterator *aiter,
                                            GtDiagbandseedKmerIterator *biter,
                                            GtUword memlimit,
                                            size_t sizeofunit,
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
    fprintf(stream, "# start calculating k-mer frequency histogram...\n");
    timer = gt_timer_new();
    gt_timer_start(timer);
  }

  /* build histogram; histogram[maxgram] := estimation for mlistlen */
  histogram = gt_calloc(maxgram + 1, sizeof *histogram);
  gt_diagbandseed_merge(NULL, /* mlist not needed: just count */
                        histogram,
                        false,
                        aiter,
                        biter,
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
                                              sizeofunit);
  *mlistlen = histogram[maxgram];
  gt_free(histogram);

  if (verbose) {
    gt_timer_show_formatted(timer,"# ... finished histogram "
                            GT_DIAGBANDSEED_FMT,stream);
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
      fprintf(stream, "# disable k-mer maximum frequency, ");
    } else {
      fprintf(stream, "# set k-mer maximum frequency to " GT_WU ", ", *maxfreq);
    }
    fprintf(stream, "expect " GT_WU " seeds.\n", *mlistlen);
  } else if (*maxfreq <= 5) {
    gt_warning("only k-mers occurring <= " GT_WU " times will be considered, "
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
                                          bool debug_seedpair,
                                          bool verbose,
                                          FILE *stream)
{
  GtTimer *timer = NULL;
  GtUword mlistlen;

  if (verbose) {
    timer = gt_timer_new();
    if (known_size > 0) {
      fprintf(stream, "# start collecting " GT_WU " seeds in %.0f MB ...\n",
                      known_size,
                      GT_MEGABYTES(known_size *
                                   gt_seedpairlist_sizeofunit(seedpairlist)));
    } else {
      fprintf(stream, "# start collecting seeds ...\n");
    }
    gt_timer_start(timer);
  }

  /* allocate mlist space according to seed count */
  gt_seedpairlist_init(seedpairlist,known_size);

  /* create mlist */
  (void) gt_diagbandseed_merge(seedpairlist,
                        NULL, /* histogram not needed: save seeds */
                        known_size > 0 ? true : false,
                        aiter,
                        biter,
                        maxfreq,
                        GT_UWORD_MAX, /* maxgram */
                        seedpairdistance,
                        selfcomp);
  mlistlen = gt_seedpairlist_length(seedpairlist);
  if (verbose) {
    fprintf(stream, "# ... collected " GT_WU " seeds ", mlistlen);
    gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
  }

  /* sort mlist */
  if (mlistlen > 0) {
    if (verbose) {
      gt_timer_start(timer);
    }
    gt_diagbandseed_seedpairlist_sort(seedpairlist);
    if (verbose) {
      fprintf(stream, "# ... sorted " GT_WU " seeds ", mlistlen);
      gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
      gt_timer_delete(timer);
    }
    if (debug_seedpair) {
      gt_diagbandseed_seedpairlist_out(stream,seedpairlist);
    }
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
    /* overlap: add positions after last counted position */
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

typedef struct
{
  GtUword extended_seeds,
          selected_seeds,
          countmatches,
          failedmatches,
          seqpairs_with_minsegment;
  bool withtiming;
  GtUword total_extension_time_usec;
} GtDiagbandseedCounts;

typedef const GtQuerymatch *(*GtExtendRelativeCoordsFunc)(void *,
                                                          const GtSeqorEncseq *,
                                                          GtUword,
                                                          GtUword,
                                                          const GtSeqorEncseq *,
                                                          bool,
                                                          GtUword,
                                                          GtUword,
                                                          GtUword,
                                                          GtReadmode);

static int gt_diagbandseed_possibly_extend(const GtQuerymatch *previousmatch,
                                           GtUword aseqnum,
                                           GtUword apos,
                                           GtUword bseqnum,
                                           GtUword bpos,
                                           bool use_apos,
                                           unsigned int seedlength,
                                           GtUword errorpercentage,
                                           GtUword userdefinedleastlength,
                                           const GtSeqorEncseq *aseqorencseq,
                                           const GtSeqorEncseq *bseqorencseq,
                                           bool same_encseq,
                                           GtProcessinfo_and_querymatchspaceptr
                                             *info_querymatch,
                                           GtReadmode query_readmode,
                                           GtExtendRelativeCoordsFunc
                                             extend_relative_coords_function,
                                           GtDiagbandseedCounts
                                             *process_seeds_counts)
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
#ifndef _WIN32
    struct timeval tvalBefore;

    if (process_seeds_counts != NULL)
    {
      gettimeofday (&tvalBefore, NULL);
    }
#endif
    ret = 1; /* perform extension */
    querymatch = extend_relative_coords_function(info_querymatch,
                                                 aseqorencseq,
                                                 aseqnum,
                                                 astart,
                                                 bseqorencseq,
                                                 same_encseq,
                                                 bseqnum,
                                                 bstart,
                                                 seedlength,
                                                 query_readmode);
#ifndef _WIN32
    if (process_seeds_counts != NULL)
    {
      struct timeval tvalAfter;
      gettimeofday (&tvalAfter, NULL);
      process_seeds_counts->total_extension_time_usec
        += (tvalAfter.tv_sec - tvalBefore.tv_sec) * 1000000L
            + tvalAfter.tv_usec - tvalBefore.tv_usec;
    }
#endif
    if (querymatch != NULL)
    {
      /* show extension results */
      if (gt_querymatch_check_final(querymatch, errorpercentage,
                                    userdefinedleastlength))
      {
        gt_querymatch_prettyprint(querymatch);
        ret = 3; /* output match */
      } else
      {
        ret = 2; /* found match, which does not satisfy length or similarity
                    constraints */
      }
    }
  }
  return ret;
}

/* * * * * SEED EXTENSION * * * * */

#define GT_DIAGBANDSEED_DIAG(AMAXLEN,APOS,BPOS)\
        (((AMAXLEN) + (GtUword) (BPOS) - (GtUword) (APOS)) \
         >> extp->logdiagbandwidth)

static void gt_diagbandseed_process_segment(
             const GtDiagbandseedExtendParams *extp,
             const GtSeqorEncseq *aseqorencseq,
             const GtSeqorEncseq *bseqorencseq,
             bool same_encseq,
             GtUword amaxlen,
             unsigned int seedlength,
             GtUword diagbands_used,
             GtUword ndiags,
             GtDiagbandseedScore *diagband_score,
             GtDiagbandseedPosition *diagband_lastpos,
             GtUword aseqnum,
             GtUword bseqnum,
             const GtSeedpairPositions *segment_positions,
             GtUword segment_length,
             GtProcessinfo_and_querymatchspaceptr *info_querymatch,
             GtReadmode query_readmode,
             GtExtendRelativeCoordsFunc extend_relative_coords_function,
             GtDiagbandseedCounts *process_seeds_counts)
{
  const GtSeedpairPositions *spp_ptr;
  bool haspreviousmatch = false, found_selected = false;

  for (spp_ptr = segment_positions;
       spp_ptr < segment_positions + segment_length;
       spp_ptr++)
  {
    const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen, spp_ptr->apos,
                                              spp_ptr->bpos);

    if ((GtUword) MAX(diagband_score[diag + 1], diagband_score[diag - 1])
        + (GtUword) diagband_score[diag]
        >= extp->mincoverage)
    {
      int ret;

      process_seeds_counts->selected_seeds++;
      found_selected = true;
      if (extp->only_selected_seqpairs)
      {
        printf("# " GT_WU "%c" GT_WU "\n",aseqnum,
               query_readmode == GT_READMODE_REVCOMPL ? '-' : '+',
               bseqnum);
        break;
      }
      ret = gt_diagbandseed_possibly_extend(
                     haspreviousmatch ? info_querymatch->querymatchspaceptr
                                      : NULL,
                     aseqnum,
                     spp_ptr->apos,
                     bseqnum,
                     spp_ptr->bpos,
                     extp->use_apos,
                     seedlength,
                     extp->errorpercentage,
                     extp->userdefinedleastlength,
                     aseqorencseq,
                     bseqorencseq,
                     same_encseq,
                     info_querymatch,
                     query_readmode,
                     extend_relative_coords_function,
                     process_seeds_counts->withtiming ? process_seeds_counts
                                                      : NULL);
      if (ret >= 1)
      {
        process_seeds_counts->extended_seeds++;
      }
      if (ret >= 2)
      {
        haspreviousmatch = true;
        if (ret == 2)
        {
          process_seeds_counts->failedmatches++;
        }
      }
      if (ret == 3)
      {
        process_seeds_counts->countmatches++;
      }
    }
  }
  if (found_selected)
  {
    process_seeds_counts->seqpairs_with_minsegment++;
  }
  /* reset diagonal band scores */
  if (diagbands_used * 3 >= ndiags) /* >= 33% of diagbands are used */
  {
    memset(diagband_score,0,sizeof *diagband_score * ndiags);
    memset(diagband_lastpos,0,sizeof *diagband_lastpos * ndiags);
  } else
  {
    /* third scan: */
    for (spp_ptr = segment_positions;
         spp_ptr < segment_positions + segment_length;
         spp_ptr++)
    {
      const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen,
                                                spp_ptr->apos,
                                                spp_ptr->bpos);
      diagband_score[diag] = 0;
      diagband_lastpos[diag] = 0;
    }
  }
}

static void gt_diagbandseed_match_header(FILE *stream,
                                         const GtDiagbandseedExtendParams *extp,
                                         unsigned int seedlength,
                                         GtUword ndiags,
                                         GtUword minsegmentlen)
{
  GtStr *add_column_header;

  fprintf(stream,"# start processing of seeds ...\n");
  fprintf(stream,"# parameters for selecting seeds: seedlength=%u, "
                 "diagonal bands=" GT_WU ", minimal segmentsize=" GT_WU
                 ", minimal coverage=" GT_WU "\n",
                 seedlength,ndiags,minsegmentlen,extp->mincoverage);
  fprintf(stream,"# columns: alen aseq astartpos strand blen bseq bstartpos "
                 "score editdist identity");
  add_column_header = gt_querymatch_column_header(extp->display_flag);
  if (gt_str_length(add_column_header) > 0)
  {
    fputs(gt_str_get(add_column_header),stream);
  }
  fputc('\n',stream);
  gt_str_delete(add_column_header);
}

static void gt_diagbandseed_info_qm_set(
                                   GtProcessinfo_and_querymatchspaceptr *ifqm,
                                   const GtDiagbandseedExtendParams *extp,
                                   const GtEncseq *aencseq,
                                   GtQuerymatchoutoptions *querymoutopt,
                                   GtReadmode query_readmode,
                                   FILE *stream,
                                   void *processinfo)
{
  ifqm->processinfo = processinfo;
  if (gt_querymatch_evalue_display(extp->display_flag) ||
      gt_querymatch_bit_score_display(extp->display_flag))
  {
    ifqm->karlin_altschul_stat = gt_karlin_altschul_stat_new_gapped();
    gt_karlin_altschul_stat_add_keyvalues(ifqm->karlin_altschul_stat,
                                          gt_encseq_total_length(aencseq),
                                          gt_encseq_num_of_sequences(aencseq));
  } else
  {
    ifqm->karlin_altschul_stat = NULL;
  }
  ifqm->querymatchspaceptr = gt_querymatch_new();
  gt_querymatch_display_set(ifqm->querymatchspaceptr,extp->display_flag);
  if (extp->verify_alignment)
  {
    gt_querymatch_verify_alignment_set(ifqm->querymatchspaceptr);
  }
  if (querymoutopt != NULL) {
    gt_querymatch_outoptions_set(ifqm->querymatchspaceptr,querymoutopt);
  }
  gt_querymatch_query_readmode_set(ifqm->querymatchspaceptr,query_readmode);
  gt_querymatch_file_set(ifqm->querymatchspaceptr, stream);
}

#define GT_USEC2SEC(TIME_IN_USEC)\
         ((GtUword) (TIME_IN_USEC)/1000000)

#define GT_USECREMAIN(TIME_IN_USEC) ((TIME_IN_USEC) -\
                                     GT_USEC2SEC(TIME_IN_USEC) * 1000000UL)

#ifndef _WIN32
static void gt_diagbandseed_process_seeds_times(
                 FILE *stream,
                 bool extendgreedy,
                 GtUword mlistlen,
                 GtUword extended_seeds,
                 GtUword total_process_seeds_usec,
                 GtUword total_extension_time_usec)
{
  GtUword process_seeds_usec;

  gt_assert(total_process_seeds_usec >= total_extension_time_usec);
  process_seeds_usec = total_process_seeds_usec - total_extension_time_usec;
  fprintf(stream,"# ... %s extension of " GT_WU " seeds in "
                 GT_WD ".%.06ld seconds.\n",
          extendgreedy ? "greedy" : "xdrop",
          extended_seeds,
          GT_USEC2SEC(total_extension_time_usec),
          GT_USECREMAIN(total_extension_time_usec));
  fprintf(stream, "# ... processed " GT_WU " seeds (excluding extensions) in "
                  GT_WD ".%06ld seconds.\n",mlistlen,
          GT_USEC2SEC(process_seeds_usec),
          GT_USECREMAIN(process_seeds_usec));
}
#endif

static void gt_diagbandseed_process_seeds_stat(FILE *stream,
                                               GtTimer *timer,
                                               const GtDiagbandseedCounts
                                                 *process_seeds_counts,
                                               GtUword mlistlen,
                                               GtUword allseqpairs,
                                               bool extendgreedy)
{
  GtWord total_process_seeds_usec
    = gt_timer_elapsed_usec(timer);

  fprintf(stream,"# total number of seeds: " GT_WU
                 "; " GT_WU ", i.e. %.2f%% of all seeds were selected"
                 "; " GT_WU ", i.e. %.2f%% of all seeds were extended)\n",
          mlistlen,
          process_seeds_counts->selected_seeds,
          100.0 * (double) process_seeds_counts->selected_seeds/mlistlen,
          process_seeds_counts->extended_seeds,
          100.0 * (double) process_seeds_counts->extended_seeds/mlistlen);
  fprintf(stream, "# number of matches output: " GT_WU "\n",
          process_seeds_counts->countmatches);
  fprintf(stream, "# number of unsuccessful extensions: " GT_WU "\n",
          process_seeds_counts->failedmatches);
  fprintf(stream, "# sequence pairs with selected seeds: " GT_WU
                  " (%.2f%% of all " GT_WU ")\n",
          process_seeds_counts->seqpairs_with_minsegment,
          100.0 * (double)
                  process_seeds_counts->seqpairs_with_minsegment/allseqpairs,
          allseqpairs);
#ifndef _WIN32
  gt_diagbandseed_process_seeds_times(
                 stream,
                 extendgreedy,
                 mlistlen,
                 process_seeds_counts->extended_seeds,
                 total_process_seeds_usec,
                 process_seeds_counts->total_extension_time_usec);
#else
  fprintf(stream, "# ... processed " GT_WU " seeds ",mlistlen);
  gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT,stream);
#endif
}

static void gt_diagbandseed_set_sequence(GtSeqorEncseq *seqorencseq,
                                         const GtSequencePartsInfo *seqranges,
                                         const GtUchar *bytesequence,
                                         GtUword first_seqstartpos,
                                         GtUword seqnum,
                                         const GtUchar *characters,
                                         GtUchar wildcardshow,
                                         bool haswildcards,
                                         const GtEncseq *encseq_for_seq_desc)
{
  const GtUword
    seqstartpos = gt_sequence_parts_info_seqstartpos(seqranges,seqnum),
    seqendpos = gt_sequence_parts_info_seqendpos(seqranges,seqnum),
    b_off = seqstartpos - first_seqstartpos;
  const char *seqdescptr;
  if (encseq_for_seq_desc != NULL)
  {
    GtUword desclen;
    seqdescptr = gt_encseq_description(encseq_for_seq_desc,
                                       &desclen,
                                       seqnum);
  } else
  {
    seqdescptr = "Unknown";
  }
  GT_SEQORENCSEQ_INIT_SEQ(seqorencseq,bytesequence + b_off,seqdescptr,
                          seqendpos - seqstartpos + 1,characters,
                          wildcardshow,
                          haswildcards);
}

typedef struct
{
  bool b_differs_from_a, a_haswildcards, b_haswildcards;
  const GtUchar *characters;
  GtUchar wildcardshow;
  GtSeqorEncseq aseqorencseq, bseqorencseq;
  GtUchar *a_byte_sequence, *b_byte_sequence;
  GtUword previous_aseqnum,
          a_first_seqnum,
          a_first_seqstartpos,
          b_first_seqnum,
          b_first_seqstartpos;
  const GtEncseq *a_encseq_for_seq_desc, *b_encseq_for_seq_desc;
} GtDiagbandSeedPlainSequence;

static void gt_diagband_seed_plainsequence_init(GtDiagbandSeedPlainSequence *ps,
                                          bool seq_desc_display,
                                          const GtEncseq *aencseq,
                                          const GtSequencePartsInfo *aseqranges,
                                          GtUword aidx,
                                          bool with_a_bytestring,
                                          const GtEncseq *bencseq,
                                          const GtSequencePartsInfo *bseqranges,
                                          GtUword bidx,
                                          bool with_b_bytestring)
{
  ps->previous_aseqnum = GT_UWORD_MAX;
  if (seq_desc_display && aencseq != NULL &&
      gt_encseq_has_description_support(aencseq))
  {
    ps->a_encseq_for_seq_desc = aencseq;
  } else
  {
    ps->a_encseq_for_seq_desc = NULL;
  }
  if (seq_desc_display && bencseq != NULL &&
      gt_encseq_has_description_support(bencseq))
  {
    ps->b_encseq_for_seq_desc = bencseq;
  } else
  {
    ps->b_encseq_for_seq_desc = NULL;
  }
  if (with_a_bytestring)
  {
    ps->a_first_seqnum = gt_sequence_parts_info_start_get(aseqranges,aidx);
    ps->a_first_seqstartpos
      = gt_sequence_parts_info_seqstartpos(aseqranges,ps->a_first_seqnum),
    ps->a_byte_sequence = gt_sequence_parts_info_seq_extract(aencseq,aseqranges,
                                                             aidx);
    ps->a_haswildcards = (gt_encseq_wildcards(aencseq) > 0) ? true : false;
  } else
  {
    ps->a_byte_sequence = NULL;
    GT_SEQORENCSEQ_INIT_ENCSEQ(&ps->aseqorencseq,aencseq);
    ps->a_haswildcards = true;
  }
  if (with_b_bytestring)
  {
    ps->b_first_seqnum = gt_sequence_parts_info_start_get(bseqranges,bidx);
    ps->b_first_seqstartpos
      = gt_sequence_parts_info_seqstartpos(bseqranges,ps->b_first_seqnum);
    if (with_a_bytestring && aencseq == bencseq && aidx == bidx)
    {
      ps->b_differs_from_a = false;
      ps->b_byte_sequence = ps->a_byte_sequence;
      ps->b_haswildcards = ps->a_haswildcards;
    } else
    {
      ps->b_differs_from_a = true;
      ps->b_byte_sequence = gt_sequence_parts_info_seq_extract(bencseq,
                                                               bseqranges,bidx);
      ps->b_haswildcards = (gt_encseq_wildcards(bencseq) > 0) ? true : false;
    }
  } else
  {
    ps->b_byte_sequence = NULL;
    GT_SEQORENCSEQ_INIT_ENCSEQ(&ps->bseqorencseq,bencseq);
    ps->b_haswildcards = true;
  }
  if (with_a_bytestring || with_b_bytestring)
  {
    ps->characters = gt_encseq_alphabetcharacters(aencseq);
    ps->wildcardshow = gt_alphabet_wildcard_show(gt_encseq_alphabet(aencseq));
  }
}

static void gt_diagband_seed_plainsequence_next_segment(
                         GtDiagbandSeedPlainSequence *ps,
                         const GtSequencePartsInfo *aseqranges,
                         GtUword currsegm_aseqnum,
                         const GtSequencePartsInfo *bseqranges,
                         GtUword currsegm_bseqnum)
{
  if (ps->previous_aseqnum == GT_UWORD_MAX ||
      ps->previous_aseqnum < currsegm_aseqnum)
  {
    if (ps->a_byte_sequence != NULL)
    {
      gt_diagbandseed_set_sequence(&ps->aseqorencseq,
                                   aseqranges,
                                   ps->a_byte_sequence,
                                   ps->a_first_seqstartpos,
                                   currsegm_aseqnum,
                                   ps->characters,
                                   ps->wildcardshow,
                                   ps->a_haswildcards,
                                   ps->a_encseq_for_seq_desc);
      ps->previous_aseqnum = currsegm_aseqnum;
    } else
    {
      GtUword seqstartpos = gt_sequence_parts_info_seqstartpos(aseqranges,
                                                             currsegm_aseqnum),
              seqendpos = gt_sequence_parts_info_seqendpos(aseqranges,
                                                           currsegm_aseqnum);

      GT_SEQORENCSEQ_ADD_SEQ_COORDS(&ps->aseqorencseq,
                                    seqstartpos,
                                    seqendpos - seqstartpos + 1);
    }
  } else
  {
    gt_assert(ps->previous_aseqnum == currsegm_aseqnum);
  }
  if (ps->b_byte_sequence != NULL)
  {
    gt_diagbandseed_set_sequence(&ps->bseqorencseq,
                                 bseqranges,
                                 ps->b_byte_sequence,
                                 ps->b_first_seqstartpos,
                                 currsegm_bseqnum,
                                 ps->characters,
                                 ps->wildcardshow,
                                 ps->b_haswildcards,
                                 ps->b_encseq_for_seq_desc);
  } else
  {
    GtUword seqstartpos = gt_sequence_parts_info_seqstartpos(bseqranges,
                                                             currsegm_bseqnum),
            seqendpos = gt_sequence_parts_info_seqendpos(bseqranges,
                                                         currsegm_bseqnum);

    GT_SEQORENCSEQ_ADD_SEQ_COORDS(&ps->bseqorencseq,
                                  seqstartpos,
                                  seqendpos - seqstartpos + 1);
  }
}

static void gt_diagband_seed_plainsequence_delete(
                         GtDiagbandSeedPlainSequence *ps)
{
  if (ps->b_byte_sequence != NULL && ps->b_differs_from_a)
  {
    gt_free(ps->b_byte_sequence);
  }
  if (ps->a_byte_sequence != NULL)
  {
    gt_free(ps->a_byte_sequence);
  }
}

/* start seed extension for seeds in mlist */
static void gt_diagbandseed_process_seeds(GtSeedpairlist *seedpairlist,
                                         const GtDiagbandseedExtendParams *extp,
                                          void *processinfo,
                                          GtQuerymatchoutoptions *querymoutopt,
                                          const GtEncseq *aencseq,
                                          const GtSequencePartsInfo *aseqranges,
                                          GtUword aidx,
                                          const GtEncseq *bencseq,
                                          const GtSequencePartsInfo *bseqranges,
                                          GtUword bidx,
                                          unsigned int seedlength,
                                          GtReadmode query_readmode,
                                          bool verbose,
                                          FILE *stream)
{
  GtDiagbandseedScore *diagband_score;
  GtDiagbandseedPosition *diagband_lastpos;
  GtExtendRelativeCoordsFunc extend_relative_coords_function = NULL;
  GtProcessinfo_and_querymatchspaceptr info_querymatch = {NULL,NULL,NULL};
  const bool same_encseq = (aencseq == bencseq) ? true : false;
  /* Although the sequences of the parts processed are shorter, we need to
     set amaxlen and bmaxlen to the maximum size of all sequences
     to get the same division into diagonal bands for all parts and thus
     obtain results independent of the number of parts chosen. */
  const GtUword amaxlen = gt_encseq_max_seq_length(aencseq),
                bmaxlen = gt_encseq_max_seq_length(bencseq),
                mlistlen = gt_seedpairlist_length(seedpairlist),
                ndiags = 1 + ((amaxlen + bmaxlen) >> extp->logdiagbandwidth),
                minsegmentlen = (extp->mincoverage - 1) / seedlength + 1;
  GtUword diagbands_used;
  GtTimer *timer = NULL;
  GtDiagbandseedCounts process_seeds_counts = {0,0,0,0,0};
  GtDiagbandSeedPlainSequence plainsequence_info;

  process_seeds_counts.withtiming = verbose;
  process_seeds_counts.total_extension_time_usec = 0;
  gt_assert(extp->mincoverage >= seedlength && minsegmentlen >= 1);
  if (mlistlen == 0 || mlistlen < minsegmentlen) {
    return;
  }
  /* select extension method */
  if (extp->extendgreedy) {
    extend_relative_coords_function = gt_greedy_extend_seed_relative;
  } else if (extp->extendxdrop) {
    extend_relative_coords_function = gt_xdrop_extend_seed_relative;
  } else { /* no seed extension */
    return;
  }
  if (verbose) {
    timer = gt_timer_new();
    gt_timer_start(timer);
  }
  gt_diagband_seed_plainsequence_init(&plainsequence_info,
                                      gt_querymatch_seq_desc_display(
                                           extp->display_flag),
                                      aencseq,
                                      aseqranges,
                                      aidx,
                                      extp->a_extend_char_access ==
                                        GT_EXTEND_CHAR_ACCESS_DIRECT ? true
                                                                     : false,
                                      bencseq,
                                      bseqranges,
                                      bidx,
                                      extp->b_extend_char_access ==
                                        GT_EXTEND_CHAR_ACCESS_DIRECT ? true
                                                                     : false
                                      );
  if (verbose)
  {
    if (plainsequence_info.a_byte_sequence != NULL ||
        plainsequence_info.b_byte_sequence != NULL)
    {
      fprintf(stream, "# ... extracted sequences ");
      gt_timer_show_formatted(timer, GT_DIAGBANDSEED_FMT, stream);
      gt_timer_start(timer);
    }
    gt_diagbandseed_match_header(stream,extp,seedlength,ndiags,minsegmentlen);
  }
  gt_diagbandseed_info_qm_set(&info_querymatch,
                              extp,
                              aencseq,
                              querymoutopt,
                              query_readmode,
                              stream,
                              processinfo);

  /* diagband_score[0] and diagband_score[ndiags+1] remain zero as boundaries */
  diagband_score = gt_calloc(ndiags + 2, sizeof *diagband_score);
  diagband_score++; /* so we need not increment the index when
                       accessing diagband_score */
  diagband_lastpos = gt_calloc(ndiags, sizeof *diagband_lastpos);
  if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_STRUCT)
  {
    const GtDiagbandseedSeedPair
      *mlist = gt_seedpairlist_mlist_struct(seedpairlist),
      *mlistend = mlist + mlistlen,
      *last_segment_start = mlistend - minsegmentlen,
      *nextsegm = mlist;

    gt_assert(nextsegm <= last_segment_start);
    /* iterate through segments of equal k-mers */
    while (nextsegm <= last_segment_start)
    {
      GtSeedpairPositions *spp_ptr, *segment_positions;
      const GtDiagbandseedSeedPair *currsegm = nextsegm;
      const GtDiagbandseedSeqnum currsegm_aseqnum = currsegm->aseqnum;
      const GtDiagbandseedSeqnum currsegm_bseqnum = currsegm->bseqnum;

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
      spp_ptr = segment_positions = (GtSeedpairPositions *) currsegm;
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
        spp_ptr->apos = nextsegm->apos;
        spp_ptr->bpos = nextsegm->bpos;
        spp_ptr++;
        nextsegm++;
      } while (nextsegm < mlistend &&
               currsegm_aseqnum == nextsegm->aseqnum &&
               currsegm_bseqnum == nextsegm->bseqnum);

      gt_diagband_seed_plainsequence_next_segment(&plainsequence_info,
                                                  aseqranges,
                                                  currsegm_aseqnum,
                                                  bseqranges,
                                                  currsegm_bseqnum);

      /* from here on we only need the apos and bpos values of the segment, as
         the segment boundaries have been identified.
         second scan: test for mincoverage and overlap to previous extension,
         based on apos and bpos values. */
      gt_diagbandseed_process_segment(extp,
                                      &plainsequence_info.aseqorencseq,
                                      &plainsequence_info.bseqorencseq,
                                      same_encseq,
                                      amaxlen,
                                      seedlength,
                                      diagbands_used,
                                      ndiags,
                                      diagband_score,
                                      diagband_lastpos,
                                      currsegm_aseqnum,
                                      currsegm_bseqnum,
                                      segment_positions,
                                      (GtUword) (nextsegm - currsegm),
                                      &info_querymatch,
                                      query_readmode,
                                      extend_relative_coords_function,
                                      &process_seeds_counts);
    }
  } else
  {
    if (seedpairlist->splt == GT_DIAGBANDSEED_SPLT_ULONG)
    {
      const GtUword *mlist = gt_seedpairlist_mlist_ulong(seedpairlist),
                    *mlistend = mlist + mlistlen,
                    *last_segment_start = mlistend - minsegmentlen,
                    *nextsegm = mlist;
     GtUword nextsegm_a_bseqnum;

      gt_assert(nextsegm <= last_segment_start);
      /* iterate through segments of equal k-mers */
      nextsegm_a_bseqnum
        = gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*mlist);
      while (nextsegm <= last_segment_start)
      {
        GtSeedpairPositions *spp_ptr, *segment_positions;
        GtUword currsegm_aseqnum, currsegm_bseqnum;
        const GtUword *currsegm = nextsegm;
        const GtUword currsegm_a_bseqnum = nextsegm_a_bseqnum;

        /* if insuffienct number of kmers in segment: skip whole segment */
        if (currsegm_a_bseqnum !=
            gt_seedpairlist_a_bseqnum_ulong (seedpairlist,
                                             currsegm[minsegmentlen-1]))
        {
          do
          {
            nextsegm++;
          } while (nextsegm < mlistend &&
                   currsegm_a_bseqnum ==
                   (nextsegm_a_bseqnum =
                   gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*nextsegm)));
          continue; /* process next segment */
        }

        /* this segment begining with nextsegm pssibly has enough seeds */
        diagbands_used = 0;
        currsegm_aseqnum = gt_seedpairlist_extract_ulong(seedpairlist,*currsegm,
                                                         idx_aseqnum) +
                           seedpairlist->aseqrange_start;
        currsegm_bseqnum = gt_seedpairlist_extract_ulong(seedpairlist,*currsegm,
                                                         idx_bseqnum) +
                           seedpairlist->bseqrange_start;
        spp_ptr = segment_positions = (GtSeedpairPositions *) currsegm;
        do
        {
          GtUword apos = gt_seedpairlist_extract_ulong(seedpairlist,*nextsegm,
                                                       idx_apos),
                  diag;

          spp_ptr->bpos
            = gt_seedpairlist_extract_ulong(seedpairlist,*nextsegm,idx_bpos);
          spp_ptr->apos = apos;
          diag = GT_DIAGBANDSEED_DIAG(amaxlen, spp_ptr->apos, spp_ptr->bpos);

          gt_assert(diag < ndiags);
          diagbands_used += gt_diagbandseed_update_dband(diagband_score,
                                                         diagband_lastpos,
                                                         diag,
                                                         spp_ptr->bpos,
                                                         seedlength);
          spp_ptr++;
          nextsegm++;
        } while (nextsegm < mlistend &&
                 currsegm_a_bseqnum ==
                 (nextsegm_a_bseqnum =
                 gt_seedpairlist_a_bseqnum_ulong (seedpairlist,*nextsegm)));

        gt_diagband_seed_plainsequence_next_segment(&plainsequence_info,
                                                    aseqranges,
                                                    currsegm_aseqnum,
                                                    bseqranges,
                                                    currsegm_bseqnum);

        /* from here on we only need the apos and bpos values of the segment, as
           the segment boundaries have been identified.
           second scan: test for mincoverage and overlap to previous extension,
           based on apos and bpos values. */
        gt_diagbandseed_process_segment(extp,
                                        &plainsequence_info.aseqorencseq,
                                        &plainsequence_info.bseqorencseq,
                                        same_encseq,
                                        amaxlen,
                                        seedlength,
                                        diagbands_used,
                                        ndiags,
                                        diagband_score,
                                        diagband_lastpos,
                                        currsegm_aseqnum,
                                        currsegm_bseqnum,
                                        segment_positions,
                                        (GtUword) (nextsegm - currsegm),
                                        &info_querymatch,
                                        query_readmode,
                                        extend_relative_coords_function,
                                        &process_seeds_counts);
      }
    } else
    {
      GtDiagbandseedSeedPair nextsegment;
      const GtUword minsegmentlen_offset = (minsegmentlen - 1) *
                                           seedpairlist->bytes_seedpair,
                    last_segment_offset = (mlistlen - minsegmentlen) *
                                          seedpairlist->bytes_seedpair,
                    mlistlen_offset = mlistlen * seedpairlist->bytes_seedpair;
      GtUword nextsegment_offset = 0;

      gt_assert(seedpairlist->splt == GT_DIAGBANDSEED_SPLT_BYTESTRING);
      /* iterate through segments of equal k-mers, segment has length > 0 */
      gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,0);
      while (nextsegment_offset <= last_segment_offset)
      {
        GtSeedpairPositions *spp_ptr, *segment_positions;
        GtDiagbandseedSeedPair endminsegment;
        GtDiagbandseedSeqnum currsegm_aseqnum = nextsegment.aseqnum;
        GtDiagbandseedSeqnum currsegm_bseqnum = nextsegment.bseqnum;

        gt_diagbandseed_decode_seedpair(&endminsegment,seedpairlist,
                                        nextsegment_offset +
                                        minsegmentlen_offset);
        if (currsegm_aseqnum != endminsegment.aseqnum ||
            currsegm_bseqnum != endminsegment.bseqnum)
        {
          /* insuffienct number of kmers in segment: skip whole segment */
          while (true)
          {
            nextsegment_offset += seedpairlist->bytes_seedpair;
            if (nextsegment_offset >= mlistlen_offset)
            {
              break;
            }
            gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,
                                            nextsegment_offset);
            if (currsegm_aseqnum != nextsegment.aseqnum ||
                currsegm_bseqnum != nextsegment.bseqnum)
            {
              break;
            }
          }
          continue; /* process next segment */
        }

        diagbands_used = 0;
        spp_ptr = segment_positions
                = (GtSeedpairPositions *)
                  (gt_seedpairlist_mlist_bytestring(seedpairlist) +
                   nextsegment_offset);
        while (true)
        {
          const GtUword diag = GT_DIAGBANDSEED_DIAG(amaxlen, nextsegment.apos,
                                                    nextsegment.bpos);

          gt_assert(diag < ndiags);
          diagbands_used += gt_diagbandseed_update_dband(diagband_score,
                                                         diagband_lastpos,
                                                         diag,
                                                         nextsegment.bpos,
                                                         seedlength);
          spp_ptr->apos = nextsegment.apos;
          spp_ptr->bpos = nextsegment.bpos;
          spp_ptr++;
          nextsegment_offset += seedpairlist->bytes_seedpair;
          if (nextsegment_offset >= mlistlen_offset)
          {
            break;
          }
          gt_diagbandseed_decode_seedpair(&nextsegment,seedpairlist,
                                          nextsegment_offset);
          if (currsegm_aseqnum != nextsegment.aseqnum ||
              currsegm_bseqnum != nextsegment.bseqnum)
          {
            break;
          }
        }

        /* from here on we only need the apos and bpos values of the segment, as
           the segment boundaries have been identified.
           second scan: test for mincoverage and overlap to previous extension,
           based on apos and bpos values. */
        currsegm_aseqnum += seedpairlist->aseqrange_start;
        currsegm_bseqnum += seedpairlist->bseqrange_start;
        gt_diagband_seed_plainsequence_next_segment(&plainsequence_info,
                                                    aseqranges,
                                                    currsegm_aseqnum,
                                                    bseqranges,
                                                    currsegm_bseqnum);
        gt_diagbandseed_process_segment(extp,
                                        &plainsequence_info.aseqorencseq,
                                        &plainsequence_info.bseqorencseq,
                                        same_encseq,
                                        amaxlen,
                                        seedlength,
                                        diagbands_used,
                                        ndiags,
                                        diagband_score,
                                        diagband_lastpos,
                                        currsegm_aseqnum,
                                        currsegm_bseqnum,
                                        segment_positions,
                                        (GtUword) (spp_ptr - segment_positions),
                                        &info_querymatch,
                                        query_readmode,
                                        extend_relative_coords_function,
                                        &process_seeds_counts);
      }
    }
  }
  diagband_score--; /* need to recover original base adress */
  gt_free(diagband_score);
  gt_free(diagband_lastpos);
  gt_querymatch_delete(info_querymatch.querymatchspaceptr);
  gt_karlin_altschul_stat_delete(info_querymatch.karlin_altschul_stat);
  gt_diagband_seed_plainsequence_delete(&plainsequence_info);
  if (verbose)
  {
    const GtUword allseqpairs = (seedpairlist->aseqrange_end -
                                 seedpairlist->aseqrange_start + 1) *
                                (seedpairlist->bseqrange_end -
                                 seedpairlist->bseqrange_start + 1);
    gt_diagbandseed_process_seeds_stat(stream,
                                       timer,
                                       &process_seeds_counts,
                                       mlistlen,
                                       allseqpairs,
                                       extp->extendgreedy);
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
                                     const GtEncseq *aencseq,
                                     const GtSequencePartsInfo *aseqranges,
                                     GtUword aidx,
                                     const GtEncseq *bencseq,
                                     const GtSequencePartsInfo *bseqranges,
                                     GtUword bidx,
                                     GtFtTrimstat *trimstat,
                                     GtError *err)
{
  GtArrayGtDiagbandseedKmerPos blist;
  GtSeedpairlist *seedpairlist = NULL;
  GtDiagbandseedKmerIterator *aiter = NULL, *biter = NULL;
  GtUword alen = 0, blen = 0, mlistlen = 0, maxfreq, len_used, alignmentwidth;
  GtRange seedpairdistance = *arg->seedpairdistance;
  char *blist_file = NULL;
  int had_err = 0;
  bool alist_blist_id, both_strands, selfcomp, equalranges, use_blist = false;
  size_t sizeofunit;
  const GtDiagbandseedExtendParams *extp = NULL;
  GtFtPolishing_info *pol_info = NULL;
  void *processinfo = NULL;
  GtQuerymatchoutoptions *querymoutopt = NULL;
  const GtUword anumseqranges = gt_sequence_parts_info_number(aseqranges),
                bnumseqranges = gt_sequence_parts_info_number(bseqranges);

  gt_assert(arg != NULL);
  maxfreq = arg->maxfreq;
  selfcomp = (arg->bencseq == arg->aencseq &&
              gt_sequence_parts_info_overlap(aseqranges,aidx,bseqranges,bidx))
              ? true : false;
  equalranges = gt_sequence_parts_info_equal(aseqranges,aidx,bseqranges,bidx);
  alist_blist_id = (selfcomp && !arg->nofwd && equalranges) ? true : false;
  both_strands = (arg->norev || arg->nofwd) ? false : true;
  if (!alist_blist_id) {
    seedpairdistance.start = 0UL;
  }

  if (arg->verbose && (anumseqranges > 1 || bnumseqranges > 1))
  {
    fprintf(stream, "# process part " GT_WU " (sequences " GT_WU "..." GT_WU
                    ") vs part " GT_WU " (sequences " GT_WU "..." GT_WU ")\n",
            aidx + 1,
            gt_sequence_parts_info_start_get(aseqranges,aidx),
            gt_sequence_parts_info_end_get(aseqranges,aidx),
            bidx + 1,
            gt_sequence_parts_info_start_get(bseqranges,bidx),
            gt_sequence_parts_info_end_get(bseqranges,bidx));
  }

  /* Create k-mer iterator for alist */
  if (alist == NULL) {
    char *alist_file;
    alist_file = gt_diagbandseed_kmer_filename(arg->aencseq, arg->seedlength,
                                               true, anumseqranges,aidx);
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
                                               !arg->nofwd, bnumseqranges,bidx);
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
    const GtReadmode readmode_kmerscan = arg->nofwd ? GT_READMODE_COMPL
                                                    : GT_READMODE_FORWARD;
    const GtUword known_size = (selfcomp && equalranges) ? alen : 0;
    blist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->seedlength,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
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
  seedpairlist = gt_seedpairlist_new(arg->splt,aseqranges,aidx,bseqranges,bidx);
  sizeofunit = gt_seedpairlist_sizeofunit(seedpairlist);
  had_err = gt_diagbandseed_get_mlistlen_maxfreq(&mlistlen,
                                                 &maxfreq,
                                                 aiter,
                                                 biter,
                                                 arg->memlimit,
                                                 sizeofunit,
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
  } else
  {
    gt_seedpairlist_delete(seedpairlist);
    seedpairlist = NULL;
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
                                               extp->a_extend_char_access,
                                               extp->b_extend_char_access,
                                               extp->cam_generic,
                                               extp->sensitivity,
                                               pol_info);
    if (extp->benchmark) {
      gt_greedy_extend_matchinfo_silent_set(grextinfo);
    }
    if (trimstat != NULL)
    {
      gt_greedy_extend_matchinfo_trimstat_set(grextinfo,trimstat);
    }
    processinfo = (void *) grextinfo;
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
    processinfo = (void *) xdropinfo;
  }
  alignmentwidth = gt_querymatch_display_alignmentwidth(extp->display_flag);
  if (extp->extendxdrop || alignmentwidth > 0 || extp->verify_alignment)
  {
    querymoutopt = gt_querymatchoutoptions_new(true,
                                               false,
                                               alignmentwidth,
                                               NULL,
                                               NULL);
    gt_assert(querymoutopt != NULL);
    if (extp->extendxdrop || extp->extendgreedy) {
      const GtUword sensitivity = extp->extendxdrop ? 100UL : extp->sensitivity;
      gt_querymatchoutoptions_extend(querymoutopt,
                                     extp->errorpercentage,
                                     extp->maxalignedlendifference,
                                     extp->history_size,
                                     extp->perc_mat_history,
                                     extp->a_extend_char_access,
                                     extp->b_extend_char_access,
                                     extp->cam_generic,
                                     extp->weakends,
                                     sensitivity,
                                     extp->matchscore_bias,
                                     extp->always_polished_ends,
                                     extp->display_flag);
    }
  }
  /* process first mlist */
  gt_assert(seedpairlist != NULL);
  gt_diagbandseed_process_seeds(seedpairlist,
                                arg->extp,
                                processinfo,
                                querymoutopt,
                                aencseq,aseqranges,aidx,
                                bencseq,bseqranges,bidx,
                                arg->seedlength,
                                arg->nofwd ? GT_READMODE_REVCOMPL
                                           : GT_READMODE_FORWARD,
                                arg->verbose,
                                stream);
  gt_seedpairlist_reset(seedpairlist);
  gt_querymatchoutoptions_reset(querymoutopt);

  /* Third (reverse) k-mer list */
  if (both_strands) {
    GtUword mrevlen = 0;
    GtArrayGtDiagbandseedKmerPos clist;

    gt_assert(blist_file == NULL && !use_blist);
    seedpairdistance.start = 0UL;
    if (arg->use_kmerfile) {
      blist_file = gt_diagbandseed_kmer_filename(arg->bencseq, arg->seedlength,
                                                 false, bnumseqranges,
                                                 bidx);
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
      const GtReadmode readmode_kmerscan = GT_READMODE_COMPL;
      clist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->seedlength,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
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
                                                     sizeofunit,
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
      gt_diagbandseed_get_seedpairs(seedpairlist,
                                    aiter,
                                    biter,
                                    maxfreq,
                                    mrevlen,
                                    &seedpairdistance,
                                    selfcomp,
                                    arg->debug_seedpair,
                                    arg->verbose,
                                    stream);
      mrevlen = gt_seedpairlist_length(seedpairlist);
      if (arg->verify && mrevlen > 0) {
        had_err = gt_diagbandseed_verify(seedpairlist,
                                         arg->aencseq,
                                         arg->bencseq,
                                         arg->seedlength,
                                         true,
                                         arg->verbose,
                                         stream,
                                         err);
        if (had_err) {
          gt_seedpairlist_delete(seedpairlist);
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
    gt_diagbandseed_process_seeds(seedpairlist,
                                  arg->extp,
                                  processinfo,
                                  querymoutopt,
                                  aencseq,aseqranges,aidx,
                                  bencseq,bseqranges,bidx,
                                  arg->seedlength,
                                  GT_READMODE_REVCOMPL,
                                  arg->verbose,
                                  stream);
  }
  gt_seedpairlist_delete(seedpairlist);

  /* Clean up */
  if (extp->extendgreedy)
  {
    polishing_info_delete(pol_info);
    gt_greedy_extend_matchinfo_delete((GtGreedyextendmatchinfo *) processinfo);
  } else
  {
    if (extp->extendxdrop)
    {
      gt_xdrop_matchinfo_delete((GtXdropmatchinfo *) processinfo);
    }
  }
  gt_querymatchoutoptions_delete(querymoutopt);
  return had_err;
}

#ifdef GT_THREADS_ENABLED
typedef struct{
  const GtDiagbandseedInfo *arg;
  const GtArrayGtDiagbandseedKmerPos *alist;
  FILE *stream;
  const GtEncseq *aencseq, *bencseq;
  const GtSequencePartsInfo *aseqranges,
                            *bseqranges;
  GtArray *combinations;
  int had_err;
  GtError *err;
} GtDiagbandseedThreadInfo;

static void gt_diagbandseed_thread_info_set(GtDiagbandseedThreadInfo *ti,
                                     const GtDiagbandseedInfo *arg,
                                     const GtArrayGtDiagbandseedKmerPos *alist,
                                     FILE *stream,
                                     const GtEncseq *aencseq,
                                     const GtSequencePartsInfo *aseqranges,
                                     const GtEncseq *bencseq,
                                     const GtSequencePartsInfo *bseqranges,
                                     GtArray *combinations,
                                     GtError *err)
{
  gt_assert(ti != NULL);
  ti->arg = arg;
  ti->alist = alist;
  ti->stream = stream;
  ti->aencseq = aencseq;
  ti->aseqranges = aseqranges;
  ti->bencseq = bencseq;
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
      info->had_err = gt_diagbandseed_algorithm(
                           info->arg,
                           info->alist,
                           info->stream,
                           info->aencseq,info->aseqranges,comb->a,
                           info->bencseq,info->bseqranges,comb->b,
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
    printf("# write " GT_WU " %u-mers to file %s\n",
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
                        const GtSequencePartsInfo *aseqranges,
                        const GtSequencePartsInfo *bseqranges,
                        const GtUwordPair *pick,
                        GtError *err)
{
  const bool self = arg->aencseq == arg->bencseq ? true : false;
  const bool apick = pick->a != GT_UWORD_MAX ? true : false;
  const bool bpick = pick->b != GT_UWORD_MAX ? true : false;
  GtArrayGtDiagbandseedKmerPos alist;
  GtUword aidx, bidx;
  int had_err = 0;
  const GtUword anumseqranges = gt_sequence_parts_info_number(aseqranges),
                bnumseqranges = gt_sequence_parts_info_number(bseqranges);
  GtFtTrimstat *trimstat = NULL;
#ifdef GT_THREADS_ENABLED
  GtDiagbandseedThreadInfo *tinfo = gt_malloc(gt_jobs * sizeof *tinfo);
  FILE **stream;
  unsigned int tidx;

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

      for (bidx = 0; !had_err && bidx < bnumseqranges; bidx++) {
        char *path;
        if (bpick && pick->b != bidx) continue;

        path = gt_diagbandseed_kmer_filename(arg->bencseq, arg->seedlength, fwd,
                                             bnumseqranges, bidx);
        if (!gt_file_exists(path)) {
          GtArrayGtDiagbandseedKmerPos blist;
          GtReadmode readmode_kmerscan = fwd ? GT_READMODE_FORWARD
                                             : GT_READMODE_COMPL;

          blist = gt_diagbandseed_get_kmers(
                              arg->bencseq,
                              arg->seedlength,
                              readmode_kmerscan,
                              gt_sequence_parts_info_start_get(bseqranges,bidx),
                              gt_sequence_parts_info_end_get(bseqranges,bidx),
                              arg->debug_kmer,
                              arg->verbose,
                              0,
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
  for (aidx = 0; !had_err && aidx < anumseqranges; aidx++) {
    /* create alist here to prevent redundant calculations */
    char *path = NULL;
    bool use_alist = false;
    if (apick && pick->a != aidx) continue;

    if (arg->use_kmerfile) {
      path = gt_diagbandseed_kmer_filename(arg->aencseq, arg->seedlength, true,
                                           anumseqranges, aidx);
    }

    if (!arg->use_kmerfile || !gt_file_exists(path)) {
      use_alist = true;
      alist = gt_diagbandseed_get_kmers(
                              arg->aencseq,
                              arg->seedlength,
                              GT_READMODE_FORWARD,
                              gt_sequence_parts_info_start_get(aseqranges,aidx),
                              gt_sequence_parts_info_end_get(aseqranges,aidx),
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
      while (!had_err && bidx < bnumseqranges) {
        if (!bpick || pick->b == bidx) {
          /* start algorithm with chosen sequence ranges */
          had_err = gt_diagbandseed_algorithm(
                           arg,
                           use_alist ? &alist : NULL,
                           stdout,
                           arg->aencseq,aseqranges,aidx,
                           arg->bencseq,bseqranges,bidx,
                           trimstat,
                           err);
        }
        bidx++;
      }
#ifdef GT_THREADS_ENABLED
    } else if (!arg->use_kmerfile) {
      const GtUword num_runs = bpick ? 1 : bnumseqranges - bidx;
      const GtUword num_runs_per_thread = (num_runs - 1) / gt_jobs + 1;
      const GtUword num_threads = (num_runs - 1) / num_runs_per_thread + 1;
      GtArray *combinations = gt_array_new(sizeof (GtUwordPair));
      GtArray *threads = gt_array_new(sizeof (GtThread *));

      gt_assert(bidx < bnumseqranges);
      gt_assert(num_threads <= gt_jobs);
      gt_assert(!bpick || num_threads == 1);

      /* start additional threads */
      for (tidx = 1; !had_err && tidx < num_threads; tidx++) {
        GtThread *thread;
        GtUword idx;
        bidx += num_runs_per_thread;
        const GtUword end = MIN(bidx + num_runs_per_thread, bnumseqranges);

        for (idx = bidx; idx < end; idx++) {
          GtUwordPair comb = {aidx, idx};
          gt_array_add(combinations, comb);
        }
        gt_diagbandseed_thread_info_set(tinfo + tidx,
                                        arg,
                                        use_alist ? &alist : NULL,
                                        stream[tidx],
                                        arg->aencseq,
                                        aseqranges,
                                        arg->bencseq,
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
             idx < MIN(bidx + num_runs_per_thread, bnumseqranges);
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
                                        arg->aencseq,
                                        aseqranges,
                                        arg->bencseq,
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
    for (aidx = 0; aidx < anumseqranges; aidx++) {
      if (apick && pick->a != aidx) continue;
      for (bidx = self ? aidx : 0; bidx < bnumseqranges; bidx++) {
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
                                      arg->aencseq,
                                      aseqranges,
                                      arg->bencseq,
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
                                      arg->aencseq,
                                      aseqranges,
                                      arg->bencseq,
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
  gt_ft_trimstat_delete(trimstat);
  return had_err;
}
