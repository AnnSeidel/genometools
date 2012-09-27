/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/intbits.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "sfx-linlcp.h"
#include "sfx-sain.h"

#define GT_SSTARLENGTH_MAX 50

typedef struct
{
  unsigned long totallength,
                specialcharacters,
                numofchars,
                *bucketsize;
  union
  {
    const GtEncseq *encseq;
    const unsigned long *array;
  } seq;
  bool hasencseq;
} GtSainseq;

static GtSainseq *gt_sain_seq_new_from_encseq(const GtEncseq *encseq)
{
  unsigned long idx;
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->hasencseq = true;
  sainseq->seq.encseq = encseq;
  sainseq->totallength = gt_encseq_total_length(encseq);
  sainseq->specialcharacters = gt_encseq_specialcharacters(encseq);
  gt_assert(sainseq->totallength >= sainseq->specialcharacters);
  sainseq->numofchars = (unsigned long) gt_encseq_alphabetnumofchars(encseq);
  sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) *
                                  sainseq->numofchars);
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = gt_encseq_charcount(encseq,(GtUchar) idx);
  }
  return sainseq;
}

static GtSainseq *gt_sain_seq_new_from_array(unsigned long *arr,
                                             unsigned long len,
                                             unsigned long numofchars)
{
  unsigned long idx;
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->hasencseq = false;
  sainseq->seq.array = arr;
  sainseq->totallength = len;
  sainseq->specialcharacters = 0;
  sainseq->numofchars = numofchars;
  sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) *
                                  sainseq->numofchars);
  for (idx = 0; idx<numofchars; idx++)
  {
    sainseq->bucketsize[idx] = 0;
  }
  for (idx = 0; idx<len; idx++)
  {
    gt_assert(arr[idx] < numofchars);
    sainseq->bucketsize[arr[idx]]++;
  }
  return sainseq;
}

static void gt_sain_seq_delete(GtSainseq *sainseq)
{
  if (sainseq != NULL)
  {
    gt_free(sainseq->bucketsize);
    gt_free(sainseq);
  }
}

static unsigned long gt_sain_seq_getchar(const GtSainseq *sainseq,
                                         unsigned long position,
                                         GT_UNUSED bool assertnospecial)
{
  gt_assert(position < sainseq->totallength);
  if (sainseq->hasencseq)
  {
    GtUchar cc = gt_encseq_get_encoded_char(sainseq->seq.encseq,
                                            position,
                                            GT_READMODE_FORWARD);
    gt_assert(!assertnospecial || ISNOTSPECIAL(cc));
    return ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (unsigned long) cc;
  } else
  {
    return sainseq->seq.array[position];
  }
}

typedef struct
{
  GtBitsequence *isStype;
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
  GtSainseq *sainseq;
} GtSaininfo;

static GtSaininfo *gt_sain_info_new(GtSainseq *sainseq)
{
  unsigned long position,
                idx,
                nextcc,
                nextSstartypepos;
  bool nextisStype = true;
  GtSaininfo *saininfo;

  saininfo = gt_malloc(sizeof *saininfo);
  saininfo->sainseq = sainseq;
  saininfo->totalSstarlength = 0;
  saininfo->countStype = 1UL;
  saininfo->countSstartype = 0;
  saininfo->longerthanmax = 0;
  GT_INITBITTAB(saininfo->isStype,saininfo->sainseq->totallength+1);
  GT_SETIBIT(saininfo->isStype,saininfo->sainseq->totallength);
  nextSstartypepos = saininfo->sainseq->totallength;
  nextcc = GT_UNIQUEINT(saininfo->sainseq->totallength);
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
  printf("%lu: S\n",saininfo->sainseq->totallength);
#endif
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    saininfo->lendist[idx] = 0;
  }
  for (position = saininfo->sainseq->totallength-1; /* Nothing */; position--)
  {
    bool currentisStype;
    unsigned long currentcc;

    currentcc = gt_sain_seq_getchar(saininfo->sainseq,position,false);
    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      saininfo->countStype++;
      currentisStype = true;
      GT_SETIBIT(saininfo->isStype,position);
#ifdef SAINSHOWSTATE
      printf("%lu: S\n",position);
#endif
    } else
    {
      currentisStype = false;
#ifdef SAINSHOWSTATE
      printf("%lu: L\n",position);
#endif
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      saininfo->countSstartype++;
      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      saininfo->totalSstarlength += currentlen;
#ifdef SAINSHOWSTATE
      printf("Sstar: %lu\n",position+1);
#endif
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        saininfo->lendist[currentlen]++;
      } else
      {
        saininfo->longerthanmax++;
      }
      nextSstartypepos = position + 1;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_assert(GT_MULT2(saininfo->countSstartype) <=
            saininfo->sainseq->totallength);
  return saininfo;
}

static void gt_sain_info_delete(GtSaininfo *saininfo)
{
  if (saininfo != NULL)
  {
    gt_free(saininfo->isStype);
    gt_free(saininfo);
  }
}

static bool gt_sain_info_isSstartype(const GtSaininfo *saininfo,
                                       unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return position == saininfo->sainseq->totallength ||
         (position > 0 &&
         GT_ISIBITSET(saininfo->isStype,position) &&
         !GT_ISIBITSET(saininfo->isStype,position-1)) ? true : false;
}

#ifdef CRITICAL
unsigned long countDcriticalsubstrings (const GtSaininfo *saininfo,
                                        unsigned long d)
{
  unsigned long int i = 0, j = 0;

  gt_assert(d >= 2UL);
  while (i < saininfo->sainseq->totallength)
  {
    unsigned long h;
    bool isLMS = false;

    for (h = 1UL; h <= d; h++)
    {
      if (gt_sain_info_isSstartype(saininfo,i+h))
      {
        isLMS = true;
        break;
      }
    }
    if (j == 0 && !isLMS)
    {
      i += d;
      continue;
    }
    i = isLMS ? i + h : i + d;
    gt_assert(i>0);
    /*printf("crititical %lu\n",i-1);*/
    j++;
  }
  return j;
}
#endif

static void gt_sain_info_show(const GtSaininfo *saininfo)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",saininfo->countStype,
                (double) saininfo->countStype/saininfo->sainseq->totallength);
  printf("Sstar-type: %lu (%.2f)\n",saininfo->countSstartype,
              (double) saininfo->countSstartype/saininfo->sainseq->totallength);
  printf("Sstar-type.length: %lu (%.2f)\n",saininfo->totalSstarlength,
              (double) saininfo->totalSstarlength/saininfo->countSstartype);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (saininfo->lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,saininfo->lendist[idx],
               (double) saininfo->lendist[idx]/saininfo->countSstartype);
    }
  }
  if (saininfo->longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,saininfo->longerthanmax,
                 (double) saininfo->longerthanmax/saininfo->countSstartype);
  }
#ifdef CRITICAL
  {
    unsigned long d;
    for (d=2UL; d<10UL; d++)
    {
      unsigned long critical = countDcriticalsubstrings (saininfo,d);
      printf("d=%lu,critical=%lu (%.2f)\n",d,critical,
                        (double) critical/saininfo->sainseq->totallength);
    }
  }
#endif
}

static bool gt_sain_info_isStype(const GtSaininfo *saininfo,
                                 unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return GT_ISIBITSET(saininfo->isStype,position) ? true : false;
}

static void gt_sain_endbuckets(unsigned long *leftborder,
                               const unsigned long *bucketsize,
                               unsigned long numofchars)
{
  unsigned long charidx;

  leftborder[0] = bucketsize[0];
  for (charidx = 1UL; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx];
  }
}

static void gt_sain_startbuckets(unsigned long *leftborder,
                                 const unsigned long *bucketsize,
                                 unsigned long numofchars)
{
  unsigned long charidx;

  leftborder[0] = 0;
  for (charidx = 1UL; charidx<numofchars; charidx++)
  {
    leftborder[charidx] = leftborder[charidx-1] + bucketsize[charidx-1];
  }
}

static void insertSstarsuffixes(const GtSaininfo *saininfo,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                GT_UNUSED unsigned long nonspecialentries)
{
  unsigned long position;

  for (position = 0; position < saininfo->sainseq->totallength; position++)
  {
    if (gt_sain_info_isSstartype(saininfo,position))
    {
      unsigned long putidx,
                    cc = gt_sain_seq_getchar(saininfo->sainseq,position,true);
      gt_assert(leftborder[cc] > 0);
      putidx = --leftborder[cc];
      gt_assert(putidx < nonspecialentries);
      suftab[putidx] = position;
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[%lu]=%lu\n",putidx,position);
#endif
    }
  }
}

static void induceLtypesuffixes(const GtSaininfo *saininfo,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                unsigned long nonspecialentries)
{
  unsigned long idx;

  for (idx = 0; idx < nonspecialentries; idx++)
  {
    unsigned long position = suftab[idx];

    if (position != ULONG_MAX && position > 0)
    {
      gt_assert(position < saininfo->sainseq->totallength);
      if (!gt_sain_info_isStype(saininfo,position-1))
      {
        unsigned long cc = gt_sain_seq_getchar(saininfo->sainseq,position-1,
                                               false);
        if (cc < saininfo->sainseq->numofchars)
        {
          unsigned long putidx = leftborder[cc]++;
          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position-1;
#ifdef SAINSHOWSTATE
          printf("L-induce: suftab[%lu]=%lu\n",putidx,position-1);
#endif
        }
      }
    }
  }
}

static void induceStypesfromspecialranges(GT_UNUSED const GtSaininfo *saininfo,
                                          const GtEncseq *encseq,
                                          unsigned long *suftab,
                                          unsigned long *leftborder,
                                          GT_UNUSED unsigned long
                                                    nonspecialentries)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    sri = gt_specialrangeiterator_new(encseq,false);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (range.start > 0)
      {
        unsigned long putidx;
        GtUchar cc;

        gt_assert(gt_sain_info_isStype(saininfo,range.start-1));
        cc = gt_encseq_get_encoded_char(encseq,range.start-1,
                                        GT_READMODE_FORWARD);
        gt_assert(ISNOTSPECIAL(cc) && leftborder[cc] > 0);
        putidx = --leftborder[cc];
        gt_assert(putidx < nonspecialentries);
        suftab[putidx] = range.start-1;
#ifdef SAINSHOWSTATE
        printf("Srange-induce: suftab[%lu]=%lu in %d-bucket\n",putidx,
                          range.start-1,(int) cc);
#endif
      }
    }
    gt_specialrangeiterator_delete(sri);
  }
}

static void induceStypesuffixes(const GtSaininfo *saininfo,
                                unsigned long *suftab,
                                unsigned long *leftborder,
                                unsigned long nonspecialentries)
{
  unsigned long idx, cc;

  cc = gt_sain_seq_getchar(saininfo->sainseq,saininfo->sainseq->totallength-1,
                           false);
  if (cc < saininfo->sainseq->numofchars)
  {
    unsigned long putidx = --leftborder[cc];
    gt_assert(putidx < nonspecialentries);
    suftab[putidx] = saininfo->sainseq->totallength-1;
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab[%lu]=%lu\n",putidx,
                                         saininfo->sainseq->totallength-1);
#endif
  }
  if (saininfo->sainseq->hasencseq)
  {
    induceStypesfromspecialranges(saininfo,saininfo->sainseq->seq.encseq,suftab,
                                  leftborder,nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (idx = nonspecialentries - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx];

    if (position != ULONG_MAX && position > 0)
    {
      if (gt_sain_info_isStype(saininfo,position-1))
      {
        cc = gt_sain_seq_getchar(saininfo->sainseq,position-1,false);
        if (cc < saininfo->sainseq->numofchars)
        {
          unsigned long putidx = --leftborder[cc];
          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position-1;
#ifdef SAINSHOWSTATE
          printf("S-induce: suftab[%lu]=%lu\n",putidx,position-1);
#endif
        }
      }
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void moveSstar2front(const GtSaininfo *saininfo,
                            unsigned long *suftab,
                            unsigned long nonspecialentries)
{
  unsigned long ridx, position;

  for (ridx = 0; ridx < nonspecialentries; ridx++)
  {
    position = suftab[ridx];
    if (!gt_sain_info_isSstartype(saininfo,position))
    {
      break;
    }
  }
  if (ridx < saininfo->countSstartype)
  {
    unsigned long widx;

    for (widx = ridx, ridx++; /* Nothing */; ridx++)
    {
      gt_assert(widx < ridx && ridx < nonspecialentries);
      position = suftab[ridx];
      gt_assert(position != ULONG_MAX);
      if (gt_sain_info_isSstartype(saininfo,position))
      {
        suftab[widx++] = position;
        suftab[ridx] = ULONG_MAX;
        if (widx == saininfo->countSstartype)
        {
          break;
        }
      }
    }
  }
}

static int gt_sain_compare_substrings(const GtSaininfo *saininfo,
                                      bool withtype,
                                      unsigned long start1,
                                      unsigned long start2)
{
  bool firstcmp = true, previousisS = true;

  gt_assert(start1 <= saininfo->sainseq->totallength &&
            start2 <= saininfo->sainseq->totallength &&
            start1 != start2);

  while (true)
  {
    unsigned long cc1, cc2;

    if (start1 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sain_seq_getchar(saininfo->sainseq,start1,false);
    cc2 = gt_sain_seq_getchar(saininfo->sainseq,start2,false);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    if (withtype)
    {
      if (gt_sain_info_isStype(saininfo,start1))
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          if (!firstcmp && !previousisS)
          {
            return 0; /* previous is L => Sstar in both stop with equality */
          }
        } else
        {
          /* S > L */
          return 1;
        }
      } else
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          /* L < S */
          return -1;
        } else
        {
          /* L == L */
          previousisS = false;
        }
      }
      firstcmp = false;
    }
    start1++;
    start2++;
  }
}

static void sain_setundefined(unsigned long *suftab,
                              unsigned long start, unsigned long end)
{
  unsigned long idx;

  for (idx = start; idx <= end; idx++)
  {
    suftab[idx] = ULONG_MAX;
  }
}

static unsigned long assignSstarnames(const GtSaininfo *saininfo,
                                      unsigned long *suftab,
                                      unsigned long suftabentries)
{
  unsigned long idx, previouspos, currentname = 0;

  gt_assert(suftabentries > 0);
  sain_setundefined(suftab,saininfo->countSstartype,suftabentries-1);
  previouspos = suftab[0];
  suftab[saininfo->countSstartype + GT_DIV2(previouspos)] = 0;
  gt_assert(gt_sain_info_isSstartype(saininfo,previouspos));
  for (idx = 1UL; idx < saininfo->countSstartype; idx++)
  {
    int cmp;
    unsigned long position = suftab[idx];

    gt_assert(gt_sain_info_isSstartype(saininfo,position));
    cmp = gt_sain_compare_substrings(saininfo,true,previouspos,position);
    gt_assert(cmp != 1);
    if (cmp == -1)
    {
      currentname++;
    }
    gt_assert(saininfo->countSstartype + GT_DIV2(position) < suftabentries);
    suftab[saininfo->countSstartype + GT_DIV2(position)] = currentname;
    previouspos = position;
  }
  return currentname+1;
}

static void movenames2front(const GtSaininfo *saininfo,
                            unsigned long *suftab,
                            GT_UNUSED unsigned long suftabentries)
{
  unsigned long ridx, widx,
                maxridx = saininfo->countSstartype +
                          GT_DIV2(saininfo->sainseq->totallength);
  for (ridx = widx = saininfo->countSstartype; ridx <= maxridx; ridx++)
  {
    if (suftab[ridx] != ULONG_MAX)
    {
      if (widx < ridx)
      {
        gt_assert(widx < suftabentries);
        suftab[widx++] = suftab[ridx];
      } else
      {
        gt_assert(widx == ridx);
        widx++;
      }
    }
  }
  gt_assert(widx == GT_MULT2(saininfo->countSstartype));
}

static void gt_sain_checkorder(const GtSaininfo *saininfo,
                               const unsigned long *suftab,
                               unsigned long start,
                               unsigned long end)
{
  unsigned long idx;

  for (idx = start+1; idx <= end; idx++)
  {
    int cmp
      = gt_sain_compare_substrings(saininfo,false,suftab[idx-1],suftab[idx]);

    gt_assert(cmp == -1);
  }
}

static void insertsortedSstarsuffixes(const GtSaininfo *saininfo,
                                      unsigned long *suftab,
                                      unsigned long *leftborder,
                                      GT_UNUSED
                                      unsigned long nonspecialsuffixes)
{
  unsigned long idx;

  if (saininfo->countSstartype == 0)
  {
    return;
  }
  for (idx = saininfo->countSstartype - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx], putidx, cc;

    cc = gt_sain_seq_getchar(saininfo->sainseq,position,false);
    gt_assert(leftborder[cc] > 0);
    putidx = --leftborder[cc];
    gt_assert(idx <= putidx);
    if (idx < putidx)
    {
      gt_assert(putidx < nonspecialsuffixes);
      suftab[putidx] = position;
      suftab[idx] = ULONG_MAX;
#ifdef SAINSHOWSTATE
      printf("insertsorted: suftab[%lu]=%lu\n",putidx,position);
      printf("insertsorted: suftab[%lu]=undef\n",idx);
#endif
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void gt_sain_filltailsuffixes(unsigned long *suftabtail,
                                     const GtEncseq *encseq)
{
  unsigned long specialcharacters = gt_encseq_specialcharacters(encseq),
                totallength = gt_encseq_total_length(encseq);

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long countspecial = 0;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      unsigned long idx;

      for (idx = range.start; idx < range.end; idx++)
      {
        gt_assert(countspecial < specialcharacters && idx < totallength);
        suftabtail[countspecial++] = idx;
      }
    }
    gt_assert(countspecial == specialcharacters);
    gt_specialrangeiterator_delete(sri);
  }
  suftabtail[specialcharacters] = totallength;
}

static void gt_sain_rec_sortsuffixes(GtSaininfo *saininfo,
                                     unsigned long *suftab,
                                     unsigned long nonspecialentries,
                                     unsigned long suftabentries,
                                     unsigned int sainmode)
{
  unsigned long idx, *leftborder, numberofnames;

  leftborder = gt_malloc(sizeof (*leftborder) * saininfo->sainseq->numofchars);
  if (saininfo->countSstartype > 0)
  {
    gt_sain_endbuckets(leftborder,saininfo->sainseq->bucketsize,
                       saininfo->sainseq->numofchars);
    insertSstarsuffixes(saininfo, suftab, leftborder, nonspecialentries);
    gt_sain_startbuckets(leftborder,saininfo->sainseq->bucketsize,
                         saininfo->sainseq->numofchars);
    induceLtypesuffixes(saininfo, suftab, leftborder, nonspecialentries);
    gt_sain_endbuckets(leftborder,saininfo->sainseq->bucketsize,
                         saininfo->sainseq->numofchars);
    induceStypesuffixes(saininfo, suftab, leftborder, nonspecialentries);
    moveSstar2front(saininfo,suftab,nonspecialentries);
    numberofnames = assignSstarnames(saininfo,suftab,suftabentries);
    movenames2front(saininfo,suftab,suftabentries);
    gt_assert(numberofnames <= saininfo->countSstartype);
    if (numberofnames < saininfo->countSstartype)
    {
    /* Now the name sequence is in the range from
       saininfo->countSstartype .. 2 * saininfo->countSstartype - 1 */
      unsigned long position, *subseq = suftab + saininfo->countSstartype;
      GtSainseq *sainseq_rec;
      GtSaininfo *saininfo_rec;

      sainseq_rec = gt_sain_seq_new_from_array(subseq,
                                               saininfo->countSstartype,
                                               numberofnames);
      saininfo_rec = gt_sain_info_new(sainseq_rec);
      gt_sain_info_show(saininfo_rec);
      printf("recursively sort the named sequence of length %lu over %lu "
             "symbols (%.2f)\n",saininfo->countSstartype,numberofnames,
                         (double) numberofnames/saininfo->countSstartype);
      sain_setundefined(suftab,0,saininfo->countSstartype-1);
      gt_sain_rec_sortsuffixes(saininfo_rec,suftab,
                               saininfo->countSstartype,
                               saininfo->countSstartype,
                               sainmode);
      gt_sain_info_delete(saininfo_rec);
      gt_sain_seq_delete(sainseq_rec);
      for (idx = 0; idx < saininfo->countSstartype; idx++)
      {
        gt_assert(saininfo->countSstartype + suftab[idx] < suftabentries);
        suftab[saininfo->countSstartype + suftab[idx]] = idx;
      }
      for (idx = saininfo->countSstartype, position = 0;
           position < saininfo->sainseq->totallength;
           position++)
      {
        if (gt_sain_info_isSstartype(saininfo,position))
        {
          unsigned long putidx = suftab[idx];

          gt_assert(putidx < nonspecialentries);
          suftab[putidx] = position;
          idx++;
        }
      }
    }
  }
  if (sainmode == 1U && saininfo->countSstartype > 0)
  {
    gt_sain_checkorder(saininfo,suftab,0,saininfo->countSstartype-1);
  }
  sain_setundefined(suftab,saininfo->countSstartype,suftabentries - 1);
  gt_sain_endbuckets(leftborder,saininfo->sainseq->bucketsize,
                     saininfo->sainseq->numofchars);
  insertsortedSstarsuffixes(saininfo,suftab,leftborder,nonspecialentries);
  gt_sain_startbuckets(leftborder,saininfo->sainseq->bucketsize,
                       saininfo->sainseq->numofchars);
  induceLtypesuffixes(saininfo, suftab, leftborder, nonspecialentries);
  gt_sain_endbuckets(leftborder,saininfo->sainseq->bucketsize,
                     saininfo->sainseq->numofchars);
  induceStypesuffixes(saininfo, suftab, leftborder, nonspecialentries);
  gt_free(leftborder);
  if (nonspecialentries > 0)
  {
    if (sainmode == 1U)
    {
      gt_sain_checkorder(saininfo,suftab,0,nonspecialentries-1);
    } else
    {
      if (saininfo->sainseq->hasencseq && sainmode == 2U)
      {
        gt_sain_filltailsuffixes(suftab + nonspecialentries,
                                 saininfo->sainseq->seq.encseq);
        gt_suftab_lightweightcheck(saininfo->sainseq->seq.encseq,
                                   GT_READMODE_FORWARD,
                                   saininfo->sainseq->totallength,
                                   suftab,
                                   NULL);
      }
    }
  }
}

/* sainmode = 1: check with own gt_sain_checkorder
   sainmode = 2: check with gt_suftab_lightweightcheck
   sainmode = 3: no check
*/

void gt_sain_sortsuffixes(const GtEncseq *encseq,unsigned int sainmode)
{
  unsigned long nonspecialentries, requiredentries, suftabentries, *suftab;
  GtSainseq *sainseq = gt_sain_seq_new_from_encseq(encseq);
  GtSaininfo *saininfo;

  saininfo = gt_sain_info_new(sainseq);
  gt_sain_info_show(saininfo);
  nonspecialentries = sainseq->totallength - sainseq->specialcharacters;
  requiredentries = saininfo->countSstartype +
                    GT_DIV2(saininfo->sainseq->totallength) + 1;
  suftabentries = MAX(nonspecialentries,requiredentries);
  if (sainmode == 2U)
  {
    suftab = gt_malloc(sizeof (*suftab) * (saininfo->sainseq->totallength+1));
  } else
  {
    suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  }
  sain_setundefined(suftab,0,suftabentries - 1);
  gt_sain_rec_sortsuffixes(saininfo,suftab,nonspecialentries,suftabentries,
                           sainmode);
  gt_sain_info_delete(saininfo);
  gt_sain_seq_delete(sainseq);
  gt_free(suftab);
}
