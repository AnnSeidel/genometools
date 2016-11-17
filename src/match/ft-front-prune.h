#ifndef FT_FRONT_PRUNE_H
#define FT_FRONT_PRUNE_H
#ifndef OUTSIDE_OF_GT
#include "match/ft-trimstat.h"
#include "match/ft-polish.h"
#include "match/ft-front-generation.h"
#include "core/encseq_api.h"
#else
#include "gt-defs.h"
#include "polish.h"
#include "front-generation.h"
#include "trimstat.h"
#endif

typedef struct
{
  void *space;
  GtUword offset, allocated;
} GtAllocatedMemory;

typedef enum
{
  GT_OUTSENSE_TRIM_ALWAYS,
  GT_OUTSENSE_TRIM_ON_NEW_PP,
  GT_OUTSENSE_TRIM_NEVER
} GtTrimmingStrategy;

#ifndef OUTSIDE_OF_GT
typedef enum
{
  GT_EXTEND_CHAR_ACCESS_ENCSEQ,
  GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER,
  GT_EXTEND_CHAR_ACCESS_DIRECT,
  GT_EXTEND_CHAR_ACCESS_ANY
} GtExtendCharAccess;

typedef struct
{
  const GtEncseq *encseq;
  GtReadmode readmode;
  GtAllocatedMemory *sequence_cache;
  GtEncseqReader *encseq_r;
  const GtUchar *bytesequence;
  GtExtendCharAccess extend_char_access;
  GtUword totallength;
} GtFTsequenceResources;
#endif

GtUword front_prune_edist_inplace(
#ifndef OUTSIDE_OF_GT
                       bool forward,
                       GtAllocatedMemory *frontspace_reservoir,
#endif
                       GtFtTrimstat *trimstat,
                       GtFtPolished_point *best_polished_point,
                       GtFronttrace *fronttrace,
                       const GtFtPolishing_info *pol_info,
                       GtTrimmingStrategy trimstrategy,
                       GtUword history,
                       GtUword minmatchnum,
                       GtUword maxalignedlendifference,
                       bool showfrontinfo,
                       GtUword seedlength,
                       GtFTsequenceResources *ufsr,
                       GtUword ustart,
                       GtUword uulen,
                       GtUword vseqstartpos,
                       GtFTsequenceResources *vfsr,
                       GtUword vstart,
                       GtUword vlen);

#endif
