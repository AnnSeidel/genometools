/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include "core/str_api.h"
#include "core/xansi_api.h"
#include "core/ma.h"
#include "core/log.h"
#include "core/fa.h"
#include "firstcodes-psbuf.h"
#include "firstcodes-spacelog.h"

GtLeftborderOutbuffer *gt_leftborderbuffer_new(GtFirstcodesspacelog *fcsl,
                                               bool useulong)
{
  GtLeftborderOutbuffer *lbbuf = gt_malloc(sizeof (*lbbuf));

  lbbuf->totalwrite = 0;
  lbbuf->outfilename = gt_str_new();
  lbbuf->fp = gt_xtmpfp(lbbuf->outfilename);
  lbbuf->nextfree = 0;
  lbbuf->allocated = 1024UL;
  lbbuf->useulong = useulong;
  if (useulong)
  {
    lbbuf->spaceulong = gt_malloc(sizeof (*lbbuf->spaceulong) *
                                  lbbuf->allocated);
    GT_FCI_ADDWORKSPACE(fcsl,"overflow_leftborderbuffer",
                        sizeof (*lbbuf->spaceulong) * lbbuf->allocated);
    lbbuf->spaceuint32_t = NULL;
  } else
  {
    lbbuf->spaceuint32_t = gt_malloc(sizeof (*lbbuf->spaceuint32_t) *
                                     lbbuf->allocated);
    GT_FCI_ADDWORKSPACE(fcsl,"leftborderbuffer",
                        sizeof (*lbbuf->spaceuint32_t) * lbbuf->allocated);
    lbbuf->spaceulong = NULL;
  }
  return lbbuf;
}

void gt_leftborderbuffer_flush(GtLeftborderOutbuffer *leftborderbuffer)
{
  if (leftborderbuffer->useulong)
  {
    gt_xfwrite(leftborderbuffer->spaceulong,
               sizeof (*leftborderbuffer->spaceulong),
               (size_t) leftborderbuffer->nextfree,
               leftborderbuffer->fp);
  } else
  {
    gt_xfwrite(leftborderbuffer->spaceuint32_t,
               sizeof (*leftborderbuffer->spaceuint32_t),
               (size_t) leftborderbuffer->nextfree,
               leftborderbuffer->fp);
  }
  leftborderbuffer->totalwrite += leftborderbuffer->nextfree;
  leftborderbuffer->nextfree = 0;
}

GtStr *gt_leftborderbuffer_delete(GtLeftborderOutbuffer *lbbuf,
                                  GtFirstcodesspacelog *fcsl,
                                  GT_UNUSED unsigned long expectedwritten)
{
  GtStr *outfilename;

  gt_assert(lbbuf != NULL);
  gt_leftborderbuffer_flush(lbbuf);
  gt_fa_fclose(lbbuf->fp);
  lbbuf->fp = NULL;
  if (lbbuf->useulong)
  {
    gt_log_log("write overflow_leftborder to file %s (%lu units of size %u)",
               gt_str_get(lbbuf->outfilename),
               lbbuf->totalwrite,(unsigned int) sizeof (*lbbuf->spaceulong));
    gt_assert(lbbuf->spaceulong != NULL);
    gt_free(lbbuf->spaceulong);
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"overflow_leftborderbuffer");
  } else
  {
    gt_log_log("write leftborder to file %s (%lu units of size %u)",
               gt_str_get(lbbuf->outfilename),
               lbbuf->totalwrite,(unsigned int) sizeof (*lbbuf->spaceuint32_t));
    gt_assert(lbbuf->spaceuint32_t != NULL);
    gt_free(lbbuf->spaceuint32_t);
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"leftborderbuffer");
  }
  gt_assert(lbbuf->totalwrite == expectedwritten);
  outfilename = lbbuf->outfilename;
  gt_free(lbbuf);
  return outfilename;
}