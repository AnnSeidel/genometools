/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef UNDEF_H
#define UNDEF_H

#include <float.h>
#include <limits.h>

#define UNDEF_BOOL         (bool) ~0
#define UNDEF_CHAR         CHAR_MAX
#define UNDEF_DOUBLE       DBL_MAX
#define UNDEF_FLOAT        FLT_MAX
#define UNDEF_INT          ~0
#define UNDEF_UCHAR        UCHAR_MAX
#define UNDEF_UINT         ~0U
#define UNDEF_LONG         LONG_MIN
#define UNDEF_ULONG        ~0UL

#define UNDEF_SCORE        UNDEF_FLOAT

#endif
