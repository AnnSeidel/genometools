/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ENUM_PATT_DEF_H
#define ENUM_PATT_DEF_H

typedef struct Enumpatterniterator Enumpatterniterator;

Enumpatterniterator *newenumpatterniterator(unsigned long minpatternlen,
                                            unsigned long maxpatternlen,
                                            const Encodedsequence *encseq,
                                            GtError *err);

const GtUchar *nextEnumpatterniterator(unsigned long *patternlen,
                                     Enumpatterniterator *epi);

void showPatterndistribution(const Enumpatterniterator *epi);

void freeEnumpatterniterator(Enumpatterniterator **epi);

#endif
