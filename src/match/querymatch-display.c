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
#include <stdbool.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "match/querymatch-display.h"

typedef enum
{
  Gt_Alignment_display,
  Gt_Cigarstring_display,
  Gt_Polishinginfo_display,
  Gt_Fstperquery_display,
  Gt_Seed_display,
  Gt_Failed_Seed_display,
  Gt_Seed_in_alignment_display,
  Gt_Seqlength_display,
  Gt_Evalue_display,
  Gt_Seqdesc_display,
  Gt_Bitscore_display
} GtSeedExtendDisplay_enum;

struct GtSeedExtendDisplayFlag
{
  unsigned int flags;
  bool a_seedpos_relative, b_seedpos_relative;
  GtUword alignmentwidth;
};

GtSeedExtendDisplayFlag *gt_querymatch_display_flag_new(void)
{
  GtSeedExtendDisplayFlag *display_flag = gt_malloc(sizeof *display_flag);

  display_flag->flags = 0;
  display_flag->alignmentwidth = 0;
  display_flag->a_seedpos_relative = true; /* as bytes is default access mode */
  display_flag->b_seedpos_relative = true;
  return display_flag;
}

static unsigned int gt_display_mask(int shift)
{
  return 1U << shift;
}

static bool gt_querymatch_display_on(const GtSeedExtendDisplayFlag
                                       *display_flag,
                                     GtSeedExtendDisplay_enum display)
{
  gt_assert((int) display <= Gt_Bitscore_display);
  return (display_flag != NULL &&
          (display_flag->flags & gt_display_mask(display))) ? true : false;
}

bool gt_querymatch_seed_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seed_display);
}

bool gt_querymatch_failed_seed_display(const GtSeedExtendDisplayFlag
                                         *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Failed_Seed_display);
}

bool gt_querymatch_cigarstring_display(const GtSeedExtendDisplayFlag
                                        *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Cigarstring_display);
}

bool gt_querymatch_seed_in_alignment_display(
                 const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seed_in_alignment_display);
}

bool gt_querymatch_evalue_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Evalue_display);
}

bool gt_querymatch_bitscore_display(const GtSeedExtendDisplayFlag
                                       *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Bitscore_display);
}

bool gt_querymatch_seqdesc_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seqdesc_display);
}

bool gt_querymatch_seqlength_display(const GtSeedExtendDisplayFlag
                                       *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Seqlength_display);
}

bool gt_querymatch_polinfo_display(const GtSeedExtendDisplayFlag *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Polishinginfo_display);
}

bool gt_querymatch_fstperquery_display(const GtSeedExtendDisplayFlag
                                            *display_flag)
{
  return gt_querymatch_display_on(display_flag,Gt_Fstperquery_display);
}

GtUword gt_querymatch_display_alignmentwidth(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag == NULL) ? 0 : display_flag->alignmentwidth;
}

bool gt_querymatch_display_alignment(const GtSeedExtendDisplayFlag
                                     *display_flag)
{
  return (display_flag != NULL &&
          display_flag->alignmentwidth > 0) ? true : false;
}

void gt_querymatch_display_seedpos_relative_set(GtSeedExtendDisplayFlag
                                                *display_flag,
                                                bool a_is_rel,
                                                bool b_is_rel)
{
  gt_assert(display_flag != NULL);
  display_flag->a_seedpos_relative = a_is_rel;
  display_flag->b_seedpos_relative = b_is_rel;
}

bool gt_querymatch_display_seedpos_a_relative(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag != NULL && display_flag->a_seedpos_relative) ? true
                                                                    : false;
}

bool gt_querymatch_display_seedpos_b_relative(const GtSeedExtendDisplayFlag
                                                *display_flag)
{
  return (display_flag != NULL && display_flag->b_seedpos_relative) ? true
                                                                    : false;
}

GtStr *gt_querymatch_column_header(const GtSeedExtendDisplayFlag *display_flag)
{
  GtStr *str = gt_str_new();

  if (gt_querymatch_seed_display(display_flag))
  {
    gt_str_append_cstr(str,", seedlen, s.seedstart, q.seedstart");
  }
  if (gt_querymatch_seqlength_display(display_flag))
  {
    gt_str_append_cstr(str,", s.seqlen, q.seqlen");
  }
  if (gt_querymatch_evalue_display(display_flag))
  {
    gt_str_append_cstr(str,", evalue");
  }
  if (gt_querymatch_bitscore_display(display_flag))
  {
    gt_str_append_cstr(str,", bit score");
  }
  if (gt_querymatch_cigarstring_display(display_flag))
  {
    gt_str_append_cstr(str,", cigarstring");
  }
  return str;
}

void gt_querymatch_fields_approx_output(const GtSeedExtendDisplayFlag
                                         *display_flag,FILE *stream)
{
  GtStr *add_column_header;

  fprintf(stream,"# Fields: s.len, s.seqnum, s.start, strand, q.len, q.seqnum, "
                 "q.start, score, editdist, identity");
  add_column_header = gt_querymatch_column_header(display_flag);
  if (gt_str_length(add_column_header) > 0)
  {
    fputs(gt_str_get(add_column_header),stream);
  }
  fputc('\n',stream);
  gt_str_delete(add_column_header);
}

void gt_querymatch_fields_exact_output(FILE *stream)
{
  fprintf(stream,"# Fields: s.len, s.seqnum, s.start, strand, q.seqnum, "
                 "q.start\n");
}

static GtStrArray *gt_se_help2display_strings(const char *helpline)
{
  const char *ptr = helpline, *laststart = NULL;
  char findchar = '\n';
  GtStrArray *display_args = gt_str_array_new();

  while (true)
  {
    ptr = strchr(ptr,findchar);
    gt_assert(ptr != NULL);
    if (findchar == '\n')
    {
      if (*(ptr+1) == '\0')
      {
        break;
      }
      ptr++;
      if (*ptr != ' ')
      {
        laststart = ptr;
        findchar = ':';
      }
    } else
    {
      GtUword len = (GtUword) (ptr - laststart);
      gt_str_array_add_cstr_nt(display_args,laststart,len);
      findchar = '\n';
    }
  }
  return display_args;
}

/* The following help line defines which keywords can be used as arguments
   to option -outfmt. Each keyword follows a newline and ends with a :. This
   is exploited by generating the list of possible keyword for
   checking the -outfmt argument ist generated by
   gt_se_help2display_strings */

const char *gt_querymatch_display_help(void)
{
  return "specify what information about the matches to display\n"
         "alignment:    display alignment (possibly followed by =<number>\n"
         "              to specify width of alignment columns)\n"
         "cigar:        show cigar string representing alignment\n"
         "polinfo:      display polishing information for displayed\n"
         "              alignment\n"
         "fstperquery:  output only the first found match per query\n"
         "seed:         display the seed of the match, i.e. the length and\n"
         "              the start position of the seed in both instances\n"
         "failed_seed:  display the seed of the match that was extended,\n"
         "              but failed (after extension) the filter conditions\n"
         "seed_in_algn: display the seed in alignment\n"
         "seqlength:    display length of sequences in which\n"
         "              the two match-instances occur\n"
         "evalue:       display evalue\n"
         "seqdesc:      display sequence description instead of numbers\n"
         "bitscore:     display bit score\n";
}

static int gt_querymatch_display_flag_set(GtWord *parameter,
                                          const GtStrArray *display_strings,
                                          GtSeedExtendDisplayFlag *display_flag,
                                          const char *arg,
                                          GtError *err)
{
  GtUword ds_idx, numofds = gt_str_array_size(display_strings);
  const GtSeedExtendDisplay_enum exclude_list[] = {Gt_Alignment_display,
                                                   Gt_Cigarstring_display};
  size_t ex_idx, numexcl = sizeof exclude_list/sizeof exclude_list[0];
  bool identifier_okay = false, parameter_found = false;
  const char *ptr;
  size_t cmplen;

  gt_assert(display_flag != NULL &&
            numofds == (size_t) Gt_Bitscore_display + 1 &&
            numexcl % 2 == 0);
  ptr = strchr(arg,'=');
  if (ptr != NULL)
  {
    cmplen = (size_t) (ptr - arg);
    if (sscanf(ptr+1,GT_WD,parameter) != 1)
    {
      gt_error_set(err,"illegal argument \"%s\" to option -outfmt: "
                       "expect integer following symbol =",arg);
      return -1;
    }
    parameter_found = true;
  } else
  {
    cmplen = 0;
  }
  for (ds_idx = 0; ds_idx < numofds; ds_idx++)
  {
    const char *dstring = gt_str_array_get(display_strings,ds_idx);
    int ret = (cmplen > 0) ? strncmp(arg,dstring,cmplen)
                           : strcmp(arg,dstring);
    if (ret == 0)
    {
      display_flag->flags |= gt_display_mask((int) ds_idx);
      identifier_okay = true;
      break;
    }
  }
  if (!identifier_okay)
  {
    GtStr *err_msg = gt_str_new();

    gt_str_append_cstr(err_msg,
                       "illegal identifier as argument of option -outfmt: "
                       "possible idenfifiers are: ");
    for (ds_idx = 0; ds_idx < numofds; ds_idx++)
    {
      const char *dstring = gt_str_array_get(display_strings,ds_idx);
      gt_str_append_cstr(err_msg,dstring);
      if (ds_idx < numofds - 1)
      {
        gt_str_append_cstr(err_msg,", ");
      }
    }
    gt_error_set(err,"%s",gt_str_get(err_msg));
    gt_str_delete(err_msg);
    return -1;
  }
  for (ex_idx = 0; ex_idx < numexcl; ex_idx+=2)
  {
    if ((display_flag->flags & gt_display_mask(exclude_list[ex_idx])) &&
        (display_flag->flags & gt_display_mask(exclude_list[ex_idx+1])))
    {
      const char *d1 = gt_str_array_get(display_strings,exclude_list[ex_idx]);
      const char *d2 = gt_str_array_get(display_strings,exclude_list[ex_idx+1]);
      gt_error_set(err,"argument \"%s\" and \"%s\" of option -outfmt exclude "
                       "each other",d1,d2);
      return -1;
    }
  }
  return parameter_found ? 1 : 0;
}

int gt_querymatch_display_flag_args_set(GtSeedExtendDisplayFlag *display_flag,
                                        const GtStrArray *display_args,
                                        GtError *err)
{
  bool haserr = false;
  GtUword da_idx;
  GtStrArray *display_strings
    = gt_se_help2display_strings(gt_querymatch_display_help());

  gt_assert(display_flag != NULL);
  for (da_idx = 0; da_idx < gt_str_array_size(display_args); da_idx++)
  {
    const char *da = gt_str_array_get(display_args,da_idx);
    GtWord parameter;
    int ret = gt_querymatch_display_flag_set(&parameter,display_strings,
                                             display_flag,da,err);
    switch (ret)
    {
      case 0:
        break;
      case 1:
        /* the only flag with a parameter is Gt_Alignment_display */
        if (parameter < 0)
        {
          gt_error_set(err,"integer following \"alignment=\" must be positive");
          haserr = true;
        } else
        {
          gt_assert(display_flag->flags &
                    gt_display_mask(Gt_Alignment_display));
          display_flag->alignmentwidth = (GtUword) parameter;
        }
        break;
      default:
        haserr = true;
    }
  }
  if (!haserr &&
      (display_flag->flags & gt_display_mask(Gt_Alignment_display)) &&
      display_flag->alignmentwidth == 0)
  {
    display_flag->alignmentwidth = 60; /* this is the default alignment width */
  }
  gt_str_array_delete(display_strings);
  return haserr ? -1 : 0;
}

void gt_querymatch_display_flag_delete(GtSeedExtendDisplayFlag *display_flag)
{
  if (display_flag != NULL)
  {
    gt_free(display_flag);
  }
}
