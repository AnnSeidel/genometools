/*
  Copyright (c) 2016 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2016 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/alphabet_api.h"
#include "core/ensure.h"
#include "core/error.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "match/karlin_altschul_stat.h"
#include "extended/scorehandler.h"

#define K_ITER_MAX 100
#define K_SUMLIMIT_DEFAULT 0.0001

/*
  this library implements calculation of karlin-altschul parameter for E-value
  of Alignments analog to NCBI tool BLAST:

  Altschul S.F., Gish W., Miller W., Myers E.W. and Lipman D.J. (1990)
  Basic local alignment search tool. J. Mol. Biol. 215: 403-410.
 */

/* stores karlin altschul parameters */
struct GtKarlinAltschulStat
{
  double lambda,
         K,
         logK,
         H,
         alpha_div_lambda,
         beta;
  GtWord matchscore, mismatchscore, gapscore;
};

typedef struct{
  double *sprob,
         score_avg;
  GtWord low_align_score,
         high_align_score;
} ScoringFrequency;

/* stores letter propabilities */
typedef struct{
  char   ch;
  double p;
} LetterProb;

/* provisional solution only dna alphabet */
static LetterProb nt_prob[] = {
  { 'A', 0.25 },
  { 'C', 0.25 },
  { 'G', 0.25 },
  { 'T', 0.25 }
};

/*
  precomputed values
  analog to BLAST

  1. Gap opening score,
  2. Gap extension score,
  3. Lambda,
  4. K,
  5. H,
  6. Alpha,
  7. Beta,
 */

typedef enum
{
  gapopidx,
  gapextdidx,
  lambdaidx,
  Kidx,
  Hidx,
  alphaidx,
  betaidx
} GT_ValuesindeX;

typedef double GA_Values[(int) (betaidx+1)];

/* matchscore = 1 && mismatchscore = -4 */
static const GA_Values ga_matrix_1_4[] = {
    { 0, -2,  1.26,  0.43, 0.90,  1.4, -1}
};

/* matchscore = 2 && mismatchscore = -7 */
static const GA_Values ga_matrix_2_7[] = {
    { 0, -4,  0.63, 0.43, 0.90,   0.7, -1}
};

/* matchscore = 1 && mismatchscore = -3 */
static const GA_Values ga_matrix_1_3[] = {
    { 0, -2,  1.25,  0.42, 0.83,  1.5, -2}
};

/* matchscore = 2 && mismatchscore = -5 */
static const GA_Values ga_matrix_2_5[] = {
    { 0, -4,  0.62, 0.39, 0.78,  0.8, -2}
};

/* matchscore = 1 && mismatchscore = -2 */
static const GA_Values ga_matrix_1_2[] = {
    { 0, -2, 1.19, 0.34, 0.66, 1.8, -3}
};

/* matchscore = 2 && mismatchscore = -3 */
static const GA_Values ga_matrix_2_3[] = {
    { 0, -4,  0.55, 0.21, 0.46,  1.2, -5}
};

void gt_karlin_altschul_stat_delete(GtKarlinAltschulStat *ka)
{
  gt_free(ka);
}

double gt_karlin_altschul_stat_get_lambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->lambda;
}

double gt_karlin_altschul_stat_get_logK(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->logK;
}

double gt_karlin_altschul_stat_get_K(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->K;
}

double gt_karlin_altschul_stat_get_alphadlambda(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->alpha_div_lambda;
}

double gt_karlin_altschul_stat_get_beta(const GtKarlinAltschulStat *ka)
{
  gt_assert(ka);
  return ka->beta;
}

/* calculate probabilities of scores */
static ScoringFrequency *gt_karlin_altschul_stat_scoring_frequency(
                                             const GtAlphabet *alphabet,
                                             const GtScoreHandler *scorehandler)
{
  unsigned int idx, jdx, numofchars;
  GtWord score, obs_min = 0, obs_max = 0, range;
  double score_avg, score_sum;

  gt_assert(alphabet && scorehandler);

  /* make generalizations of alphabet probabilities, for now nt_prob */
  gt_assert(gt_alphabet_is_dna(alphabet));

  ScoringFrequency *sf = gt_malloc(sizeof (*sf));
  gt_assert(sf);

  numofchars = gt_alphabet_num_of_chars(alphabet);
  for (idx = 0; idx < numofchars; idx++)
  {
    for (jdx = 0; jdx < numofchars; jdx++)
    {
      score = gt_scorehandler_get_replacement(scorehandler, idx, jdx);
      obs_min = MIN(obs_min, score);
      obs_max = MAX(obs_max, score);
    }
  }

  /* for theoretically valid scoring systems */
  gt_assert(obs_min <= 0 && obs_max >= 0);
  sf->low_align_score = obs_min;
  sf->high_align_score = obs_max;

  range = obs_max - obs_min + 1;
  sf->sprob = gt_calloc(range, sizeof (*sf->sprob));

  for (idx = 0; idx < numofchars; idx++)
  {
    for (jdx = 0; jdx < numofchars; jdx++)
    {
      score = gt_scorehandler_get_replacement(scorehandler, idx, jdx);

      if (score >= sf->low_align_score)
        sf->sprob[score-sf->low_align_score] += nt_prob[idx].p * nt_prob[jdx].p;
    }
  }

  score_sum = 0.0;
  for (score = obs_min; score <= obs_max; score++)
  {
    if (sf->sprob[score-obs_min] > 0)
      score_sum += sf->sprob[score-obs_min];
  }

  score_avg = 0.0;
  for (score = obs_min; score <= obs_max; score++)
  {
    sf->sprob[score-obs_min] /= score_sum;
    score_avg += score * sf->sprob[score-obs_min];
  }
  sf->score_avg = score_avg;

  return sf;
}

static double gt_karlin_altschul_stat_calculate_ungapped_lambda(
                                                     const ScoringFrequency *sf)
{
  double x0, x, lambda, tolerance, q, dq;
  GtWord low, high;
  GtUword k_max_iterations, idx, jdx;

  /* solve phi(lambda) = -1 + sum_{i=l}^{u} sprob(i)*exp(i*lambda) = 0 */

  x0 = 0.5; /* x0 in (0,1) */
  tolerance = 1.e-5;
  k_max_iterations = 20;

  low = sf->low_align_score;
  high = sf->high_align_score;

  /* write phi as phi(lambda) = exp(u*lambda) * q(exp(-lambda)) and solve the
    polynomial q by apply newton's method
    q(x) = -x^u + sum_{k=0}^{u-l} sprob(u-k)* x^k */
  for (idx = 0; idx < k_max_iterations; idx++)
  {
    q = -pow(x0,high);
    for (jdx = 0; jdx <= high-low; jdx++)
      q += sf->sprob[high-low-jdx] * pow(x0,jdx);

    dq = -high*pow(x0,(high-1));
    for (jdx = 1; jdx <= high-low; jdx++)
      dq += sf->sprob[high-low-jdx] * jdx * pow(x0,jdx-1);

    x = x0 - (q/dq);

    if (fabs(x-x0) < tolerance)
      break;

    x0 = x;
  }

  lambda = -log(x);

  /* better solution would be to apply Horner's rule for evaluating a
     polynomial and its derivative (s. BLAST), but for the moment it works */

  return lambda;
}

static double gt_karlin_altschul_stat_calculate_H(const ScoringFrequency *sf,
                                                  double lambda)
{
  double H, sum, etonlambda;
  GtWord idx, low, high, scale;
  gt_assert(sf->sprob);

  low = sf->low_align_score;
  high = sf->high_align_score;

  etonlambda = exp(-lambda);
  sum = low * sf->sprob[0];
  for (idx = low + 1; idx <= high; idx++)
  {
    sum = idx * sf->sprob[idx-low] + etonlambda * sum;
  }

  scale = pow(etonlambda,high);
  if (scale > 0.0)
    H = lambda * sum/scale;
  else /* case underflow */
    H = lambda * exp(lambda * high + log(sum));

  return H;
}

static GtWord gt_karlin_altschul_stat_gcd(const ScoringFrequency *sf)
{
  GtWord idx, range, div;

  range = sf->high_align_score - sf->low_align_score+1;
  div = -sf->low_align_score;
  for (idx = 1; idx < range && div > 1; idx++)
  {
    if (sf->sprob[idx] != 0.0)
    {
      GtWord val = labs(idx+sf->low_align_score);
      if (val > div)
      {
        GtWord tmp = div;
        div = val;
        val = tmp;
      }
      while (val != 0)
      {
        GtWord tmp = div % val;
        div = val;
        val = tmp;
      }
    }
  }
  return div;
}

static double gt_karlin_altschul_stat_calculate_ungapped_K(
                                                     const ScoringFrequency *sf,
                                                     double lambda,
                                                     double H)
{
  GtWord  low,
          high,
          div,
          low_align_score,
          high_align_score,
          count, idx, jdx, firstidx, lastidx, secondidx,
          first, last;
  GtUword range,
          iterlimit,
          size,
          sigma = 0;
  double  score_avg,
          score_avg_div,
          one_minus_expnlambda,
          *alignnment_score_probs,
          expnlambda,
          K,
          sumlimit,
          inner_sum;

  gt_assert(lambda > 0 && H > 0);

  score_avg = sf->score_avg;
  gt_assert(score_avg < 0.0);

  /* greatest common divisor */
  div = gt_karlin_altschul_stat_gcd(sf);

  low = sf->low_align_score/div;
  high = sf->high_align_score/div;
  lambda *= div;

  range = high - low;
  expnlambda = exp(-lambda);

  if (low == -1 && high == 1)
  {
    K = (sf->sprob[0] - sf->sprob[sf->high_align_score-sf->low_align_score]) *
        (sf->sprob[0] - sf->sprob[sf->high_align_score-sf->low_align_score])/
        sf->sprob[0];
  }
  else if (low == -1 || high == 1)
  {
    one_minus_expnlambda = 1-exp(-lambda);
    if (high != 1)
    {
      score_avg_div = score_avg / div;
      K = lambda * one_minus_expnlambda / H *(score_avg_div * score_avg_div);
    }
    else
    {
      K = H/lambda * one_minus_expnlambda;
    }
  }
  else
  {
    /* K = lambda*exp(-2*sigma)/(H*(1-exp(-lambda)) */

    sumlimit = K_SUMLIMIT_DEFAULT;
    iterlimit = K_ITER_MAX;

    size = iterlimit * range + 1;
    alignnment_score_probs = gt_calloc(size, sizeof (*alignnment_score_probs));
    gt_assert(alignnment_score_probs);

    low_align_score = 0;
    high_align_score = 0;
    inner_sum = 1.0;
    alignnment_score_probs[0] = 1.0;

    for (count = 0; count < iterlimit && inner_sum > sumlimit; count++)
    {
      if (count > 0)
      {
        inner_sum /= count;
        sigma += inner_sum;
      }

      first = last = range;
      low_align_score += low;
      high_align_score += high;
      for (idx = high_align_score-low_align_score; idx >= 0; idx--)
      {
        firstidx = idx-first;
        lastidx = idx-last;
        secondidx = sf->sprob[low] + first;
        for (inner_sum = 0.; firstidx >= lastidx; )
        {
          inner_sum += alignnment_score_probs[firstidx] *
                      alignnment_score_probs[secondidx];
          firstidx--;
          secondidx++;
        }

        if (first > 0)
          --first;
        if (idx <= range)
          --last;

        alignnment_score_probs[idx] = inner_sum;
        if (idx == 0)
          break;
      }

      inner_sum = alignnment_score_probs[++idx];
      for (jdx = low_align_score + 1; jdx < 0; jdx++)
      {
        inner_sum = alignnment_score_probs[++idx] + inner_sum * expnlambda;
      }
      inner_sum *= expnlambda;
      for (/*Nothing*/; jdx <= high_align_score; ++jdx)
        inner_sum += alignnment_score_probs[++jdx];

      gt_free(alignnment_score_probs);

    }
    /* no terms of geometric progression, check to add these terms for
       correction in future */

    K = -exp(-2.0*sigma)/(H/lambda*expnlambda);
  }

  return K;
}

static int get_values_from_matrix(GtKarlinAltschulStat *ka,
                                  GA_Values *matrix,
                                  GtUword length,
                                  GtWord gap_extension,
                                  GT_UNUSED GtWord gap_open,
                                  GtError *err)
{
  GtUword idx;

  for (idx = 0; idx < length; idx++)
  {
    if (matrix[idx][gapextdidx] == gap_extension)
    {
      if (matrix[idx][gapopidx] != 0)
      {
        gt_assert(false); /* not implemented: linear scores only */
      }
      ka->lambda = matrix[idx][lambdaidx];
      ka->K = matrix[idx][Kidx];
      ka->logK = log(ka->K);
      ka->H = matrix[idx][Hidx];

      gt_assert(ka->lambda != 0.0);
      ka->alpha_div_lambda = matrix[idx][alphaidx]/ka->lambda;
      ka->beta = matrix[idx][betaidx];
      return 0;
    }
  }

  /* not found */
  gt_error_set(err, "no precomputed values for ungapped alignment parameters "
                    "were found\n");
  return 1;
}

static int gt_karlin_altschul_stat_get_gapped_params(GtKarlinAltschulStat *ka,
                                                     const GtScoreHandler
                                                           *scorehandler,
                                                     GtError *err)
{
  GtWord gap_open, gap_extension, matchscore, mismatchscore;
  GA_Values *ga_matrix = NULL;
  GtUword length = 0;

  gt_assert(ka && scorehandler);
  matchscore = gt_scorehandler_get_matchscore(scorehandler);
  mismatchscore = gt_scorehandler_get_mismatchscore(scorehandler);

  if (matchscore == 1 && mismatchscore == -4)
  {
    length = sizeof (ga_matrix_1_4) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_1_4;
  }
  else if (matchscore == 2 && mismatchscore == -7)
  {
    length = sizeof (ga_matrix_2_7) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_2_7;
  }
  else if (matchscore == 1 && mismatchscore == -3)
  {
    length = sizeof (ga_matrix_1_3) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_1_3;
  }
  else if (matchscore == 2 && mismatchscore == -5)
  {
    length = sizeof (ga_matrix_2_5) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_2_5;
  }
  else if (matchscore == 1 && mismatchscore == -2)
  {
    length = sizeof (ga_matrix_1_2) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_1_2;
  }
  else if (matchscore == 2 && mismatchscore == -3)
  {
    length = sizeof (ga_matrix_2_3) / sizeof (GA_Values);
    ga_matrix = (GA_Values*) ga_matrix_2_3;
  }
  else
  {
    gt_error_set(err, "no precomputed values for combination matchscore "
                      GT_WD" and mismatchscore "GT_WD" in evalue calculation "
                      "of gapped alignments\n",
                      matchscore, mismatchscore);
    return 1;
  }

  gap_extension = gt_scorehandler_get_gapscore(scorehandler);
  gap_open = gt_scorehandler_get_gap_opening(scorehandler);
  if (get_values_from_matrix(ka,
                         ga_matrix,
                         length,
                         gap_extension,
                         gap_open,
                         err))
  {
    return 1;
  }

  return 0;
}

GtKarlinAltschulStat *gt_karlin_altschul_stat_new(bool gapped_alignment,
                                                  const GtAlphabet *alphabet,
                                                  const GtScoreHandler
                                                        *scorehandler,
                                                  GtError *err)
{
  GtKarlinAltschulStat *ka;
  ka = gt_malloc(sizeof (GtKarlinAltschulStat));
  gt_assert(ka);

  ka->lambda = 0;
  ka->K = 0;
  ka->logK = 0;
  ka->H = 0;

  gt_assert(gt_scorehandler_get_gap_opening(scorehandler) == 0);
  /* only implemented for linear */
  if (!gapped_alignment)
  {
    /* New ScoringFrequency */
    ScoringFrequency *sf
      = gt_karlin_altschul_stat_scoring_frequency(alphabet,scorehandler);
    gt_assert(sf != NULL);

    /* karlin altschul parameters for ungapped alignments */
    ka->lambda = gt_karlin_altschul_stat_calculate_ungapped_lambda(sf);
    ka->H = gt_karlin_altschul_stat_calculate_H(sf, ka->lambda);
    ka->K = gt_karlin_altschul_stat_calculate_ungapped_K(sf, ka->lambda, ka->H);
    ka->logK = log(ka->K);
    gt_assert(ka->H != 0.0);
    ka->alpha_div_lambda = 1/ka->H;
    ka->beta = 0;

    gt_free(sf->sprob);
    gt_free(sf);

  }
  else /* gapped alignments */
  {
    if (gt_karlin_altschul_stat_get_gapped_params(ka, scorehandler, err))
    {
      gt_karlin_altschul_stat_delete(ka);
      return NULL;
    }
  }
  ka->matchscore = gt_scorehandler_get_matchscore(scorehandler);
  ka->mismatchscore = gt_scorehandler_get_mismatchscore(scorehandler);
  ka->gapscore = gt_scorehandler_get_gapscore(scorehandler);
  return ka;
}

GtWord gt_karlin_altschul_stat_mismatchscore(const GtKarlinAltschulStat *ka)
{
  return ka->mismatchscore;
}

GtWord gt_karlin_altschul_stat_matchscore(const GtKarlinAltschulStat *ka)
{
  return ka->matchscore;
}

GtWord gt_karlin_altschul_stat_gapscore(const GtKarlinAltschulStat *ka)
{
  return ka->gapscore;
}

int gt_karlin_altschul_stat_unit_test(GtError *err)
{
  GtKarlinAltschulStat *ka;
  GtScoreHandler *scorehandler;
  GtAlphabet *alphabet;
  double q;

  int had_err = 0;
  gt_error_check(err);

  scorehandler = gt_scorehandler_new(1,-2,0,-2);

  /* check function for gapped alignments */
  ka = gt_karlin_altschul_stat_new(true, NULL, scorehandler, err);
  gt_ensure(ka->lambda == 1.19);
  gt_ensure(ka->H == 0.66);
  gt_ensure(ka->K == 0.34);
  gt_karlin_altschul_stat_delete(ka);

  /* check function for ungapped alignments */
  alphabet = gt_alphabet_new_dna();
  ka = gt_karlin_altschul_stat_new(false, alphabet, scorehandler, err);
  q = ka->lambda/1.33; /* lambda = 1.33 */
  gt_ensure(0.99 < q && q < 1.01);
  q = ka->H/1.12; /* H = 1.12 */
  gt_ensure(0.99 < q && q < 1.01);
  q = ka->K/0.621; /* K = 0.621 */
  gt_ensure(0.99 < q && q < 1.01);
  gt_karlin_altschul_stat_delete(ka);

  gt_scorehandler_delete(scorehandler);
  gt_alphabet_delete(alphabet);
  return had_err;
}
