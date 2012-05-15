/*Library providing C support to poppy

Copyright 2010 Los Alamos National Security, LLC.
*/

#include "randomkit.h"

/*Performs inner loop of association analysis
 *p2: list of times in process 2 (count: n_p2)
 *p1: list of times in process 1 (count: n_p1)
 *lags: lag times between p1 and p2 (n_lags)
 *n_assoc: 2D array of resulting association numbers (n_lags, n_p1)
 *winhalf: window half-size between p1, p2
 */
void assoc(double* p2, double* p1, double* lags, long* n_assoc,
	   double winhalf, long n_p2, long n_p1, long n_lags)
{
  long j, curr_start;
  double laggedp1;
  double *currp1, *maxp1;
  double* lagend = lags + n_lags;
  maxp1 = p1 + n_p1;
  for(; lags < lagend; lags++) {
    curr_start = 0;
    for(currp1 = p1; currp1 < maxp1; currp1++, n_assoc++) {
      laggedp1 = *currp1 + *lags;
      for (;
	   curr_start < n_p2 && p2[curr_start] < laggedp1 - winhalf;
	   curr_start++);
      for(j = curr_start;
	  j < n_p2 && p2[j] <= laggedp1 + winhalf;
	  j++);
      *n_assoc = j - curr_start;
    }
  }
}

/*Performs inner loop of bootstrapping
 *Bootstraps n_els elements of data, n times
 *surrogates: list of surrogates (n * n_els)
 *data: list of original data (n_els)
 *n: number of bootstraps
 *n_els: number of elements
 *seed: number to seed RNG
 *clock_seed: if true, seed RNG from the clock, ignore seed
 */
void boots(double* surrogates, double* data,
	   unsigned long n, unsigned long n_els,
	   unsigned long sample_size,
	   unsigned long seed, int clock_seed)
{
  rk_state state;
  double *curr_surr, *last_surr;
  unsigned long max = n_els - 1;
  if(clock_seed)
    rk_randomseed(&state);
  else
    rk_seed(seed, &state);
  curr_surr = surrogates;
  last_surr = surrogates + (n * sample_size);
  for(curr_surr=surrogates; curr_surr < last_surr; curr_surr++)
    *curr_surr = data[rk_interval(max, &state)];
}

/*Performs inner loop of association analysis confidence interval, including
 *bootstrapping and summing on the bootstraps.
 *n_assoc: number of associations for each lag, each time in process1 
 *         (n_lags * n_p1)
 *surr_assoc_total: (out) total number of associations for each lag, each
 *                  surrogate series (n_lags * n_surr)
 *n_lags: number of time lags
 *n_p1: number of times in process 1, i.e. number of association numbers for
 *      each lag
 *n_surr: number of surrogate series, i.e. number of times to execute the
 *        bootstrap for each lag
 *seeds: numbers to seed RNG, reseeded each lag (n_lags)
 *clock_seed: if true, seed RNG from the clock, ignore seeds
 */
void aa_ci(unsigned long *n_assoc, unsigned long *surr_assoc_total,
	   unsigned long n_lags, unsigned long n_p1, unsigned long n_surr,
	   unsigned long *seeds, int clock_seed)
{
  unsigned long sample; /*Sample number*/
  unsigned long *curr_assoc; /*Association numbers for the current lag*/
  unsigned long *max_assoc; /*Just past last value of association numbers*/
  unsigned long *totalp; /*Point to total assoc. no.; current lag, surrogate*/
  unsigned long *max_total; /*Just past last value of total for this lag*/
  rk_state state;
  unsigned long max = n_p1 - 1;
  if(clock_seed)
    rk_randomseed(&state);
  
  max_assoc = n_assoc + n_p1 * n_lags;
  totalp = surr_assoc_total;
  /*Looping over every lag*/
  for(curr_assoc=n_assoc; curr_assoc<max_assoc; curr_assoc += n_p1, seeds++) {
    if(!clock_seed)
      rk_seed(*seeds, &state);
    max_total = totalp + n_surr;
    /*Looping over every desired surrogate for this lag*/
    for(; totalp<max_total; totalp++) {
      *totalp = 0;
      /*Sample with replacement the original distribution, to the same size,
       *and calculate the total of this resampling
       */
      for(sample=0; sample<n_p1; sample++)
	*totalp += curr_assoc[rk_interval(max, &state)];
    }
  }
}
