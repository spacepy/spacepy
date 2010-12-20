/*JTN 20101213
 *Library providing C support to poppy
 */

#include <stdio.h>
#ifdef THREADED
#include <pthread.h>
#include <string.h>
#endif

#include "randomkit.h"

/*Returns index of leftmost value in a list which is greater than or equal
 *to the comparison value.
 *If entire list is less than the comparison, returns length of list.
 *list: set of inputs to compare to
 *n: number of values in list
 *value: the value to compare list against
 *Largely based on implementation in Python 2.6.6
 */
long bisect_left(double* list, long n, double value)
{
  long lo = 0;
  long hi = n;
  long mid;
  while(lo < hi) {
    mid = (lo + hi) / 2;
    if (list[mid] < value)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

/*Returns index of righttmost value in a list which is less than or equal
 *to the comparison value.
 *If entire list is less than the comparison, returns length of list.
 *list: set of inputs to compare to
 *n: number of values in list
 *value: the value to compare list against
 *Largely based on implementation in Python 2.6.6
 */
long bisect_right(double* list, long n, double value)
{
  long lo = 0;
  long hi = n;
  long mid;
  while(lo < hi) {
    mid = (lo + hi) / 2;
    if (value < list[mid])
      hi = mid;
    else
      lo = mid + 1;
  }
  return lo;
}

/*Performs inner loop of association analysis
 *p2: list of times in process 2 (count: n_p2)
 *starts: list of times in p1 at the beginning of the window (n_p1)
 *stops: times in p1 at the end of the window (n_p1)
 *lags: lag times between p1 and p2 (n_lags)
 *n_assoc: 2D arry of resulting association numbers (n_lags, n_p1)
 */
void assoc(double* p2, double* starts, double* stops, double* lags,
	   long* n_assoc, long n_p2, long n_p1, long n_lags)
{
  long l_idx, i;
  long* curr_assoc;
  for(l_idx = 0; l_idx < n_lags; l_idx++) {
    curr_assoc = n_assoc + (l_idx * n_p1);
    for(i = 0; i < n_p1; i++) {
      curr_assoc[i] = bisect_right(p2, n_p2, stops[i] + lags[l_idx]) - 
	bisect_left(p2, n_p2, starts[i] + lags[l_idx]);
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
  last_surr = surrogates + (n * n_els);
  for(curr_surr=surrogates; curr_surr < last_surr; curr_surr++)
    *curr_surr = data[rk_interval(max, &state)];
}

/*Perfoms inner loop of association analysis confidence interval, including
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
 *seed: number to seed RNG
 *clock_seed: if true, seed RNG from the clock, ignore seed
 */
void aa_ci(unsigned long *n_assoc, unsigned long *surr_assoc_total,
	   unsigned long n_lags, unsigned long n_p1, unsigned long n_surr,
	   unsigned long seed, int clock_seed)
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
  else
    rk_seed(seed, &state);

  max_assoc = n_assoc + n_p1 * n_lags;
  totalp = surr_assoc_total;
  /*Looping over every lag*/
  for(curr_assoc=n_assoc; curr_assoc<max_assoc; curr_assoc += n_p1) {
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

#ifdef THREADED
#define MAXTHREADS 32

/*Structure to pass information to a subjob of aa_ci*/
struct aa_ci_data {
  unsigned long *n_assoc;
  unsigned long *surr_assoc_total;
  unsigned long n_lags;
  unsigned long n_p1;
  unsigned long n_surr;
  rk_state state; /*State of the RNG*/
};

/*A subjob of aa_ci
 *Parameters passed via threadarg, pointer to aa_ci_data
 *Pointers (and equivalent counters) passed in should be updated to include
 *only the slice of the total job that this subjob should be doing.
 */
void *aa_ci_subjob(void *threadarg)
{
  struct aa_ci_data *jobdata;
  unsigned long sample; /*Sample number*/
  unsigned long *curr_assoc; /*Association numbers for the current lag*/
  unsigned long *max_assoc; /*Just past last value of association numbers*/
  unsigned long *totalp; /*Point to total assoc. no.; current lag, surrogate*/
  unsigned long *max_total; /*Just past last value of total for this lag*/
  unsigned long max; /*Maximum index from RNG*/

  jobdata = (struct aa_ci_data *) threadarg;
  max = jobdata->n_p1 - 1;
  max_assoc = jobdata->n_assoc + jobdata->n_p1 * jobdata->n_lags;
  totalp = jobdata->surr_assoc_total;
  /*Looping over every lag*/
  for(curr_assoc=jobdata->n_assoc; curr_assoc<max_assoc;
      curr_assoc += jobdata->n_p1) {
    max_total = totalp + jobdata->n_surr;
    /*Looping over every desired surrogate for this lag*/
    for(; totalp<max_total; totalp++) {
      *totalp = 0;
      /*Sample with replacement the original distribution, to the same size,
       *and calculate the total of this resampling
       */
      for(sample=0; sample<jobdata->n_p1; sample++)
	*totalp += curr_assoc[rk_interval(max, &(jobdata->state))];
    }
  }
  return 0;
}

/*Multithreaded version: lags are split across threads.
 *Perfoms inner loop of association analysis confidence interval, including
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
 *seed: number to seed RNG
 *clock_seed: if true, seed RNG from the clock, ignore seed
 *            NB: If clock_seed is NOT specified (i.e. seed is used),
 *            threads will run SEQUENTIALLY (not concurrently),
 *            to ensure determinism.
 *threadcount: number of worker threads to spawn
 */
void aa_ci_threaded(unsigned long *n_assoc, unsigned long *surr_assoc_total,
		    unsigned long n_lags, unsigned long n_p1,
		    unsigned long n_surr,
		    unsigned long seed, int clock_seed, int threadcount)
{
  int threadno, rc;
  void *status;
  unsigned long offsets[MAXTHREADS+1]; /*Start/end offset for each thread*/
  double count; /*Number of lags handled by each thread*/
  struct aa_ci_data jobdata[MAXTHREADS];
  pthread_t threads[MAXTHREADS];
  pthread_attr_t attr;

  if(threadcount > MAXTHREADS)
    threadcount = MAXTHREADS;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  count = (double)(n_lags) / threadcount;
  offsets[0] = 0;
  offsets[threadcount] = n_lags;
  for(threadno=1; threadno<threadcount; threadno++)
    offsets[threadno] = (unsigned long)((double)(offsets[threadno-1]) + count);
  for(threadno=0; threadno<threadcount; threadno++) {
    jobdata[threadno].n_assoc = n_assoc + (offsets[threadno] * n_p1);
    jobdata[threadno].surr_assoc_total = surr_assoc_total +
      (offsets[threadno] * n_surr);
    jobdata[threadno].n_lags = offsets[threadno + 1] - offsets[threadno];
    jobdata[threadno].n_p1 = n_p1;
    jobdata[threadno].n_surr = n_surr;
    if(clock_seed)
      rk_randomseed(&jobdata[threadno].state);
    else {
      if(threadno > 0)
	memcpy(&(jobdata[threadno].state), &jobdata[threadno - 1].state,
	       sizeof(rk_state));
      else
	rk_seed(seed, &jobdata[threadno].state);
    }
    rc = pthread_create(&threads[threadno], &attr, aa_ci_subjob,
			(void*)&jobdata[threadno]);
    if(rc) {
      fprintf(stderr, "Unable to make thread: %d\n", rc);
      return;
    }

    if(!clock_seed) {
      rc = pthread_join(threads[threadno], &status);
      if(rc) {
	fprintf(stderr, "Unable to join thread: %d\n", rc);
	return;
      }
    }
  }
  if(clock_seed) {
    for(threadno=0; threadno<threadcount; threadno++) {
      rc = pthread_join(threads[threadno], &status);
      if(rc)
	fprintf(stderr, "Unable to join thread: %d\n", rc);
    }
  }
}
#endif
