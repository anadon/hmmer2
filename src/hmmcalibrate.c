/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2006 HHMI Janelia Farm
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* hmmcalibrate.c
 *
 * Score an HMM against random sequence data sets;
 * set histogram fitting parameters.
 */

#include "config.h"    /* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>

#include "squid.h"    /* general sequence analysis library    */
#include "stopwatch.h"    /* process timings                      */
#include "structs.h"    /* data structures, macros, #define's   */
#include "funcs.h"    /* function declarations                */
#include "globals.h"    /* alphabet global variables            */


static char banner[] = "hmmcalibrate -- calibrate HMM search statistics";

static char usage[] = "\
Usage: hmmcalibrate [-options] <hmmfile>\n\
Available options are:\n\
  -h             : print short usage and version info, then exit\n\
";

static char experts[] = "\
  --cpu <n>      : run <n> threads in parallel (if threaded)\n\
  --fixed <n>    : fix random sequence length at <n>\n\
  --histfile <f> : save histogram(s) to file <f>\n\
  --mean <x>     : set random seq length mean at <x> [350]\n\
  --num <n>      : set number of sampled seqs to <n> [5000]\n\
  --sd <x>       : set random seq length std. dev to <x> [350]\n\
  --seed <n>     : set random seed to <n> [time()]\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",         TRUE,  sqdARG_NONE  },
  { "--cpu",      FALSE, sqdARG_INT },
  { "--fixed",    FALSE, sqdARG_INT   },
  { "--histfile", FALSE, sqdARG_STRING },
  { "--mean",     FALSE, sqdARG_FLOAT },
  { "--num",      FALSE, sqdARG_INT   },
  { "--sd",       FALSE, sqdARG_FLOAT },
  { "--seed",     FALSE, sqdARG_INT},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static void main_loop_serial(struct plan7_s *hmm, int seed, int nsample,
                             float lenmean, float lensd, int fixedlen,
                             struct histogram_s **ret_hist, float *ret_max);


/* A structure of this type is shared by worker threads in the POSIX
 * threads parallel version.
 */
struct workpool_s {
  /* Static configuration:
   */
  struct plan7_s  *hmm;    /* ptr to single HMM to search with    */
  int    fixedlen;    /* if >0, fix random seq len to this   */
  float  lenmean;    /* mean of Gaussian for random seq len */
  float  lensd;      /* s.d. of Gaussian for random seq len */
  float *randomseq;             /* 0..Alphabet_size-1 i.i.d. probs     */
  int    nsample;    /* number of random seqs to do         */

  /* Shared (mutex-protected) input:
   */
  int    nseq;      /* current number of seqs searched     */

  /* Shared (mutex-protected) output:
   */
  struct histogram_s *hist;     /* histogram          */
  float          max_score;     /* maximum score seen */
  Stopwatch_t    watch;    /* Timings accumulated for threads */

  /* Thread pool information:
   */
  pthread_t      *thread;       /* our pool of threads */
  int             num_threads;  /* number of threads   */
  pthread_mutex_t input_lock;  /* a mutex protecting input fields */
  pthread_mutex_t output_lock;  /* a mutex protecting output fields */
};
static void main_loop_threaded(struct plan7_s *hmm, int seed, int nsample,
                               float lenmean, float lensd, int fixedlen,
                               int nthreads,
                               struct histogram_s **ret_hist, float *ret_max,
                               Stopwatch_t *twatch);
static struct workpool_s *workpool_start(struct plan7_s *hmm,
    float lenmean, float lensd, int fixedlen,
    float *randomseq, int nsample,
    struct histogram_s *hist,
    int num_threads);
static void  workpool_stop(struct workpool_s *wpool);
static void  workpool_free(struct workpool_s *wpool);
static void *worker_thread(void *ptr);


int
main(int argc, char **argv) {
  char    *hmmfile;             /* HMM file to open                */
  char    *tmpfile;             /* temporary calibrated HMM file   */
  HMMFILE *hmmfp;               /* opened hmm file pointer         */
  FILE    *outfp;               /* for writing HMM(s) into tmpfile */
  char    *mode;                /* write mode, "w" or "wb"         */
  struct plan7_s *hmm;          /* the hidden Markov model         */
  int     idx;      /* counter over sequences          */
  sigset_t blocksigs;    /* list of signals to protect from */
  int     nhmm;      /* number of HMMs calibrated       */

  struct histogram_s *hist = NULL; /* a resulting histogram           */
  float   max = 0;    /* maximum score from an HMM       */
  char   *histfile;             /* histogram save file             */
  FILE   *hfp;                  /* open file pointer for histfile  */

  Stopwatch_t stopwatch;  /* main stopwatch for process      */
  Stopwatch_t extrawatch;  /* stopwatch for threads           */

  float  *mu;      /* array of EVD mu's for HMMs      */
  float  *lambda;    /* array of EVD lambda's for HMMs  */
  int     mu_lumpsize;    /* allocation lumpsize for mu, lambda */

  int     nsample;    /* number of random seqs to sample */
  int     seed;      /* random number seed              */
  int     fixedlen;    /* fixed length, or 0 if unused    */
  float   lenmean;    /* mean of length distribution     */
  float   lensd;    /* std dev of length distribution  */


  char *optname;    /* name of option found by Getopt() */
  char *optarg;      /* argument found by Getopt()       */
  int   optind;            /* index in argv[]                  */

  int   num_threads;            /* number of worker threads */


  /***********************************************
   * Parse the command line
   ***********************************************/
  StopwatchStart(&stopwatch);
  StopwatchZero(&extrawatch);

  nsample      = 5000;
  fixedlen     = 0;
  lenmean      = 325.;
  lensd        = 200.;
  seed         = (int) time ((time_t *) NULL);
  histfile     = NULL;
  mu_lumpsize  = 100;
  num_threads  = sysconf(_SC_NPROCESSORS_ONLN);

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg)) {
    if      (strcmp(optname, "--cpu")      == 0) num_threads  = atoi(optarg);
    else if (strcmp(optname, "--fixed")    == 0) fixedlen = atoi(optarg);
    else if (strcmp(optname, "--histfile") == 0) histfile = optarg;
    else if (strcmp(optname, "--mean")     == 0) lenmean  = atof(optarg);
    else if (strcmp(optname, "--num")      == 0) nsample  = atoi(optarg);
    else if (strcmp(optname, "--sd")       == 0) lensd    = atof(optarg);
    else if (strcmp(optname, "--seed")     == 0) seed     = atoi(optarg);
    else if (strcmp(optname, "-h") == 0) {
      HMMERBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(0);
    }
  }

  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  hmmfile = argv[optind++];


  /***********************************************
   * Open our i/o file pointers, make sure all is well
   ***********************************************/

  /* HMM file */
  if ((hmmfp = HMMFileOpen(hmmfile, NULL)) == NULL)
    Die("failed to open HMM file %s for reading.", hmmfile);

  /* histogram file */
  hfp = NULL;
  if (histfile != NULL) {
    if ((hfp = fopen(histfile, "w")) == NULL)
      Die("Failed to open histogram save file %s for writing\n", histfile);
  }

  /* Generate calibrated HMM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * HMM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original HMM
   * file remains uncorrupted. tmpnam() doesn't work portably here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  tmpfile = MallocOrDie(strlen(hmmfile) + 5);
  strcpy(tmpfile, hmmfile);
  strcat(tmpfile, ".xxx");  /* could be more inventive here... */
  if (FileExists(tmpfile))
    Die("temporary file %s already exists; please delete it first", tmpfile);
  if (hmmfp->is_binary) mode = "wb";
  else                  mode = "w";

  /***********************************************
   * Show the banner
   ***********************************************/

  HMMERBanner(stdout, banner);
  printf("HMM file:                 %s\n", hmmfile);
  if (fixedlen)
    printf("Length fixed to:          %d\n", fixedlen);
  else {
    printf("Length distribution mean: %.0f\n", lenmean);
    printf("Length distribution s.d.: %.0f\n", lensd);
  }
  printf("Number of samples:        %d\n", nsample);
  printf("random seed:              %d\n", seed);
  printf("histogram(s) saved to:    %s\n",
         histfile != NULL ? histfile : "[not saved]");
  if (num_threads > 0)
    printf("POSIX threads:            %d\n", num_threads);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /***********************************************
   * Read the HMMs one at a time, and send them off
   * in probability form to one of the main loops.
   * The main loop functions are responsible for
   * synthesizing random sequences and returning
   * a score histogram for each HMM.
   ***********************************************/

  nhmm = 0;
  mu     = MallocOrDie(sizeof(float) * mu_lumpsize);
  lambda = MallocOrDie(sizeof(float) * mu_lumpsize);

  while (HMMFileRead(hmmfp, &hmm)) {
    if (hmm == NULL)
      Die("HMM file may be corrupt or in incorrect format; parse failed");

    if (num_threads <= 1)
      main_loop_serial(hmm, seed, nsample, lenmean, lensd, fixedlen,
                       &hist, &max);
    else if (num_threads > 0)
      main_loop_threaded(hmm, seed, nsample, lenmean, lensd, fixedlen,
                         num_threads, &hist, &max, &extrawatch);
    else
      Die("wait. that can't happen. I didn't do anything.");


    /* Fit an EVD to the observed histogram.
     * The TRUE left-censors and fits only the right slope of the histogram.
     * The 9999. is an arbitrary high number that means we won't trim
     * outliers on the right.
     */
    if (! ExtremeValueFitHistogram(hist, TRUE, 9999.))
      Die("fit failed; --num may be set too small?\n");

    mu[nhmm]     = hist->param[EVD_MU];
    lambda[nhmm] = hist->param[EVD_LAMBDA];
    nhmm++;
    if (nhmm % 100 == 0) {
      mu     = ReallocOrDie(mu,     sizeof(float) * (nhmm+mu_lumpsize));
      lambda = ReallocOrDie(lambda, sizeof(float) * (nhmm+mu_lumpsize));
    }

    /* Output
     */
    printf("HMM    : %s\n",   hmm->name);
    printf("mu     : %12f\n", hist->param[EVD_MU]);
    printf("lambda : %12f\n", hist->param[EVD_LAMBDA]);
    printf("max    : %12f\n", max);
    printf("//\n");

    if (hfp != NULL) {
      fprintf(hfp, "HMM: %s\n", hmm->name);
      PrintASCIIHistogram(hfp, hist);
      fprintf(hfp, "//\n");
    }

    FreeHistogram(hist);
    FreePlan7(hmm);
  }
  SQD_DPRINTF1(("Main body believes it has calibrations for %d HMMs\n", nhmm));

  /*****************************************************************
   * Rewind the HMM file for a second pass.
   * Write a temporary HMM file with new mu, lambda values in it
   *****************************************************************/

  HMMFileRewind(hmmfp);
  if (FileExists(tmpfile))
    Die("Ouch. Temporary file %s appeared during the run.", tmpfile);
  if ((outfp = fopen(tmpfile, mode)) == NULL)
    Die("Ouch. Temporary file %s couldn't be opened for writing.", tmpfile);

  for (idx = 0; idx < nhmm; idx++) {
    /* Sanity checks
     */
    if (!HMMFileRead(hmmfp, &hmm))
      Die("Ran out of HMMs too early in pass 2");
    if (hmm == NULL)
      Die("HMM file %s was corrupted? Parse failed in pass 2", hmmfile);

    /* Put results in HMM
     */
    hmm->mu     = mu[idx];
    hmm->lambda = lambda[idx];
    hmm->flags |= PLAN7_STATS;
    Plan7ComlogAppend(hmm, argc, argv);

    /* Save HMM to tmpfile
     */
    if (hmmfp->is_binary) WriteBinHMM(outfp, hmm);
    else                  WriteAscHMM(outfp, hmm);

    FreePlan7(hmm);
  }

  /*****************************************************************
   * Now, carefully remove original file and replace it
   * with the tmpfile. Note the protection from signals;
   * we wouldn't want a user to ctrl-C just as we've deleted
   * their HMM file but before the new one is moved.
   *****************************************************************/

  HMMFileClose(hmmfp);
  if (fclose(outfp)   != 0) PANIC;

  if (sigemptyset(&blocksigs) != 0) PANIC;
  if (sigaddset(&blocksigs, SIGINT) != 0) PANIC;
  if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   PANIC;
  if (remove(hmmfile) != 0)                            PANIC;
  if (rename(tmpfile, hmmfile) != 0)                   PANIC;
  if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) PANIC;
  free(tmpfile);

  /***********************************************
   * Exit
   ***********************************************/

  StopwatchStop(&stopwatch);

  free(mu);
  free(lambda);
  if (hfp != NULL) fclose(hfp);
  SqdClean();
  return 0;
}

/* Function: main_loop_serial()
 * Date:     SRE, Tue Aug 18 16:18:28 1998 [St. Louis]
 *
 * Purpose:  Given an HMM and parameters for synthesizing random
 *           sequences; return a histogram of scores.
 *           (Serial version)
 *
 * Args:     hmm      - an HMM to calibrate.
 *           seed     - random number seed
 *           nsample  - number of seqs to synthesize
 *           lenmean  - mean length of random sequence
 *           lensd    - std dev of random seq length
 *           fixedlen - if nonzero, override lenmean, always this len
 *           ret_hist - RETURN: the score histogram
 *           ret_max  - RETURN: highest score seen in simulation
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_serial(struct plan7_s *hmm, int seed, int nsample,
                 float lenmean, float lensd, int fixedlen,
                 struct histogram_s **ret_hist, float *ret_max) {
  struct histogram_s *hist;
  struct dpmatrix_s  *mx;
  float  randomseq[MAXABET];
  float  p1;
  float  max;
  int    sqlen;

  /* Initialize.
   * We assume we've already set the alphabet (safe, because
   * HMM input sets the alphabet).
   */
  sre_srandom(seed);
  P7Logoddsify(hmm, TRUE);
  P7DefaultNullModel(randomseq, &p1);
  hist = AllocHistogram(-200, 200, 100);
  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);
  max = -FLT_MAX;

  for (int idx = 0; idx < nsample; idx++) {
    float  score;
    /* choose length of random sequence */
    if (fixedlen) sqlen = fixedlen;
    else do sqlen = (int) Gaussrandom(lenmean, lensd);
      while (sqlen < 1);
    /* generate it */
    char  *seq;
    unsigned char  *dsq;
    seq = RandomSequence(Alphabet, randomseq, Alphabet_size, sqlen);
    dsq = DigitizeSequence(seq, sqlen);

#ifdef ALTIVEC
    /* We only need the score here (not the trace), so we can just
     * call the fast Altivec routine directly. The memory needs in this
     * routine is only proportional to the model length (hmm->M), and
     * preallocated, so we don't have to consider the low-memory alternative.
     */
    score = P7ViterbiNoTrace(dsq, sqlen, hmm, mx);
#else
    if (P7ViterbiSpaceOK(sqlen, hmm->M, mx))
      score = P7Viterbi(dsq, sqlen, hmm, mx, NULL);
    else
      score = P7SmallViterbi(dsq, sqlen, hmm, mx, NULL);
#endif

    AddToHistogram(hist, score);
    if (score > max) max = score;

    free(dsq);
    free(seq);
  }

  FreePlan7Matrix(mx);
  *ret_hist   = hist;
  *ret_max    = max;
  return;
}


/* Function: main_loop_threaded()
 * Date:     SRE, Wed Dec  1 12:43:09 1999 [St. Louis]
 *
 * Purpose:  Given an HMM and parameters for synthesizing random
 *           sequences; return a histogram of scores.
 *           (Threaded version.)
 *
 * Args:     hmm      - an HMM to calibrate.
 *           seed     - random number seed
 *           nsample  - number of seqs to synthesize
 *           lenmean  - mean length of random sequence
 *           lensd    - std dev of random seq length
 *           fixedlen - if nonzero, override lenmean, always this len
 *           nthreads - number of threads to start
 *           ret_hist - RETURN: the score histogram
 *           ret_max  - RETURN: highest score seen in simulation
 *           twatch   - RETURN: accumulation of thread times
 *
 * Returns:  (void)
 *           hist is alloc'ed here, and must be free'd by caller.
 */
static void
main_loop_threaded(struct plan7_s *hmm, int seed, int nsample,
                   float lenmean, float lensd, int fixedlen,
                   int nthreads,
                   struct histogram_s **ret_hist, float *ret_max,
                   Stopwatch_t *twatch) {
  struct histogram_s *hist;
  float  randomseq[MAXABET];
  float  p1;
  struct workpool_s *wpool;     /* pool of worker threads  */

  /* Initialize.
   * We assume we've already set the alphabet (safe, because
   * HMM input sets the alphabet).
   */
  sre_srandom(seed);
  P7Logoddsify(hmm, TRUE);
  P7DefaultNullModel(randomseq, &p1);
  hist = AllocHistogram(-200, 200, 100);

  wpool = workpool_start(hmm, lenmean, lensd, fixedlen, randomseq, nsample,
                         hist, nthreads);
  workpool_stop(wpool);

  *ret_hist = hist;
  *ret_max  = wpool->max_score;
  StopwatchInclude(twatch, &(wpool->watch));

  workpool_free(wpool);
  return;
}

/*****************************************************************
 * POSIX threads implementation.
 * API:
 *      workpool_start()   (makes a workpool_s structure. Starts calculations.)
 *      workpool_stop()    (waits for threads to finish.)
 *      [process histogram]
 *      workpool_free()    (destroys the structure)
 *
 * Threads:
 *      worker_thread()    (the actual parallelized worker thread).
 *****************************************************************/

/* Function: workpool_start()
 * Date:     SRE, Thu Jul 16 11:09:05 1998 [St. Louis]
 *
 * Purpose:  Initialize a workpool_s structure, and return it.
 *
 * Args:     hmm      - the HMM to calibrate
 *           fixedlen - 0, or a fixed length for seqs (bypass of Gaussian)
 *           lenmean  - mean sequence length
 *           lensd    - std. dev. for sequence length
 *           randomseq- i.i.d. frequencies for residues, 0..Alphabet_size-1
 *           nsample  - how many seqs to calibrate on
 *           hist     - histogram structure for storing results
 *           num_threads - how many processors to run on
 *
 * Returns:  ptr to struct workpool_s.
 *           Caller must wait for threads to finish with workpool_stop(),
 *           then free the structure with workpool_free().
 */
static struct workpool_s *
workpool_start(struct plan7_s *hmm, float lenmean, float lensd, int fixedlen,
               float *randomseq, int nsample, struct histogram_s *hist,
               int num_threads) {
  struct workpool_s *wpool;
  pthread_attr_t    attr;
  int i;
  int rtn;

  wpool         = MallocOrDie(sizeof(struct workpool_s));
  wpool->thread = MallocOrDie(num_threads * sizeof(pthread_t));
  wpool->hmm        = hmm;
  wpool->fixedlen   = fixedlen;
  wpool->lenmean    = lenmean;
  wpool->lensd      = lensd;
  wpool->randomseq  = randomseq;
  wpool->nsample    = nsample;

  wpool->nseq       = 0;
  wpool->hist       = hist;
  wpool->max_score  = -FLT_MAX;
  wpool->num_threads= num_threads;

  StopwatchZero(&(wpool->watch));

  if ((rtn = pthread_mutex_init(&(wpool->input_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));
  if ((rtn = pthread_mutex_init(&(wpool->output_lock), NULL)) != 0)
    Die("pthread_mutex_init FAILED; %s\n", strerror(rtn));

  /* Create slave threads.
   */
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  for (i = 0; i < num_threads; i++)
    if ((rtn = pthread_create(&(wpool->thread[i]), &attr,
                              worker_thread, (void *) wpool)) != 0)
      Die("Failed to create thread %d; return code %d\n", i, rtn);

  pthread_attr_destroy(&attr);

  return wpool;
}

/* Function: workpool_stop()
 * Date:     SRE, Thu Jul 16 11:20:16 1998 [St. Louis]
 *
 * Purpose:  Waits for threads in a workpool to finish.
 *
 * Args:     wpool -- ptr to the workpool structure
 *
 * Returns:  (void)
 */
static void
workpool_stop(struct workpool_s *wpool) {
  int i;
  /* wait for threads to stop */
  for (i = 0; i < wpool->num_threads; i++)
    if (pthread_join(wpool->thread[i],NULL) != 0)
      Die("pthread_join failed");
  return;
}

/* Function: workpool_free()
 * Date:     SRE, Thu Jul 16 11:26:27 1998 [St. Louis]
 *
 * Purpose:  Free a workpool_s structure, after the threads
 *           have finished.
 *
 * Args:     wpool -- ptr to the workpool.
 *
 * Returns:  (void)
 */
static void
workpool_free(struct workpool_s *wpool) {
  free(wpool->thread);
  free(wpool);
  return;
}

/* Function: worker_thread()
 * Date:     SRE, Thu Jul 16 10:41:02 1998 [St. Louis]
 *
 * Purpose:  The procedure executed by the worker threads.
 *
 * Args:     ptr  - (void *) that is recast to a pointer to
 *                  the workpool.
 *
 * Returns:  (void *)
 */
void *
worker_thread(void *ptr) {
  struct plan7_s    *hmm;
  struct dpmatrix_s *mx;
  struct workpool_s *wpool;
  int         len;
  int         rtn;
  Stopwatch_t thread_watch;

  StopwatchStart(&thread_watch);
  wpool = (struct workpool_s *) ptr;
  hmm   = wpool->hmm;
  mx    = CreatePlan7Matrix(1, hmm->M, 25, 0);
  for (;;) {
    /* 1. Synthesize a random sequence.
     *    The input sequence number is a shared resource,
     *    and sre_random() isn't thread-safe, so protect
     *    the whole section with mutex.
     */
    /* acquire a lock */
    if ((rtn = pthread_mutex_lock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    /* generate a sequence */
    wpool->nseq++;
    if (wpool->nseq > wpool->nsample) {
      /* we're done; release input lock, break loop */
      if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
        Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
      break;
    }
    if (wpool->fixedlen) len = wpool->fixedlen;
    else do len = (int) Gaussrandom(wpool->lenmean, wpool->lensd);
      while (len < 1);
    char *seq;
    seq = RandomSequence(Alphabet, wpool->randomseq, Alphabet_size, len);

    /* release the lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->input_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

    /* 2. Score the sequence against the model.
     */
    unsigned char *dsq;
    dsq = DigitizeSequence(seq, len);

    float       sc;
#ifdef ALTIVEC
    /* We only need the score here (not the trace), so we can just
     * call the fast Altivec routine directly. The memory needs in this
     * routine is only proportional to the modem length (hmm->M), and
     * preallocated, so we don't have to consider the low-memory alternative.
     */
    sc = P7ViterbiNoTrace(dsq, len, hmm, mx);
#else
    if (P7ViterbiSpaceOK(len, hmm->M, mx))
      sc = P7Viterbi(dsq, len, hmm, mx, NULL);
    else
      sc = P7SmallViterbi(dsq, len, hmm, mx, NULL);
#endif

    free(dsq);
    free(seq);

    /* 3. Save the output; hist and max_score are shared,
     *    so protect this section with the output mutex.
     */
    /* acquire lock on the output queue */
    if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
    /* save output */
    AddToHistogram(wpool->hist, sc);
    if (sc > wpool->max_score) wpool->max_score = sc;
    /* release our lock */
    if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
      Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));
  }

  StopwatchStop(&thread_watch);
  /* acquire lock on the output queue */
  if ((rtn = pthread_mutex_lock(&(wpool->output_lock))) != 0)
    Die("pthread_mutex_lock failure: %s\n", strerror(rtn));
  /* accumulate cpu time into main stopwatch */
  StopwatchInclude(&(wpool->watch), &thread_watch);
  /* release our lock */
  if ((rtn = pthread_mutex_unlock(&(wpool->output_lock))) != 0)
    Die("pthread_mutex_unlock failure: %s\n", strerror(rtn));

  FreePlan7Matrix(mx);
  pthread_exit(NULL);
  return NULL; /* solely to silence compiler warnings */
}
