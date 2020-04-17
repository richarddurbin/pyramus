/*  File: pyramus.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2020
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 17 01:57 2020 (rd109)
 * Created: Wed Apr  8 07:55:58 2020 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"

typedef struct {
  char *seq ;
  int len ;
  int K, K4 ;           // kmer size K and 4^K
  U32 *ki ;		// size len: kmer at position i 
  int *kc ;             // size 2^K: count of kmer k
  int *ko ;             // size 2^K: sum of kc[0..(k-1)]
  int *ks ;             // size len: location in sequence of i'th sorted kmer
} Reference ;

typedef struct {
  int start, end ;
  int n ;
  int off ;
} Interval ;

static int intervalCompareN (const void *a, const void *b) // sort on decreasing n
{ return ((Interval*)b)->n - ((Interval*)a)->n ; }

static int intervalCompareStart (const void *a, const void *b) // sort on increasing start
{ return ((Interval*)a)->start - ((Interval*)b)->start ; }

/* GLOBALS */

static int K = 8 ;

/*****************************************************/
// Kmer package: main iteration is inlined for speed.  Usage:
//   Kmer *k = kmerCreate (s, n, K) ;
//   while (k->n) { U32 x = kmerNext (k) ; <do something with x> }
//   kmerDestroy (k) ;
// NB does not copy the sequence - requires it to stay present during use

typedef struct {
  char *s ;
  int n ;			/* remaining bases: decreases as read */
  U32 x ;
  U32 mask ;
} Kmer ;

Kmer *kmerCreate (char *s, int len, int K)
{
  if (len < K) die ("kmerCreate: len %d , K %d\n", len, K) ;
  Kmer *km = new0 (1, Kmer) ;
  km->mask = (1 << (2*K)) - 1 ;
  --K ;
  km->n = len - K ;
  while (K--) km->x = (km->x << 2) + *s++ ;
  km->s = s ;
  return km ;
}

static inline U32 kmerNext (Kmer *km)
{ km->x = ((km->x << 2) + *km->s++) & km->mask ;
  --km->n ;
  return km->x ;
}

void kmerDestroy (Kmer *km) { free (km) ; }

/*****************************************************/
/* utility to reverse complement a sequence */

void reverseComplement (char *s, int len)
{
  int i = 0, j = len-1 ;
  while (i <= j)
    { char t = s[j] ; s[j] = 3 - s[i] ; s[i] = 3 - t ; // a little wasteful if i == j, but OK
      ++i ; --j ;
    }
}

/*****************************************************/

void examine1 (Reference *ref)
{
  if (!ref) die ("must load a reference before --examine 1") ;

  int i, *his = new0(ref->len, int) ;
  for (i = 0 ; i < ref->K4 ; ++i) ++his[ref->kc[i]] ;
  for (i = 1 ; i < ref->len ; ++i)
    if (his[i]) printf ("  %d kmers with count %d\n", his[i], i) ;
  free (his) ;
}

void examine2 (Reference *ref, SeqIO *rio) // show where reads match
{
  if (!ref) die ("must load a reference before --examine 2") ;
  if (!rio) die ("must call --reads before --examine 2") ;
  
  int i ;
  int *loc1 = new(150,int), *loc2 = new(150,int) ;
  int *n1 = new0(150,int), *n2 = new0(150,int) ;
  for (i = 0 ; i < ref->len-K ; ++i)
    if (ref->kc[ref->ki[i]] == 1) ++n1[i/200] ;
    else ++n2[i/200] ;
  int nSeq = 0 ;
  while (seqIOread (rio))
    { bzero (loc1, 150*sizeof(int)) ; bzero (loc2, 150*sizeof(int)) ;
      Kmer *km = kmerCreate (sqioSeq(rio), rio->seqLen, K) ;
      for (i = 0 ; km->n ; ++i)
	{ U32 k = kmerNext (km) ;
	  int c = ref->kc[k] ;
	  int *s = &ref->ks[ref->ko[k]] ;
	  if (c == 1) ++loc1[*s/200] ;
	  else while (c--) ++loc2[*s++/200] ;
	}

      printf ("read length %llu\n", rio->seqLen) ;

      for (i = 0 ; i < 150 ; ++i)
	if (!loc1[i]) putchar ('.') ;
	else if (loc1[i] > 9) putchar ('*') ;
	else putchar ('0'+loc1[i]) ;
      for (i = 0 ; i < 150 ; ++i) if (loc1[i] > 9) printf (" %d/%d", loc1[i], n1[i]) ;
      putchar ('\n') ;
      for (i = 0 ; i < 150 ; ++i)
	if (!loc2[i]) putchar ('.') ;
	else if (loc2[i] > 9) putchar ('*') ;
	else putchar ('0'+loc2[i]) ;
      for (i = 0 ; i < 150 ; ++i) if (loc2[i] > 9) printf (" %d/%d", loc2[i], n2[i]) ;
      putchar ('\n') ;
      kmerDestroy (km) ;
      if (++nSeq > 100) exit (0) ;
    }
  free (loc1) ; free (loc2) ;
}

/*******************************************/

void findIntervals (Reference *ref, SeqIO *rio, int binSize)
{
  if (!ref) die ("must call --reference before --findIntervals") ;
  if (!rio) die ("must call --reads before --findIntervals") ;
  if (binSize <= 0) die ("must give binSize > 0 for --findIntervals") ;

  int i, j ;

  int nBin = ref->len/binSize ;
  int **counts = new(nBin,int*) ;
  for (i = 0 ; i < nBin ; ++i) counts[i] = new0(nBin,int) ;
  int *cumu = new(ref->len,int) ; // cumulative count of unique kmers in reference
  cumu[0] = ref->kc[ref->ki[0]] == 1 ? 1 : 0 ;
  for (i = 1 ; i < ref->len-K ; ++i)
    if (ref->kc[ref->ki[i]] == 1) cumu[i] = cumu[i-1]+1 ;
    else cumu[i] = cumu[i-1] ;
  while (i < ref->len) { cumu[i] = cumu[i-1]+1 ; ++i ; }

  Array aicF = arrayCreate (16, Interval) ;
  Array aicR = arrayCreate (16, Interval) ;
  int icBestF, icBestR ;
  int nSeq = 0 ;
  while (seqIOread (rio))
    { int icBest, strand = 0 ;
      Interval *ic ;
      while (strand < 2)
	{ if (strand == 1) reverseComplement (sqioSeq(rio), rio->seqLen) ;
	  Array aic = strand ? aicR : aicF ; arrayMax(aic) = 0 ; // restart candidate list
	  Kmer *km = kmerCreate (sqioSeq(rio), rio->seqLen, K) ;
	  for (i = 0 ; km->n ; ++i)
	    { U32 k = kmerNext (km) ;
	      int c = ref->kc[k] ;
	      int *s = &ref->ks[ref->ko[k]] ;
	      if (c == 1)
		{ if (arrayMax(aic))
		    { ic = arrp(aic,icBest,Interval) ;
		      int xDiff = *s - (ic->start + i + ic->off) ;
		      if (xDiff > -20 && xDiff < 20)
			{ ++ic->n ; ic->off += xDiff ; ic->end = *s + rio->seqLen - i ;
			  continue ;
			}
		      else
			for (j = 0 ; j < arrayMax(aic) ; ++j)
			  { ic = arrp(aic,j,Interval) ;
			    xDiff = *s - (ic->start + i + ic->off) ;
			    if (xDiff > -20 && xDiff < 20)
			      { ++ic->n ; ic->off += xDiff ; ic->end = *s + rio->seqLen - i ;
				if (ic->n > arrp(aic,icBest,Interval)->n)
				  icBest = j ;
				continue ;
			      }
			  }
		    }
		  // only get here if there was no candidate match - start a new candidate
		  ic = arrayp(aic,arrayMax(aic),Interval) ;
		  ic->start = *s - i ; ic->end = *s + rio->seqLen - i ;
		  ic->n = 1 ; ic->off = 0 ;
		  if (arrayMax(aic) == 1) icBest = 0 ; // this was the first candidate
		}
	    }
	  kmerDestroy (km) ;
	  if (strand) icBestR = icBest ; else icBestF = icBest ;
	  ++strand ;
	}

      if (arrp(aicF,icBestF,Interval)->n >= arrp(aicR,icBestR,Interval)->n)
	{ ic = arrp(aicF,icBestF,Interval) ; strand = 0 ; }
      else
	{ ic = arrp(aicR,icBestR,Interval) ; strand = 1 ; }
      if (ic->start < 0) ic->start = 0 ; else if (strand) ic->start += 15 ;
      if (ic->end >= ref->len) ic->end = ref->len - 1 ; else if (!strand) ic->end -= 15 ;
      int nTot = 0 ;
      for (j = 0 ; j < arrayMax(aicF) ; ++j) nTot += arrp(aicF,j,Interval)->n ;
      for (j = 0 ; j < arrayMax(aicR) ; ++j) nTot += arrp(aicR,j,Interval)->n ;
//        int nPoss = cumu[ic->end - K] - cumu[ic->start] ;
//	  printf ("read %d len %llu match %s %d - %d len %d n %d / %d off %d, else %d / %d\n",
//		  nSeq, rio->seqLen,
//		  strand ? "reverse" : "forward", ic->start, ic->end, (ic->end - ic->start),
//		  ic->n, nPoss, ic->off, nTot - ic->n, 2*cumu[ref->len-1] - nPoss) ;
      ++counts[ic->start/binSize][ic->end/binSize] ;
      ++nSeq ;
    }
  fprintf (stderr, "finished mapping endpoints of %d reads: ", nSeq) ; timeUpdate (stderr) ;

  /* Strategy is: for each start find the median end, and for each of these the median start.
     Then output pairs which are mutual medians.
  */
  Array a = arrayCreate (1024, Interval) ;
  for (i = nBin ; i-- ; )
    { int m = 0 ; for (j = i ; j < nBin ; ++j) m += counts[i][j] ;
      int m2 = 0 ; for (j = i ; m2 < (m+1)/2 && j < nBin ; ++j) m2 += counts[i][j] ; 
      int jj ; for (jj = i ; jj < nBin ; ++jj) counts[i][jj] = 0 ;
      --j ; counts[i][j] = m ;
    }
  for (j = 0 ; j < nBin ; ++j)
    { int m = 0 ; for (i = j+1 ; i-- ; ) m += counts[i][j] ;
      if (!m) continue ;
      int m2 = 0 ; for (i = j+1 ; i-- && m2 < (m+1)/2 ; ) m2 += counts[i][j] ;
      Interval *ic = arrayp (a, arrayMax(a), Interval) ;
      ic->start = ++i ; ic->end = j ; ic->n = m ;
    }
  /* now greedily select a covering by eliminating overlaps from selected intervals */
  arraySort (a, intervalCompareN) ;
  Interval *ic, *jc ;
  BOOL *cover = new0(nBin,BOOL) ;
  for (i = 0, ic = arrp(a,i,Interval) ; i < arrayMax(a) ; ++i, ++ic)
    if (ic->n)
      { for (j = ic->start ; j <= ic->end ; ++j) cover[j] = TRUE ;
	for (j = i+1, jc = arrp(a,j,Interval) ; j < arrayMax(a) ; ++j, ++jc)
	  if (jc->n)
	    { if (jc->start >= ic->start-50/binSize && jc->end <= ic->end+50/binSize) // NB hardcoded 50
		{ ic->n += jc->n ; jc->n = 0 ; }
	      BOOL covered = TRUE ;
	      int k ; for (k = jc->start ; k <= jc->end ; ++k) covered &= cover[k] ;
	      if (covered) { ic->n += jc->n ; jc->n = 0 ; }
	    }
      }
  free (cover) ;
  arraySort (a, intervalCompareStart) ;

  int nTile = 0, lastEnd = 0 ;
  for (i = 0, ic = arrp(a,i,Interval) ; i < arrayMax(a) ; ++i, ++ic)
    if (ic->n)
      { printf (" start %d end %d length %d overlap %d count %d\n", ic->start*binSize, ic->end*binSize,
		(ic->end-ic->start)*binSize, (lastEnd-ic->start)*binSize, ic->n) ;
	lastEnd = ic->end ;
	++nTile ;
      }
  printf ("selected %d intervals from %d candidates from %d reads\n",
	  nTile, arrayMax(a), nSeq) ;
  arrayDestroy (a) ;

  arrayDestroy (aicF) ; arrayDestroy (aicR) ; free (cumu) ;
  for (i = 0 ; i < ref->len/binSize ; ++i) free (counts[i]) ; free (counts) ;
  timeUpdate (stderr) ;
}

/*******************************************/

void usage (void)
{ fprintf (stderr, "Usage: pyramus <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them\n") ;
  fprintf (stderr, "  -R | --reference <refFileName>\n") ;
  fprintf (stderr, "  -r | --reads <readFileName>       // must call before -I, -x\n") ;
  fprintf (stderr, "  -I | --findIntervals <binSize>\n") ;
  fprintf (stderr, "  -K | --Kmer <K> [%d]\n", K) ;
  fprintf (stderr, "  -x | --explore <number>           // various exploratory functions\n") ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */
  if (!argc) usage () ;

  timeUpdate (stderr) ;		/* initialise timer */

  int i ;			/* generically useful variable */

  Reference *ref = 0 ;		/* key state-holding variables */
  SeqIO *rio = 0 ;

  while (argc > 0)
    { if (**argv != '-')
	die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
      fprintf (stderr, "COMMAND %s", *argv) ;
      for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
      fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
      if (ARGMATCH("-R","--reference",2))
	{ SeqIO *si = seqIOopenRead (argv[-1], dna2indexConv, FALSE) ;
	  if (!si) die ("failed to open reference sequence %s", argv[-1]) ;
	  if (!seqIOread (si)) die ("failed to read reference sequence %s", argv[-1]) ;
	  ref = new0 (1, Reference) ;
	  ref->len = si->seqLen ;
	  ref->seq = new(ref->len,char) ;
	  memcpy (ref->seq, sqioSeq(si), ref->len) ;
	  seqIOclose (si) ;
	
	  ref->ki = new(ref->len,U32) ;	// kmer at position i
	  ref->K = K ; ref->K4 = 1 << (2*K) ;
	  ref->kc = new0(ref->K4,int) ;	// counts of each kmer
	  { Kmer *km = kmerCreate (ref->seq, ref->len, K) ;
	    for (i = 0 ; km->n ; ++i) { ref->ki[i] = kmerNext (km) ; ++ref->kc[ref->ki[i]] ; }
	    kmerDestroy (km) ;
	  }
	  ref->ko = new(ref->K4,int) ;	// offsets of location lists in ks
	  ref->ko[0] = 0 ;
	  for (i = 1 ; i < ref->K4 ; ++i) ref->ko[i] = ref->ko[i-1] + ref->kc[i-1] ;
	  ref->ks = new(ref->len,int) ;	// locations of kmer hits
	  bzero (ref->kc, ref->K4*sizeof(int)) ;
	  for (i = 0 ; i < ref->len - ref->K + 1 ; ++i)
	    { int k = ref->ki[i] ; ref->ks[ref->ko[k] + ref->kc[k]++] = i ; }
	  
	  printf ("K = %d, K4 = %d, reference sequence has %d bases\n", ref->K, ref->K4, ref->len) ;
	}
      else if (ARGMATCH("-r","--reads", 2))
	{ rio = seqIOopenRead (argv[-1], dna2indexConv, FALSE) ;
	  if (!rio) die ("failed to open reads file %s", argv[-1]) ;
	}
      else if (ARGMATCH("-I","--findIntervals",2))
	{ findIntervals (ref, rio, atoi(argv[-1])) ;
	  seqIOclose (rio) ; rio = 0 ;
	}
      else if (ARGMATCH("-x","--examine",2))
	switch (atoi (argv[-1]))
	  {
	  case 1:
	    examine1 (ref) ;
	    break ;
	    
	  case 2:
	    examine2 (ref, rio) ;
	    seqIOclose (rio) ; rio = 0 ;
	    break ;

	  default: die ("unknown examine option %s", argv[-1]) ;
	  }
      else
	die ("unknown command %s or not enough arguments", *argv) ;
    }

  if (ref)
    { free (ref->seq) ; free (ref->ki) ; free (ref->kc) ;
      free (ref->ko) ; free (ref->ks) ; free (ref) ;
    }
  if (rio)
    seqIOclose (rio) ;

  timeTotal (stderr) ;
  exit (0) ;
}
