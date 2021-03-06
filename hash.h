/*  File: hash.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *-------------------------------------------------------------------
 * Description:
      - can hash integers, floats, pointers
      - returns an integer to be used as index into a local array
      - index of first item to be stored is 1, second is 2 etc.
      - support removal of objects, with free-list to reuse indices of removed objects
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 24 07:17 2019 (rd109)
 * * Jan 24 07:01 2019 (rd109): restructured so hashAdd() and hashFind() return BOOL 
         and the index value is returned via a pointer, as for DICT.
 * Created: Thu Jan 13 10:57:20 2011 (rd)
 *-------------------------------------------------------------------
 */

#ifndef HASH_DEFINED
#define HASH_DEFINED

typedef void* HASH ;
typedef union {long int i ; float f ; void* p ;} HASHKEY ;

/* define keys so as to allow ourselves to hash value 0 */
static inline HASHKEY HASH_INT   (long int i) { HASHKEY hk ; hk.i = i^INT_MAX ; return hk ; }
static inline HASHKEY HASH_FLOAT (float f) { HASHKEY hk ; hk.f = f ; hk.i ^= INT_MAX ; return hk ; }
static inline HASHKEY HASH_PTR   (void *p) { HASHKEY hk ; hk.p = p ; return hk ; }

HASH hashCreate (int n) ;
void hashDestroy (HASH h) ;
void hashClear (HASH h) ;
BOOL hashAdd  (HASH h, HASHKEY k, int *index) ; /* return TRUE if added, fill index if non-zero */
BOOL hashFind (HASH h, HASHKEY k, int *index) ; /* return TRUE if found, fill index if non-zero */
BOOL hashRemove (HASH h, HASHKEY k) ; /* if found, remove and return TRUE, else FALSE */
int  hashCount (HASH h) ;	      /* number of keys stored in hash */

/* iterator to get all key-value pairs, in arbitrary order */
/* note that if nothing is removed then values are incrementing from 1 to n_added */
/* note also that adding while iterating can invalidate the iterator */
void hashInitIterator (HASH h) ;
int  hashNextKeyValue (HASH h, HASHKEY *kp, int *ip) ; /* returns value, 0 if done */

void hashStats (void) ;		/* overall stats on package/performance */

#endif

/*********** end of file ************/
