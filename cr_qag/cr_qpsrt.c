/**
 *  \file  cr_qpsrt.c
 *  \brief  Extension of qpsrt.c to accomodate complex-valued functions with real-valued free parameters.
 *
 *  <+DETAILED+>
 *
 *  \author  Sandy Pratama
 *
 *  \internal
 *       Created:  04-08-14
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 2014 Sandy Pratama
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

static inline void
cr_qpsrt_normal(qag_workspace *workspace);

static inline
void cr_qpsrt_normal(qag_workspace *workspace)
//Used if i_split is picked based on the same workspace, in which case elist[last] always contains the smaller error.
{
  const size_t last = workspace->size - 1;		//qpsrt is called after size is updated
  //const size_t limit = workspace->limit;

  double *elist = workspace->elist;
  size_t *order = workspace->order;

  size_t i_nrmax = workspace->nrmax;
  const size_t i_maxerr = order[i_nrmax];

	if (last < 2)
	{
		order[0] = 0;
		order[1] = 1;
		workspace->i = i_maxerr;
		return;
	}

  const double errmax = elist[i_maxerr];

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]])
	{
		order[i_nrmax] = order[i_nrmax - 1];
		i_nrmax--;
	}

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
/*
  int top;
  if (last < (limit/2 + 2))
  	top = last;
  else
  	top = limit - last + 1;
*/

  const int top = last;

  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */

  int i = i_nrmax + 1;

  /* The order of the tests in the following line is important to
     prevent a segmentation fault */

  while (i < top && errmax < elist[order[i]])
	{
		order[i-1] = order[i];
		i++;
	}

  order[i-1] = i_maxerr;

  /* Insert errmin by traversing the list bottom-up */

  const double errmin = elist[last];

  int k = top - 1;

  while (k > i - 2 && errmin >= elist[order[k]])
	{
		order[k+1] = order[k];
		k--;
	}

  order[k+1] = last;

  /* Set i_max and e_max */

  workspace->i = order[i_nrmax];
  workspace->nrmax = i_nrmax;
}


static inline void
cr_qpsrt_special(qag_workspace *workspace);

static inline
void cr_qpsrt_special(qag_workspace *workspace)
//Used if i_split is picked based on the other workspace, in which case elist[last] does NOT always contains the
//smaller error.
{
  const size_t last = workspace->size - 1;		//qpsrt is called after size is updated
  //const size_t limit = workspace->limit;

  double *elist = workspace->elist;
  size_t *order = workspace->order;

  size_t i_nrmax = workspace->nrmax;
  const size_t i_maxerr = order[i_nrmax];

  if (last < 2)  //order[0] and order[1] are already assigned values
	{
		workspace->i = i_maxerr;
		return;
	}

  const double errmax = elist[i_maxerr];

  /* This part of the routine is only executed if, due to a difficult
     integrand, subdivision increased the error estimate. In the normal
     case the insert procedure should start after the nrmax-th largest
     error estimate. */

  while (i_nrmax > 0 && errmax > elist[order[i_nrmax - 1]])
	{
		order[i_nrmax] = order[i_nrmax - 1];
		i_nrmax--;
	}

  /* Compute the number of elements in the list to be maintained in
     descending order. This number depends on the number of
     subdivisions still allowed. */
  /*
  int top;
  if (last < (limit/2 + 2))
  	top = last;
  else
  	top = limit - last + 1;

  top = last;
  */

  const int top = last;

  /* Insert errmax by traversing the list top-down, starting
     comparison from the element elist(order(i_nrmax+1)). */

  int i = i_nrmax + 1;

  /* The order of the tests in the following line is important to
     prevent a segmentation fault */

  while (i < top && errmax < elist[order[i]])
	{
		order[i-1] = order[i];
		i++;
	}

  order[i-1] = i_maxerr;

  /* Insert errmin by traversing the list bottom-up */
  const int i_small_err = order[last];
  const double errmin = elist[i_small_err];

  int k = top - 1;

  while (k > i - 2 && errmin >= elist[order[k]])
	{
		order[k+1] = order[k];
		k--;
	}

  order[k+1] = i_small_err;

  /* Set i_max and e_max */

  workspace->i = order[i_nrmax];
  workspace->nrmax = i_nrmax;
}

