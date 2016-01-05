/**
 *  \file  qag_workspace.c
 *  \brief  Allocate and free workspace for QAG adaptive quadrature from GSL library.
 *
 *  Based on integration/workspace.c of the source code of GSL library.
 *
 *  \internal
 *       Created:  04-08-14
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

#include "cr_quad1d.h"
#include <gsl/gsl_errno.h>

#include <stdlib.h>

#define RETURN_IF_NULL(x) if (x == NULL) { return ; }

qag_workspace *
qag_workspace_alloc(const size_t n)
{
  qag_workspace *w;

  if (n == 0)
	{
		GSL_ERROR_VAL("workspace length n must be positive integer",
									GSL_EDOM, NULL);
	}

  w = malloc(sizeof *w);

  if (w == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for workspace struct",
									GSL_ENOMEM, NULL);
	}

  w->alist = malloc(n * sizeof *w->alist);

  if (w->alist == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for alist ranges",
									GSL_ENOMEM, NULL);
	}

  w->blist = malloc(n * sizeof *w->blist);

  if (w->blist == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for blist ranges",
									GSL_ENOMEM, NULL);
	}

  w->rlist = malloc(n * sizeof *w->rlist);

  if (w->rlist == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for rlist ranges",
									GSL_ENOMEM, NULL);
	}


  w->elist = malloc(n * sizeof *w->elist);

  if (w->elist == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for elist ranges",
									GSL_ENOMEM, NULL);
	}

  w->order = malloc(n * sizeof *w->order);

  if (w->order == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for order ranges",
									GSL_ENOMEM, NULL);
	}

  w->level = malloc(n * sizeof *w->level);

  if (w->level == NULL)
	{
		GSL_ERROR_VAL("failed to allocate space for level ranges",
									GSL_ENOMEM, NULL);
	}

  w->size = 0;
  w->limit = n;
  w->maximum_level = 0;

  return w;
}

void
qag_workspace_free(const qag_workspace *w)
{
  RETURN_IF_NULL(w);
  free (w->level);
  free (w->order);
  free (w->elist);
  free (w->rlist);
  free (w->blist);
  free (w->alist);
  free ((void *)w);
}
