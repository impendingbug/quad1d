/**
 *  \file  initialise.c
 *  \brief  Copied from integration/initialise.c of the GSL source code.
 *
 *  <+DETAILED+>
 *
 *  \internal
 *       Created:  04-08-14
 *      Revision:  none
 *      Compiler:  gcc -std=c11 -pedantic
 *  Organization:  DINS, Utrecht
 *
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 *  This source code is released for free distribution under the terms of the
 *  GNU General Public License as published by the Free Software Foundation.
 */

static inline
void initialise(qag_workspace *workspace, double a, double b);

static inline
void initialise(qag_workspace *workspace, double a, double b)
{
  workspace->size = 0;
  workspace->nrmax = 0;
  workspace->i = 0;
  workspace->alist[0] = a;
  workspace->blist[0] = b;
  workspace->rlist[0] = 0.0;
  workspace->elist[0] = 0.0;
  workspace->order[0] = 0;
  workspace->level[0] = 0;

  workspace->maximum_level = 0;
}
