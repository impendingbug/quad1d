/**
 *  \file  set_initial.c
 *  \brief  Copied from integration/set_initial.c of the GSL source code.
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
void set_initial_result(qag_workspace *workspace,
                        double result, double error);

static inline
void set_initial_result(qag_workspace *workspace,
                        double result, double error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}
