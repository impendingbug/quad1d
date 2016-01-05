/**
 *  \file  cr_util.c
 *  \brief  Extension of integration/util.c from GSL source code, to accomodate complex-valued functions with
 *          real-valued free parameters.
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

static inline
void cr_update_re(qag_workspace *workspace_re,
                  qag_workspace *workspace_im,
                  double a1, double b1, double complex area1, double complex error1,
                  double a2, double b2, double complex area2, double complex error2);

static inline
void cr_update_im(qag_workspace *workspace_re,
                  qag_workspace *workspace_im,
                  double a1, double b1, double complex area1, double complex error1,
                  double a2, double b2, double complex area2, double complex error2);

static inline
void cr_update_1st(qag_workspace *workspace_re,
                   qag_workspace *workspace_im,
                   double a1, double b1, double complex area1, double complex error1,
                   double a2, double b2, double complex area2, double complex error2);

static inline
void cr_update_2nd(qag_workspace *workspace_re,
                   qag_workspace *workspace_im,
                   double a1, double b1, double complex area1, double complex error1,
                   double a2, double b2, double complex area2, double complex error2);


static inline void
cr_retrieve_re(const qag_workspace *workspace_re,
               const qag_workspace *workspace_im,
               double *a, double *b,
               double complex *r, double complex *e);


static inline void
cr_retrieve_im(const qag_workspace *workspace_re,
               const qag_workspace *workspace_im,
               double *a, double *b,
               double complex *r, double complex *e);

static inline
void cr_update_re(qag_workspace *workspace_re,
                  qag_workspace *workspace_im,
                  double a1, double b1, double complex area1, double complex error1,
                  double a2, double b2, double complex area2, double complex error2)
{
    double *alist_re = workspace_re->alist;
    double *blist_re = workspace_re->blist;
    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;
    size_t *level_re = workspace_re->level;
    const size_t i_split = workspace_re->i;   //the interval to bisect is based on workspace_re
    const size_t i_new = workspace_re->size;
    const size_t new_level = workspace_re->level[i_split] + 1;

    double *alist_im = workspace_im->alist;
    double *blist_im = workspace_im->blist;
    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;
    size_t *level_im = workspace_im->level;
    size_t *order_im = workspace_im->order;
    size_t nrmax_im;

    if (i_new > 1)
        nrmax_im = workspace_im->nrmax;
    else
        nrmax_im = 0;   //At the first subdivision, i_split is order_im[0] = 0 instead of order_im[nrmax_im]

    if (order_im[nrmax_im] != i_split)       //find nrsplit_im and then pull it through the order list to nrmax_im
    {                                  //nrsplit_m is the position of the split interval in the order_im list
		int nrsplit_im = 0;

		while (order_im[nrsplit_im] != i_split)
			nrsplit_im++;

		if (nrsplit_im > nrmax_im)
		{
			do {
				order_im[nrsplit_im] = order_im[nrsplit_im - 1];
				nrsplit_im--;
			} while (nrsplit_im > nrmax_im);
		}
		else
		{
			do {
				order_im[nrsplit_im] = order_im[nrsplit_im + 1];
				nrsplit_im++;
			} while (nrsplit_im < nrmax_im);
		}
	}
	//now we have inserted an empty cell to the "beginning" of the order_im list; the position of the split interval
	//within the list not yet defined

  /* append the newly-created intervals to the list */
	if (creal(error2) > creal(error1))
	{
		rlist_re[i_split] = creal(area2);
		elist_re[i_split] = creal(error2);
		level_im[i_split] = level_re[i_split] = new_level;
		alist_im[i_split] = alist_re[i_split] = a2;        /* blist[i_split] is already == b2 */
		rlist_im[i_split] = cimag(area2);
		elist_im[i_split] = cimag(error2);

		rlist_im[i_new] = cimag(area1);
		elist_im[i_new] = cimag(error1);
		level_re[i_new] = level_im[i_new] = new_level;
		alist_re[i_new] = alist_im[i_new] = a1;
		blist_re[i_new] = blist_im[i_new] = b1;
		rlist_re[i_new] = creal(area1);
		elist_re[i_new] = creal(error1);

		if (cimag(error2) > cimag(error1))   //make order[nrmax_im] point to the bigger error and order[i_new] to the smaller
		{
			order_im[nrmax_im] = i_split;
			order_im[i_new] = i_new;
		}
		else
		{
			order_im[nrmax_im] = i_new;
			order_im[i_new] = i_split;
		}
	}

	else
	{
		rlist_re[i_split] = creal(area1);
		elist_re[i_split] = creal(error1);
		level_im[i_split] = level_re[i_split] = new_level;
		blist_im[i_split] = blist_re[i_split] = b1;        /* alist[i_split] is already == a1 */
		rlist_im[i_split] = cimag(area1);
		elist_im[i_split] = cimag(error1);

		rlist_im[i_new] = cimag(area2);
		elist_im[i_new] = cimag(error2);
		level_re[i_new] = level_im[i_new] = new_level;
		alist_re[i_new] = alist_im[i_new] = a2;
		blist_re[i_new] = blist_im[i_new] = b2;
		rlist_re[i_new] = creal(area2);
		elist_re[i_new] = creal(error2);

		if (cimag(error2) > cimag(error1))   //make order[nrmax_im] point to the bigger error and order[i_new] to the smaller
		{
			order_im[nrmax_im] = i_new;
			order_im[i_new] = i_split;
		}
		else
		{
			order_im[nrmax_im] = i_split;
			order_im[i_new] = i_new;
		}
	}

  if (new_level > workspace_re->maximum_level)
	{
		workspace_re->maximum_level = new_level;
		workspace_im->maximum_level = new_level;
	}
  workspace_re->size++;
  cr_qpsrt_normal(workspace_re);
  workspace_im->size++;
  cr_qpsrt_special(workspace_im);
}

static inline
void cr_update_im(qag_workspace *workspace_re,
                  qag_workspace *workspace_im,
                  double a1, double b1, double complex area1, double complex error1,
                  double a2, double b2, double complex area2, double complex error2)
{
    double *alist_im = workspace_im->alist;
    double *blist_im = workspace_im->blist;
    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;
    size_t *level_im = workspace_im->level;
    const size_t i_split = workspace_im->i;   //the interval to bisect is based on workspace_im
    const size_t i_new = workspace_im->size;
    const size_t new_level = workspace_im->level[i_split] + 1;

    double *alist_re = workspace_re->alist;
    double *blist_re = workspace_re->blist;
    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;
    size_t *level_re = workspace_re->level;
	size_t *order_re = workspace_re->order;
	size_t nrmax_re;
	if (i_new > 1)
		nrmax_re = workspace_re->nrmax;
	else
		nrmax_re = 0;   //At the first subdivision, i_split is order_re[0] instead of order_re[nrmax_re]

	if (order_re[nrmax_re] != i_split)       //find nrsplit_re and then pull it through the order list to nrmax_re
	{                                        //nrsplit_m is the position of the split interval in the order_re list
		int nrsplit_re = 0;
		while (order_re[nrsplit_re] != i_split)
			nrsplit_re++;
		if (nrsplit_re > nrmax_re)
		{
			do {
				order_re[nrsplit_re] = order_re[nrsplit_re - 1];
				nrsplit_re--;
			} while (nrsplit_re > nrmax_re);
		}
		else
		{
			do {
				order_re[nrsplit_re] = order_re[nrsplit_re + 1];
				nrsplit_re++;
			} while (nrsplit_re < nrmax_re);
		}
	}
	//now we have inserted an empty cell to the "beginning" of the order_re list; the position of the split interval
	//within the list not yet defined

  /* append the newly-created intervals to the list */
	if (cimag(error2) > cimag(error1))
	{
		rlist_re[i_split] = creal(area2);
		elist_re[i_split] = creal(error2);
		level_im[i_split] = level_re[i_split] = new_level;
		alist_im[i_split] = alist_re[i_split] = a2;        /* blist[i_split] is already == b2 */
		rlist_im[i_split] = cimag(area2);
		elist_im[i_split] = cimag(error2);

		rlist_im[i_new] = cimag(area1);
		elist_im[i_new] = cimag(error1);
		level_re[i_new] = level_im[i_new] = new_level;
		alist_re[i_new] = alist_im[i_new] = a1;
		blist_re[i_new] = blist_im[i_new] = b1;
		rlist_re[i_new] = creal(area1);
		elist_re[i_new] = creal(error1);

        //make order[nrmax_re] point to the bigger error and order[i_new] to the smaller
		if (creal(error2) > creal(error1))
		{
			order_re[nrmax_re] = i_split;
			order_re[i_new] = i_new;
		}
		else
		{
			order_re[nrmax_re] = i_new;
			order_re[i_new] = i_split;
		}
	}
	else
	{
		rlist_re[i_split] = creal(area1);
		elist_re[i_split] = creal(error1);
		level_im[i_split] = level_re[i_split] = new_level;
		blist_im[i_split] = blist_re[i_split] = b1;        /* alist[i_split] is already == a1 */
		rlist_im[i_split] = cimag(area1);
		elist_im[i_split] = cimag(error1);

		rlist_im[i_new] = cimag(area2);
		elist_im[i_new] = cimag(error2);
		level_re[i_new] = level_im[i_new] = new_level;
		alist_re[i_new] = alist_im[i_new] = a2;
		blist_re[i_new] = blist_im[i_new] = b2;
		rlist_re[i_new] = creal(area2);
		elist_re[i_new] = creal(error2);

        //make order[nrmax_re] point to the bigger error and order[i_new] to the smaller
		if (creal(error2) > creal(error1))
		{
			order_re[nrmax_re] = i_new;
			order_re[i_new] = i_split;
		}
		else
		{
			order_re[nrmax_re] = i_split;
			order_re[i_new] = i_new;
		}
	}

  if (new_level > workspace_re->maximum_level)
	{
		workspace_im->maximum_level = new_level;
		workspace_re->maximum_level = new_level;
	}

    workspace_re->size++;
    cr_qpsrt_special(workspace_re);
    workspace_im->size++;
    cr_qpsrt_normal(workspace_im);
}

static inline
void cr_update_1st(qag_workspace *workspace_re,
                   qag_workspace *workspace_im,
                   double a1, double b1, double complex area1, double complex error1,
                   double a2, double b2, double complex area2, double complex error2)
{
    double *alist_re = workspace_re->alist;
    double *blist_re = workspace_re->blist;
    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;
    size_t *level_re = workspace_re->level;
    const size_t i_split = workspace_re->i;   //the interval to bisect is based on workspace_re
    const size_t i_new = workspace_re->size;
    const size_t new_level = workspace_re->level[i_split] + 1;

    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;

  /* append the newly-created intervals to the list */
	if (creal(error2) > creal(error1))
	{
		rlist_re[i_split] = creal(area2);
		elist_re[i_split] = creal(error2);
		level_re[i_split] = new_level;
		alist_re[i_split] = a2;        /* blist[i_split] is already == b2 */
		rlist_im[i_split] = cimag(area2);
		elist_im[i_split] = cimag(error2);

		rlist_im[i_new] = cimag(area1);
		elist_im[i_new] = cimag(error1);
		level_re[i_new] = new_level;
		alist_re[i_new] = a1;
		blist_re[i_new] = b1;
		rlist_re[i_new] = creal(area1);
		elist_re[i_new] = creal(error1);
	}

    else
    {
        rlist_re[i_split] = creal(area1);
        elist_re[i_split] = creal(error1);
        level_re[i_split] = new_level;
        blist_re[i_split] = b1;        /* alist[i_split] is already == a1 */
        rlist_im[i_split] = cimag(area1);
        elist_im[i_split] = cimag(error1);

        rlist_im[i_new] = cimag(area2);
        elist_im[i_new] = cimag(error2);
        level_re[i_new] = new_level;
        alist_re[i_new] = a2;
        blist_re[i_new] = b2;
        rlist_re[i_new] = creal(area2);
        elist_re[i_new] = creal(error2);
    }

    if (new_level > workspace_re->maximum_level)
        workspace_re->maximum_level = new_level;
    workspace_re->size++;
    qpsrt(workspace_re);
    workspace_im->size++;			//imaginary part is not sorted
}

static inline
void cr_update_2nd(qag_workspace *workspace_re,
                   qag_workspace *workspace_im,
                   double a1, double b1, double complex area1, double complex error1,
                   double a2, double b2, double complex area2, double complex error2)
{
    double *alist_im = workspace_im->alist;
    double *blist_im = workspace_im->blist;
    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;
    size_t *level_im = workspace_im->level;
    const size_t i_split = workspace_im->i;   //the interval to bisect is based on workspace_im
    const size_t i_new = workspace_im->size;
    const size_t new_level = workspace_im->level[i_split] + 1;

    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;

  /* append the newly-created intervals to the list */
	if (cimag(error2) > cimag(error1))
	{
		rlist_im[i_split] = cimag(area2);
		elist_im[i_split] = cimag(error2);
		level_im[i_split] = new_level;
		alist_im[i_split] = a2;        /* blist[i_split] is already == b2 */
		rlist_re[i_split] = creal(area2);
		elist_re[i_split] = creal(error2);

		rlist_re[i_new] = creal(area1);
		elist_re[i_new] = creal(error1);
		level_im[i_new] = new_level;
		alist_im[i_new] = a1;
		blist_im[i_new] = b1;
		rlist_im[i_new] = cimag(area1);
		elist_im[i_new] = cimag(error1);
	}

	else
	{
        rlist_im[i_split] = cimag(area1);
        elist_im[i_split] = cimag(error1);
        level_im[i_split] = new_level;
        blist_im[i_split] = b1;        /* alist[i_split] is already == a1 */
        rlist_re[i_split] = creal(area1);
        elist_re[i_split] = creal(error1);

        rlist_re[i_new] = creal(area2);
        elist_re[i_new] = creal(error2);
        level_im[i_new] = new_level;
        alist_im[i_new] = a2;
        blist_im[i_new] = b2;
        rlist_im[i_new] = cimag(area2);
        elist_im[i_new] = cimag(error2);
	}

    if (new_level > workspace_im->maximum_level)
        workspace_im->maximum_level = new_level;
    workspace_im->size++;
    qpsrt(workspace_im);
    workspace_re->size++;		//real part is not sorted
}

static inline void
cr_retrieve_re(const qag_workspace *workspace_re,
               const qag_workspace *workspace_im,
               double *a, double *b,
               double complex *r, double complex *e)
{
    const size_t i = workspace_re->i;
    double *alist = workspace_re->alist;
    double *blist = workspace_re->blist;
    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;
    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;

    *a = alist[i];
    *b = blist[i];

    *r = rlist_re[i] + I * rlist_im[i];
    *e = elist_re[i] + I * elist_im[i];
}

static inline void
cr_retrieve_im(const qag_workspace *workspace_re,
               const qag_workspace *workspace_im,
               double *a, double *b,
               double complex *r, double complex *e)
{
    const size_t i = workspace_im->i;
    double *alist = workspace_im->alist;
    double *blist = workspace_im->blist;
    double *rlist_im = workspace_im->rlist;
    double *elist_im = workspace_im->elist;
    double *rlist_re = workspace_re->rlist;
    double *elist_re = workspace_re->elist;

    *a = alist[i];
    *b = blist[i];

    *r = rlist_re[i] + I * rlist_im[i];
    *e = elist_re[i] + I * elist_im[i];
}

static inline double
sum_results(const qag_workspace * workspace);

static inline double
sum_results(const qag_workspace * workspace)
{
    const double *const rlist = workspace->rlist ;
    const size_t n = workspace->size;

    double result_sum = 0;

    for (int k = 0; k < n; k++)
        result_sum += rlist[k];

    return result_sum;
}

static inline int
subinterval_too_small(double a1, double a2, double b2);

static inline int
subinterval_too_small(double a1, double a2, double b2)
{
    const double e = DBL_EPSILON;
    const double u = DBL_MIN;

    double tmp = (1 + 100 * e) * (fabs (a2) + 1000 * u);

    int status = fabs (a1) <= tmp && fabs (b2) <= tmp;

    return status;
}

