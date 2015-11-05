/* poly/solve_quadratic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* solve_quadratic.c - finds the real roots of a x^2 + b x + c = 0 */

#include <math.h>


#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <string.h>



#ifndef _types_h
#define _types_h
#include "types.h"
#endif



#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif



#include "ld.h"



int 
gsl_poly_solve_quadratic (double a, double b, double c, 
                          double *x0, double *x1)
{
  double disc = b * b - 4 * a * c;

  if (a == 0) /* Handle linear case */
    {
      if (b == 0)
        {
          return 0;
        }
      else
        {
          *x0 = -c / b;
          return 1;
        };
    }

  if (disc > 0)
    {
      if (b == 0)
        {
          double r = fabs (0.5 * sqrt (disc) / a);
          *x0 = -r;
          *x1 =  r;
        }
      else
        {
          double sgnb = (b > 0 ? 1 : -1);
          double temp = -0.5 * (b + sgnb * sqrt (disc));
          double r1 = temp / a ;
          double r2 = c / temp ;

          if (r1 < r2) 
            {
              *x0 = r1 ;
              *x1 = r2 ;
            } 
          else 
            {
              *x0 = r2 ;
              *x1 = r1 ;
            }
        }
      return 2;
    }
  else if (disc == 0) 
    {
      *x0 = -0.5 * b / a ;
      *x1 = -0.5 * b / a ;
      return 2 ;
    }
  else
    {
      return 0;
    }
}








/* poly/solve_cubic.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */


#include <math.h>
/* "-std=c99" skips this constant! */
#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int 
gsl_poly_solve_cubic (double a, double b, double c, 
                      double *x0, double *x1, double *x2)
{
  double q = (a * a - 3 * b);
  double r = (2 * a * a * a - 9 * a * b + 27 * c);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;
  double xx=CR2 - CQ3; //anders
  double tol=0.000000001;
  if (R == 0 && Q == 0)
    {
      *x0 = - a / 3 ;
      *x1 = - a / 3 ;
      *x2 = - a / 3 ;
      return 3 ;
    }
  else if ((xx<tol&&xx>-tol)|(xx<-tol&&xx>tol)) //anders
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          *x0 = -2 * sqrtQ  - a / 3;
          *x1 = sqrtQ - a / 3;
          *x2 = sqrtQ - a / 3;
        }
      else
        {
          *x0 = - sqrtQ  - a / 3;
          *x1 = - sqrtQ - a / 3;
          *x2 = 2 * sqrtQ - a / 3;
        }
      return 3 ;
    }
  else if (CR2 < CQ3) /* equivalent to R2 < Q3 */
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      *x0 = norm * cos (theta / 3) - a / 3;
      *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
      *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;
      
      /* Sort *x0, *x1, *x2 into increasing order */

      if (*x0 > *x1)
        SWAP(*x0, *x1) ;
      
      if (*x1 > *x2)
        {
          SWAP(*x1, *x2) ;
          
          if (*x0 > *x1)
            SWAP(*x0, *x1) ;
        }
      
      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      *x0 = A + B - a / 3;
      return 1;
    }
}





#define my_warning if(0) printf

#define MIN(x,y) ((x) < (y) ? (x):(y))


//return the transposed iMatrix
iMatrix *getTransposed(iMatrix *mat){
  iMatrix *retMat = allocIntMatrix(mat->y,mat->x);
  for (int i=0;i<retMat->x;i++)
    for (int j=0;j<retMat->y;j++)
      retMat->matrix[i][j]=mat->matrix[j][i];
  return retMat;
}





static double *get_expt(geno_cptr res, double p)
{
  double u, v, w, x;
  double *expt =(double *)calloc(9, sizeof(double));
  
  u = res->j + p * res->e;
  v = res->k + (1.0 - p) * res->e;
  w = res->l + (1.0 - p) * res->e;
  x = res->m + p * res->e;
  
  expt[0] = u * u;
  expt[1] = 2. * u * v;
  expt[2] = v * v;
  
  expt[3] = 2. * u * w;
  expt[4] = 2. * u * x + 2. * v * w;
  expt[5] = 2. * v * x;
  
  expt[6] = w * w;
  expt[7] = 2. * w * x;
  expt[8] = x * x;

  return expt;
}



static void pick_best_p(geno_cptr res, double *p, int root_count)
{
  int i;
  my_warning("root count2 %d\n", root_count);
  for (i = 0; i < root_count; i++)
    {
      my_warning(" tryingg %f\n", p[i]);
      //      pick_best_p(res, x, 3);
      if ((p[i] > 1.0) || (p[i] < 0.0))
	{
	  my_warning("skiping %f\n", p[i]);
	}
      else
	{
	  int j;
	  double *expt = get_expt(res, p[i]);
	  double sum = 0;
	  for (j = 0; j < 9; j++)
	    {
	      if(res->count[j])
		{
		  if(expt[j] > 0.0)
		    {
		      sum += res->count[j] * log(expt[j]);
		    }
		  else
		    {
		      my_warning("why is this %f?\n", expt[j]);
		    }
		}
	    }
	  if (sum > res->post_p)
	    {
	      res->post_p = sum;
	      res->p = p[i];

	      if (res->expt)
		free(res->expt);
	      res->expt = expt;
	    }
	  else
	    {
	      free(expt);
	    }
	}
    }
}


geno_cptr do_geno_cal(geno_cptr res)
{
  int *count = res->count;
  int uv =0, uw =0, vx =0, wx = 0;
  int root_count = 3;
  double *expt =NULL;

  uv = res->j + res->k +res->e;
  wx = res->l + res->m +res->e;
  uw = res->j + res->l +res->e;
  vx = res->k + res->m +res->e;

  res->sign_of_r = 0;
  /* */
  if(res->e)
    {
      /* e != 0 */
      int tmp1 = res->e * (res->j + res->m - res->k - res->l);
      int e2 = res->e * res->e;
      int jm = (res->j) * (res->m);
      int kl = (res->k) * (res->l);
      double x[3];

      double a0 = - jm ;
      double a1 = e2 - tmp1 + jm + kl;
      double a2 = tmp1 - 3 * e2;
      double a3 = 2 * e2;
      
      if ((jm == 0) || (kl == 0))
	{
	  /* one of the roots are 0 or 1 */
	  if ((jm == 0) && (kl == 0))
	    {
	      x[0] = 0.0;
	      x[1] = 1.0;
	      x[2] = ((double) (res->e + res->k + res->l - res->j -res->m)) / res->e / 2.0 ;
	      /* haploview does something very different with this case */
	      my_warning("case 1\n");
	    }
	  else if (jm == 0)
	    {
	      x[0] = 0.0;
	      root_count = gsl_poly_solve_quadratic(a3,a2,a1,
						    &(x[1]), &(x[2])) + 1;
	      my_warning("case 2\n");
	    }
	  else
	    {
	      /* kl == 0, jm != 0, */
	      x[0] = 1.0;
	      root_count = gsl_poly_solve_quadratic(a3, (a2 - 2 * tmp1),
						     (a1 + 2 * tmp1),
						    &(x[1]), &(x[2])) + 1;
	      x[1] = 1.0 - x[1];
	      x[2] = 1.0 - x[2];
	      my_warning("case 3\n");
	    }

	}
      else 
	{ 
	  root_count = gsl_poly_solve_cubic (a2/a3, 
					     (a1)/a3,
					     a0/a3,
					     &(x[0]), 
					     &(x[1]), &(x[2]));
	  my_warning("case 4\n");
	}
	  pick_best_p(res, x, root_count);
    }
  else
    {
      /* e == 0 */
      double x[3];
      int bottom  = (res->j * res->m + res->k * res->l);
      my_warning("e = 0\n");
      if(bottom != 0)
	x[0] = (res->j * res->m)*1.0 /  bottom;
      else
	{
	  x[0] = 0;
	}
      /* just try the edge cases for completeness */
      //thorfinn start
      //x[0]=0;x[1]=1;x[2]=0;
      //x[2] = 1;
      //x[3] = 0;
      //thorfinn slut
      //following as a fix for snpMatrix 1.4.2
      x[1] = 1;
      x[2] = 0;
      
      pick_best_p(res, x, 3);
      /* haploview skip this case */
      my_warning("case 5\n");
    }

  my_warning("p = %f\n", res->p);

  res->u = res->j + res->p * res->e;
  res->v = res->k + (1.0 - res->p) * res->e;
  res->w = res->l + (1.0 - res->p) * res->e;
  res->x = res->m + res->p * res->e;
    
  expt = res->expt;
  if (res->total)
    {
      int jj;
      for (jj = 0; jj < 9; jj ++)
	expt[jj] /= (4. * res->total);
    }

  my_warning(" %i\t%i\t%i\t%f\t%f\t%f\n", count[0], count[1], count[2], expt[0], expt[1], expt[2]);
  my_warning(" %i\t%i\t%i\t%f\t%f\t%f\n", count[3], count[4], count[5], expt[3], expt[4], expt[5]);
  my_warning(" %i\t%i\t%i\t%f\t%f\t%f\n", count[6], count[7], count[8], expt[6], expt[7], expt[8]);
  
  if(res->total)
    {
      res->bigD = res->u * res->x - res->w * res->v;
      if (res->bigD > 0) 
	{
	  int bottom = MIN(vx *uv , wx * uw);
	  res->sign_of_r = +1;
	  if(bottom !=0)
	    res->dprime = res->bigD / bottom;
	  else
	    res->dprime = -1.0;
	}
      else    
	{
	  int bottom = MIN (uw * uv, wx * vx);
	  res->sign_of_r = -1;
	  if (bottom !=0)
	    res->dprime = - res->bigD / bottom ;
	  else
	    res->dprime = -1.0;
	}
      
      {
	/* integer over flow if multiplying all 4 first */
	if ((uw) && (uv) && (wx) && (vx))
	  res->rsq2 = (res->bigD /uw/uv) * (res->bigD /wx/vx);
	else
	  res->rsq2 = -1.0;
      }
    }
  else
    {
      res->rsq2   = -1.0;
      res->dprime = -1.0;
      res->bigD   =  0.0;
    }

  res->lod = 0.0;
  if(res->j) 
    {
      my_warning("%d %d %f %d %d\n", res->total , res->j, res->u , uw , uv);
      res->lod += res->j * log ( res->total * 2.0 * res->u / uw / uv);
    }
  if(res->k) 
    res->lod += res->k * log ( res->total * 2.0 * res->v / vx / uv);  
  if(res->l) 
    res->lod += res->l * log ( res->total * 2.0 * res->w / wx / uw);
  if(res->m) 
    res->lod += res->m * log ( res->total * 2.0 * res->x / wx / vx);
  if(res->e) 
    res->lod += res->e * log ( res->total * res->total * 2.0 /uv /uw * (res->u * res->x + res->v * res->w) / wx / vx ); /* *4 /2.0 */

  my_warning("%d %d %d %d %d, %d %d %d %d\n",
	     res->j, res->k, res->l, res->m, res->e, 
	     uv, uw , wx , vx);

#ifndef M_LOG10E
  /* /usr/include/math.h on linux if certain glibc extension macros are defined */ 
# define M_LOG10E   0.43429448190325182765  /* log_10 e */
#endif
  res->lod *= M_LOG10E;

  my_warning("%f %f %f %f\n",
	     res->u /res->total * 0.5, res->v/res->total * 0.5, res->w/res->total * 0.5, res->x/res->total * 0.5);

  my_warning("d\' = %f , r^2 = %f, lod= %f\n", res->dprime, res->rsq2, res->lod);
  return res;
}

geno_cptr get_geno_count(int *pos1, int *pos2, int length)
{
  int run = 0;
  int *count = NULL;
  geno_cptr res = (geno_cptr)calloc(1, sizeof(geno_count));
  res->post_p = - 1.0e308; /* large negative number */
  res->expt = NULL;  
  count = res->count;
  memset(count, 0, sizeof(count));

  while (run < length)
    
    {

      if ((*pos1) && (*pos2))
	{
	  count[ ((*pos1) -1) * 3 + (*pos2) -1] ++;
	}
      pos1++;
      pos2++;
      run++;
    }
  res->a = count[0];
  res->b = count[1];
  res->c = count[2];
  res->d = count[3];
  res->e = count[4];
  res->f = count[5];
  res->g = count[6];
  res->h = count[7];
  res->i = count[8];
  //  for(int i=0;i<9;i++)
  //printf(" %d \n",count[i]);

  res->total = count[0] + count[1] + count[2]
    + count[3] + count[4] + count[5]
    + count[6] + count[7] + count[8];

  /* */
  res->j = res->a * 2 + res->b + res->d;
  res->k = res->c * 2 + res->b + res->f;
  res->l = res->g * 2 + res->d + res->h;
  res->m = res->i * 2 + res->f + res->h;

 return do_geno_cal(res);
}
