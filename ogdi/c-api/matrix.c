/******************************************************************************
 *
 * Component: OGDI Core C API
 * Purpose: Double precision matrix inverse, allocation, multiplication,
 *          and print routines.
 * 
 ******************************************************************************
 */
#include <math.h>

#include "ecs.h"
#include "matrix.h"

#define EPSILON 1.0e-16

/*
 * initialize values within a matrix.
 */
void mat_init(double **matrix, int height, int width) 
{
  int i=0, j=0;
      for (i=0; i < height; i++) {
         for (j=0; j < width; j++) {
	    matrix[i][j]=0;
          } 
     }	
}

/**
 * mat_mult_direct - multiply directly two 2D matrixes
 *
 */
void mat_mul_direct (matrix_a, height_a, width_a, matrix_b, height_b, width_b, matrix_result)
    double **matrix_a; 
    int height_a; 
    int width_a;
    double **matrix_b; 
    int height_b; 
    int width_b;
    double **matrix_result; 
{

 /* initialize loop iterators */
 int i=0, j=0, k=0;
  
 /* compute multiply of each elements into destination matrix  */
 for ( i = 0 ; i< height_a; i++ )
    for ( j = 0 ; j< width_b; j++ )
      for ( k = 0 ; k< width_a; k++ )
	matrix_result[i][j] += matrix_a[i][k]*matrix_b[k][j];
  return;
}

/**
 * mat_mult_traspose - multiply traspose of first 2D matrix 
 *                     transposed with second 2D matrix
 */
void mat_mul_transposed (matrix_a, height_a, width_a, matrix_b, height_b, width_b, matrix_result)
    double **matrix_a; 
    int height_a; 
    int width_a;
    double **matrix_b; 
    int height_b; 
    int width_b;
    double **matrix_result; 
{
 /* initialize loop iterators */
  int i=0, j=0, k=0;

 /* compute multiply of each elements into destination matrix */
 for ( i = 0 ; i< height_a; i++ )
    for ( j = 0 ; j< width_b; j++ )
      for ( k = 0 ; k< width_a; k++ )
	matrix_result[i][j] += matrix_a[i][k]*matrix_b[j][k];
  return;
}

/*
 * inverse: invert a square matrix (puts pivot elements on main diagonal).
 *          returns arg2 as the inverse of arg1.
 *
 *  This routine is based on a routine found in Andrei Rogers, "Matrix
 *  Methods in Urban and Regional Analysis", (1971), pp. 143-153.
 */
int mat_inverse (double **matrix,int n)
{
    int res = -1;
    int i, j, k, l, ir=0, ic=0 ;
    int *ipivot = calloc(n, sizeof(int)), **itemp = calloc(n, 2*sizeof(int));
    double *pivot = calloc(n,sizeof(double)), t;
    double fabs();

    /* initialization */
    for (i = 0; i < n; i++)
        ipivot[i] = 0;

    for (i = 0; i < n; i++)
    {
        t = 0.0;  /* search for pivot element */

        for (j = 0; j < n; j++)
        {
            if (ipivot[j] == 1) /* found pivot */
                continue;

            for (k = 0; k < n; k++)
                switch (ipivot[k]-1)
                {
                    case  0:
                        break;
                    case -1:
                        if (fabs (t) < fabs (matrix[j][k]))
                        {
                            ir = j;
                            ic = k;
                            t = matrix[j][k];
                        }
                        break;
                    case  1:
		    	goto exit;
                        break;
                    default: /* shouldn't get here */
		    	goto exit;
                        break;
                }
        }

        ipivot[ic] += 1;
        if (ipivot[ic] > 1) /* check for dependency */
                {
            goto exit;
                }

        /* interchange rows to put pivot element on diagonal */
        if (ir != ic)
            for (l = 0; l < n; l++)
            {
                t = matrix[ir][l];
                matrix[ir][l] = matrix[ic][l];
                matrix[ic][l] = t;
            }

        itemp[i][0] = ir;
        itemp[i][1] = ic;
        pivot[i] = matrix[ic][ic];

        /* check for zero pivot */
        if (fabs (pivot[i]) < EPSILON)
                {
            goto exit;
                }

        /* divide pivot row by pivot element */
        matrix[ic][ic] = 1.0;

        for (j = 0; j < n; j++)
            matrix[ic][j] /= pivot[i];

        /* reduce nonpivot rows */
        for (k = 0; k < n; k++)
            if (k != ic)
            {
                t = matrix[k][ic];
                matrix[k][ic] = 0.0;

                for (l = 0; l < n; l++)
                    matrix[k][l] -= (matrix[ic][l] * t);
            }
    }

    /* interchange columns */
    for (i = 0; i < n; i++)
    {
        l = n - i - 1;
        if (itemp[l][0] == itemp[l][1])
            continue;

        ir = itemp[l][0];
        ic = itemp[l][1];

        for (k = 0; k < n; k++)
        {
            t = matrix[k][ir];
            matrix[k][ir] = matrix[k][ic];
            matrix[k][ic] = t;
        }
    }

    res = 1;

exit:
    free(ipivot);
    free(itemp);
    free(pivot);
    return res;
}

/*
 * allocate memory for a 2D matrix with type of double elements.
 * returns NULL on failure.
 */
double **mat_malloc(int height,int width)
{
  unsigned i=0;
  double **matrix=NULL;

    matrix = (double **)malloc(height*sizeof(double *));
	  if(matrix == NULL) {
	        printf("Allocating memory for matrix computation pointers failed. \n");
		free(matrix);
		exit(1);
		}

    for (i=0;i < width;i++) {
        matrix[i] = (double *)malloc(width * sizeof(double));
 	  if(matrix == NULL) {
	        printf("Allocating memory for matrix data failed. \n");
		free(matrix);
		exit(1);
	                     }
                    }
  /* initialize values in matrix */
  mat_init(matrix,height,width);

  /* matrix created, return it */
  return matrix;
}

/*
 * unallocate memory for a 2D matrix.
 */
void mat_free(double **matrix, int height) 
{
  int i=0;
  for (i=0; i<height; i++)
    free(matrix[i]);
  /* free root on witdh > 1 type matrix only.*/
  if (matrix)  
       free(matrix);
}
