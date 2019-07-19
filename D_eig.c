/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: eigen.c,v 1.99 2013/01/06 15:30:39 zlb Exp $ */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define check phgPrintf("HELLO I AM HERE LINE = %d\n",__LINE__)

FLOAT alpha = 0.0;

FLOAT mass = 0.0;
/* analytic solution of the first eigen vector in the case of the unit cube */
static int nev_plus, (*modes)[3], *mode;


static void
func_phi(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if(x==0&&y==0&&z==0)
    {
       *value = 1.e12;
    }
    else
    {
       *value = -2/sqrt(x*x+y*y+z*z);
    }
} 

static void
func_A1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.0;
} 

static void
func_A2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.0;
} 

static void
func_A3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = 0.0;
} 

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Sin(mode[0] * M_PI * x) *
	     Sin(mode[1] * M_PI * y) *
	     Sin(mode[2] * M_PI * z);
} 

FLOAT
compute_error(DOF *u_h, int which)
{
    /* current range is [start, start + n) */
    static int start = -1, n = 0, *pvt = NULL;
    static FLOAT *A = NULL, *b;
    static DOF **v = NULL;
    DOF *u;
    INT i, k;

    if (u_h == NULL) {
	/* free buffers */
	if (v != NULL) {
	    for (i = 0; i < n; i++)
		phgDofFree(v + i);
	    FreeAtExit(A);
	    FreeAtExit(pvt);
	    FreeAtExit(v);
	    A = NULL;
	    pvt = NULL;
	    v = NULL;
	}
	start = -1;
	n = 0;
	return 0.;
    }

    mode = (int *)(modes + which);
    i = mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2];
    for (k = which; k >= 0; k--) {
	mode = (int *)(modes + k);
	if (i != mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2])
	    break;
    }
    if ((++k) != start) {
	if (start >= 0 && n > 0) {
	    assert(v != NULL);
	    for (start = 0; start < n; start++)
		phgDofFree(v + start);
	    phgFree(A);
	    phgFree(pvt);
	    phgFree(v);
	}
	start = k;
	for (n = which + 1; n < nev_plus; n++) {
	    mode = (int *)(modes + n);
	    if (i != mode[0] * mode[0] + mode[1] * mode[1] + mode[2] * mode[2])
		break;
	}
	n -= start;
	A = phgAlloc((n + 1) * n * sizeof(*A));
	b = A + n * n;
	pvt = phgAlloc(n * sizeof(*pvt));
	v = phgAlloc(n * sizeof(*v));
	for (i = 0; i < n; i++) {
	    mode = (int *)(modes + start + i);
	    v[i] = phgDofNew(u_h->g, u_h->type, 1, "v", func_u);
	}
	for (i = 0; i < n; i++) {
	    A[i * n + i] = phgDofDotL2Vec(v[i], v[i]);
	    for (k = i + 1; k < n; k++)
		A[i * n + k] = A[k * n + i] = phgDofDotL2Vec(v[i], v[k]);
	}
	phgSolverDenseLU(n, A, pvt);
    }

    for (i = 0; i < n; i++)
	b[i] = phgDofDotL2Vec(v[i], u_h);
    phgSolverDenseSV(n, A, pvt, 1, b);
    u = NULL;
    for (i = 0; i < n; i++)
	phgDofAXPY(b[i], v[i], &u);
    phgDofAXPY(-1.0, u_h, &u);
    b[0] = phgDofNormL2Vec(u) / phgDofNormL2Vec(u_h);
    phgDofFree(&u);

    return b[0];
}

FLOAT 
phgQuadBasAGradient1Bas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v, int m,
		       int order)
{
    int i;
    FLOAT d;
    const FLOAT *g1, *g2, *w;
    QUAD *quad;

//check;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    
	quad = phgQuadGetQuad3D(order);

//check;

    d = 0.;

    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;
    if (NULL == A) {
	for (i = 0; i < quad->npoints; i++,g1++,g2 += 3,w++) {
	    d += g1[0]*g2[0]*(*w);
/*	    d[1]+ = g1[0]*g2[1]*(*(w));*/
/*	    d[2]+ = g1[0]*g2[2]*(*(w));*/
	}
    }
    	  return d * phgGeomGetVolume(u->g, e);
}

FLOAT 
phgQuadBasAGradient2Bas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v, int m,
		       int order)
{
    int i;
    FLOAT d;
    const FLOAT *g1, *g2, *w;
    QUAD *quad;
//check;
    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    
	quad = phgQuadGetQuad3D(order);

    d = 0.;
//check;
    g1 = phgQuadGetBasisValues(e, u, n, quad);
//check;
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
//check;
    w = quad->weights;
    if (NULL == A) {
	for (i = 0; i < quad->npoints; i++,g1++,g2 +=3,w++) {
	    d += g1[0]*g2[1]*(*w);
	}
    }
    	    return d * phgGeomGetVolume(u->g, e);
}

FLOAT 
phgQuadBasAGradient3Bas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v, int m,
		       int order)
{
    int i;
    FLOAT d;
    const FLOAT *g1, *g2, *w;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    
	quad = phgQuadGetQuad3D(order);

    d = 0.;

    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;
    if (NULL == A) {
	for (i = 0; i < quad->npoints; i++,g1++,g2 +=3,w++) {
	    d += g1[0]*g2[2]*(*w);
	}
    }
    	    return d * phgGeomGetVolume(u->g, e);
}

static void
build_matrices(MAT *matA, MAT *matB, DOF *u_h1, DOF *phi, DOF *A1, DOF *A2, DOF *A3, FLOAT kx, FLOAT ky, FLOAT kz)
{
    int N1 = u_h1->type->nbas;	/* number of basis functions in an element */
    int i, j;
    GRID *g = u_h1->g;
    ELEMENT *e;
    int N = 8*N1;
    FLOAT A[N][N], B[N][N];
    INT I[N];

assert(u_h1->dim == 1);
ForAllElements(g, e){
	bzero(A, sizeof(A));
	bzero(B, sizeof(B));
//check;
/*	for(i=0;i<N1;i++)*/
/*		for(j=0;j<N1;j++)*/
/*			A[i][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass;*/
/*	for(i=0;i<N1;i++)*/
/*		for(j=0;j<N1;j++)*/
/*			A[i][j+3*N1] = -phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5);*/
/*	for(i=0;i<N1;i++)*/
/*		for(j=0;j<N1;j++)*/
/*			A[i][j+6*N1] = phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);*/
/*	for(i=0;i<N1;i++)*/
/*		for(j=0;j<N1;j++)*/
/*			A[i][j+7*N1] = phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5);*/
/*			*/
	for(i=0;i<N1;i++)
	{
	  for(j=0;j<N1;j++)
	    {
//check;
              A[i][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
              //+ alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      A[i][j+2*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz - phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i][j+3*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx-phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i][j+6*N1] = phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i][j+7*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky + phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      
	      
	      A[i+N1][j+N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //+ alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      A[i+N1][j+2*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx + phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+N1][j+3*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz + phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+N1][j+6*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky + phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5) + phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+N1][j+7*N1] = -phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      
	      
	      A[i+2*N1][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz - phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+2*N1][j+N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx - phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+2*N1][j+2*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      A[i+2*N1][j+4*N1] = phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+2*N1][j+5*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky + phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      
	      
	      A[i+3*N1][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx + phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+3*N1][j+N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz + phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+3*N1][j+3*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      A[i+3*N1][j+4*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky + phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      A[i+3*N1][j+5*N1] = -phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      
	      
	      A[i+4*N1][j+2*N1] = -phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+4*N1][j+3*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky - phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) + phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      A[i+4*N1][j+4*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      A[i+4*N1][j+6*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz - phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+4*N1][j+7*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx - phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      
	      
	      A[i+5*N1][j+2*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky - phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      A[i+5*N1][j+3*N1] = phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+5*N1][j+5*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //+ alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      A[i+5*N1][j+6*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx + phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+5*N1][j+7*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz + phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      
	      
	      A[i+6*N1][j] = -phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+6*N1][j+N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky - phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) + phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      A[i+6*N1][j+4*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz - phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+6*N1][j+5*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx - phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+6*N1][j+6*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //+ alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      
	      
	      A[i+7*N1][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*ky - phgQuadBasAGradient1Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A2,u_h1,j,5);
	      A[i+7*N1][j+N1] = phgQuadBasAGradient3Bas(e,u_h1,i,NULL,u_h1,j,5);
	      A[i+7*N1][j+4*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kx + phgQuadBasAGradient2Bas(e,u_h1,i,NULL,u_h1,j,5) - phgQuadBasABas(e,u_h1,i,A1,u_h1,j,5);
	      A[i+7*N1][j+4*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*kz + phgQuadBasABas(e,u_h1,i,A3,u_h1,j,5);
	      A[i+7*N1][j+7*N1] = -phgQuadBasDotBas(e,u_h1,i,u_h1,j,5)*mass + phgQuadBasABas(e,u_h1,i,phi,u_h1,j,5);
	      //+ alpha*phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      
	      
	      
	      B[i][j] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+N1][j+N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+2*N1][j+2*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+3*N1][j+3*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+4*N1][j+4*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+5*N1][j+5*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+6*N1][j+6*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	      B[i+7*N1][j+7*N1] = phgQuadBasDotBas(e,u_h1,i,u_h1,j,5);
	    }
	}
//check;	
// loop on basis functions
	for(i=0;i<N1;i++)
	  {
	    I[i] = phgMapE2L(matA->cmap,0,e,i);
//check;
	    I[i+N1] = phgMapE2L(matA->cmap,1,e,i);

//check;
	    I[i+2*N1] = phgMapE2L(matA->cmap,2,e,i);
//check;
	    I[i+3*N1] = phgMapE2L(matA->cmap,3,e,i);
//check;
	    I[i+4*N1] = phgMapE2L(matA->cmap,4,e,i);
	    I[i+5*N1] = phgMapE2L(matA->cmap,5,e,i);
	    I[i+6*N1] = phgMapE2L(matA->cmap,6,e,i);
	    I[i+7*N1] = phgMapE2L(matA->cmap,7,e,i);
	  }
//check;	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i,N,I,A[i]);
	  phgMatAddEntries(matB,1,I+i,N,I,B[i]);
	}
	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+N1,N,I,A[i+N1]);
	  phgMatAddEntries(matB,1,I+i+N1,N,I,B[i+N1]);
	}
	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+2*N1,N,I,A[i+2*N1]);
	  phgMatAddEntries(matB,1,I+i+2*N1,N,I,B[i+2*N1]);
	}

	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+3*N1,N,I,A[i+3*N1]);
	  phgMatAddEntries(matB,1,I+i+3*N1,N,I,B[i+3*N1]);
	}
	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+4*N1,N,I,A[i+4*N1]);
	  phgMatAddEntries(matB,1,I+i+4*N1,N,I,B[i+4*N1]);
	}
	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+5*N1,N,I,A[i+5*N1]);
	  phgMatAddEntries(matB,1,I+i+5*N1,N,I,B[i+5*N1]);
	}
	
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
	  phgMatAddEntries(matA,1,I+i+6*N1,N,I,A[i+6*N1]);
	  phgMatAddEntries(matB,1,I+i+6*N1,N,I,B[i+6*N1]);
	}
	  
	for(i=0;i<N1;i++)
	{
//	  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
//	  continue;
  	  phgMatAddEntries(matA,1,I+i+7*N1,N,I,A[i+7*N1]);
	  phgMatAddEntries(matB,1,I+i+7*N1,N,I,B[i+7*N1]);
	}

	//for(i=0;i<N1;i++)
	//{
	//  if(phgDofDirichletBC(u_h1,e,i,NULL,NULL,NULL,DOF_PROJ_NONE))
	//  continue;
	//  phgMatAddEntries(matA,1,I+i,N,I,A[i]);
	//  phgMatAddEntries(matA,1,I+i+N1,N,I,A[i+N1]);
	//  phgMatAddEntries(matA,1,I+i+2*N1,N,I,A[i+2*N1]);
	//  phgMatAddEntries(matA,1,I+i+3*N1,N,I,A[i+3*N1]);
	//  phgMatAddEntries(matA,1,I+i+4*N1,N,I,A[i+4*N1]);
	//  phgMatAddEntries(matA,1,I+i+5*N1,N,I,A[i+5*N1]);
	//  phgMatAddEntries(matA,1,I+i+6*N1,N,I,A[i+6*N1]);
	//  phgMatAddEntries(matA,1,I+i+7*N1,N,I,A[i+7*N1]);
	//  
	//  phgMatAddEntries(matB,1,I+i,N,I,B[i]);
	//  phgMatAddEntries(matB,1,I+i+N1,N,I,B[i+N1]);
	//  phgMatAddEntries(matB,1,I+i+2*N1,N,I,B[i+2*N1]);
	//  phgMatAddEntries(matB,1,I+i+3*N1,N,I,B[i+3*N1]);
	//  phgMatAddEntries(matB,1,I+i+4*N1,N,I,B[i+4*N1]);
	//  phgMatAddEntries(matB,1,I+i+5*N1,N,I,B[i+5*N1]);
	//  phgMatAddEntries(matB,1,I+i+6*N1,N,I,B[i+6*N1]);
	//  phgMatAddEntries(matB,1,I+i+7*N1,N,I,B[i+7*N1]);
        //}
//check;
  }
check;
}


static FLOAT
estimate_error(FLOAT lambda, DOF *u_h, DOF *error)
/* compute H1 error indicator */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *grad_u, *jump, *residual, *tmp;
    FLOAT PEoo = 0.0;

    grad_u = phgDofGradient(u_h, NULL, NULL, NULL);
    tmp = phgDofDivergence(grad_u, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(u_h, &residual, NULL, NULL);
    phgDofAXPBY(1.0, tmp, lambda, &residual);
    phgDofFree(&tmp);
    jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	eta = 0.0;
	/* for each face F compute [grad_u \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    eta += *DofFaceData(jump, e->faces[i]) * h;
	}
	eta = eta * .5 + diam*diam*phgQuadDofDotDof(e, residual, residual, -1);
	if (*DofElementData(error, e->index) < eta)
	    *DofElementData(error, e->index) = eta;
	if (PEoo < eta)
	    PEoo = eta;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);
    phgDofFree(&grad_u);

    return Sqrt(PEoo);
}

static int
bc_map(int bctype)
{
    return DIRICHLET;	/* set Dirichlet BC on all boundaries */
}
//
/* Note: make USER_CFLAGS="-DMATRIX_FREE_TEST=1" eigen */
#ifndef MATRIX_FREE_TEST
# define MATRIX_FREE_TEST 0
#endif	/* !defined(MATRIX_FREE_TEST) */

#if MATRIX_FREE_TEST
/* callback function for matrix A */
static int
mat_vec(MAT_OP op, MAT *mat, VEC *x, VEC *y)
/* computes y = op(A) * x */
{
    phgMatVec(op, 1.0, ((MAT **)mat->mv_data)[0], x, 0.0, &y);
    return 0;
}
#endif	/* MATRIX_FREE_TEST */


///////////////////////////////////////////////////////////////////////////////////////////////main

int
main(int argc, char *argv[])
{
    //static char *fn = "./cell.mesh";
    static char *fn = "./cube4.dat";
    /*static char *fn = "../test/model.mesh";*/
   // INT periodicity = X_MASK|Y_MASK|Z_MASK ;
    static INT mem_max = 400000;
    size_t mem, mem_peak;
    int i, j, k, n, nit;
    INT nev = 12;
    INT pre_refines = 8;
    GRID *g;
    DOF **u_h1,**u_h2,**u_h3,**u_h4,**u_h5,**u_h6,**u_h7,**u_h8,*error;
    DOF *phi, *A1, *A2, *A3;
    MAP *map;
    MAT *A, *B;
    
    FLOAT tol = 1.e-3, tau = 1.0, PEoo, thres;
    phgPrintf("\n mem_max = %d,nev = %d, tau = %0.9lf \n",mem_max,nev,tau);
    FLOAT *evals;
    double wtime;
#if MATRIX_FREE_TEST
    MAT *C;

    phgOptionsPreset("-arpack_solver pcg");
#endif	/* MATRIX_FREE_TEST */

    phgOptionsPreset("-dof_type P1");

   // phgOptionsRegisterInt("periodicity", "Set periodicity", &periodicity);
    phgOptionsRegisterFilename("-mesh_file", "Mesh filename", &fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-nev", "Number of eigenvalues", &nev);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);
    phgOptionsRegisterFloat("-tau", "The shift", &tau);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);

    phgInit(&argc, &argv);
    //phgPrintf("\n mem_max = %d,nev = %d, tau = %0.2lf \n",mem_max,nev,tau);

    if (DOF_DEFAULT->mass_lumping == NULL)
	phgPrintf("Order of FE bases: %d\n", DOF_DEFAULT->order);
    else
	phgPrintf("Order of FE bases: %d\n", DOF_DEFAULT->mass_lumping->order0);
//check;
    g = phgNewGrid(-1);
    //phgSetPeriodicity(g,periodicity);
    //phgImportSetBdryMapFunc(bc_map);
   phgSetPeriodicity(g, X_MASK | Y_MASK | Z_MASK);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    /* pre-refinement */
    phgRefineAllElements(g, pre_refines);
//check;
    phi = phgDofNew(g, DOF_DEFAULT, 1, "phi", func_phi);
    A1 = phgDofNew(g, DOF_DEFAULT, 1, "A1", func_A1);
    A2 = phgDofNew(g, DOF_DEFAULT, 1, "A2", func_A2);
    A3 = phgDofNew(g, DOF_DEFAULT, 1, "A3", func_A3);

    //u_h = phgAlloc(nev * sizeof(*u_h));
    u_h1 = phgAlloc(nev * sizeof(*u_h1));
    u_h2 = phgAlloc(nev * sizeof(*u_h2));
    u_h3 = phgAlloc(nev * sizeof(*u_h3));
    u_h4 = phgAlloc(nev * sizeof(*u_h4));
    u_h5 = phgAlloc(nev * sizeof(*u_h5));
    u_h6 = phgAlloc(nev * sizeof(*u_h6));
    u_h7 = phgAlloc(nev * sizeof(*u_h7));
    u_h8 = phgAlloc(nev * sizeof(*u_h8));
    
    evals = phgCalloc(nev, sizeof(*evals));
    for (i = 0; i < nev; i++) {
	//u_h[i] = phgDofNew(g, DOF_DEFAULT, 8, "u_h", DofInterpolation);
	u_h1[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h1", DofInterpolation);
	u_h2[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h2", DofInterpolation);
	u_h3[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h3", DofInterpolation);
	u_h4[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h4", DofInterpolation);
	u_h5[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h5", DofInterpolation);
	u_h6[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h6", DofInterpolation);
	u_h7[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h7", DofInterpolation);
	u_h8[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h8", DofInterpolation);
#if 0
	/* All-Neumann BC */
	phgDofSetDirichletBoundaryMask(u_h[i], 0);
#else
	/* All-Dirichlet BC */
	//phgDofSetDirichletBoundaryMask(u_h[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h1[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h2[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h3[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h4[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h5[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h6[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h7[i], BDRY_MASK);
//	phgDofSetDirichletBoundaryMask(u_h8[i], BDRY_MASK);
#endif
	/* set random initial values for the eigenvectors */
	//phgDofRandomize(u_h[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h1[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h2[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h3[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h4[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h5[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h6[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h7[i], i == 0 ? 123 : 0);
	phgDofRandomize(u_h8[i], i == 0 ? 123 : 0);
	
    }
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

check;

    nev_plus = nev + 1000;
    modes = phgAlloc(nev_plus * sizeof(modes[0]));
    for (n = 0, nit = 3; n < nev_plus; nit++) {
	/* find i, j, k such that i^2 + j^2 + k^2 == nit */
	for (i = 1; n < nev_plus && i * i <= nit; i++) {
	    for (j = 1; n < nev_plus && j * j <= nit - i * i; j++) {
		for (k = 1; k * k < nit - i * i - j * j; k++);
		if (i * i + j * j + k * k != nit)
		    continue;
		modes[n][0] = i;
		modes[n][1] = j;
		modes[n][2] = k;
		n++;
	    }
	}
    }

check;

#if 0 * USE_BLOPEX
const char *pc1 = phgOptionsGetKeyword("-blopex_pc_solver1");
/* disable local PC for the first loop */
phgOptionsSetKeyword("-blopex_pc_solver1", "none");
#else
const char *pc1 = NULL;
#endif


	FLOAT kx = 0.0, ky = 0.0, kz = 0.0;   
	INT kxk = 0, kxN = 40;
	FLOAT eigs0[kxN], eigs1[kxN];
	for(kxk = 0; kxk < kxN; kxk++)
	{
	  kx = kxk*2*M_PI/kxN; 


    //while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh\n");
	phgPrintf("%"dFMT" DOF, %"dFMT
		  " elements, %d submesh%s, load imbalance: %lg\n",
			DofGetDataCountGlobal(u_h1[0])*8, g->nleaf_global,
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	wtime = phgGetTime(NULL);
	map = phgMapCreate(u_h1[0],u_h2[0],u_h3[0],u_h4[0],u_h5[0],u_h6[0],u_h7[0],u_h8[0], NULL);
check;
#if 1 
	A = phgMapCreateMat(map, map);
	B = phgMapCreateMat(map, map);
	build_matrices(A, B,  u_h1[0], phi, A1, A2, A3, kx, ky, kz);
	phgMatRemoveBoundaryEntries(A);
	phgMatRemoveBoundaryEntries(B);
#else
	{
	    MAP *m = phgMapRemoveBoundaryEntries(map);
	    A = phgMapCreateMat(m, m);
	    B = phgMapCreateMat(m, m);
	    phgMapDestroy(&m);
	}
check;	
	build_matrices(A, B, u_h1[0], phi, A1, A2, A3, kx, ky, kz);

#endif
	phgPrintf("Build matrices, matrix size: %"dFMT", wtime: %0.2lfs\n",
			A->rmap->nglobal, phgGetTime(NULL) - wtime);
	wtime = phgGetTime(NULL);

#if MATRIX_FREE_TEST
	C = phgMapCreateMatrixFreeMat(A->rmap, A->cmap, mat_vec, A, NULL);
	n = phgDofEigenSolve(C, B, nev, EIGEN_CLOSEST, tau, &nit, evals,
         		     map, u_h1,u_h2,u_h3,u_h4,u_h5,u_h6,u_h7,u_h8, NULL);
	phgMatDestroy(&C);
#else	/* MATRIX_FREE_TEST */
	n = phgDofEigenSolve(A, B, nev, EIGEN_SMALLEST, tau, &nit, evals,
			     map, u_h1,u_h2,u_h3,u_h4,u_h5,u_h6,u_h7,u_h8, NULL);
#endif	/* MATRIX_FREE_TEST */
	phgPrintf("%d iterations, converged eigenvalues: %d, wtime: %0.2lfs\n",
			nit, n, phgGetTime(NULL) - wtime);

if (pc1 != NULL && n == nev) phgOptionsSetKeyword("-blopex_pc_solver1", pc1);
	for (i = 0; i < n; i++) {
	    FLOAT err_tau;
	    nit = modes[i][0] * modes[i][0] +
		  modes[i][1] * modes[i][1] +
		  modes[i][2] * modes[i][2];
	    err_tau = Fabs(evals[i]*evals[i] - nit * Pow(M_PI,2)) / (evals[i]*evals[i]);
	    phgPrintf("  tau[%d]=%0.12e, error(tau)=%0.2e, error(u_h1)=%0.2e\n",
			i, (double)evals[i], (double)err_tau,
			(double)compute_error(u_h1[i], i));
	}
	*(eigs0 + kxk) = evals[0];
	*(eigs1 + kxk) = evals[4];
	
	compute_error(NULL, 0);		/* reset compute_error() */
	for ( ; i < nev; i++){
            phgDofRandomize(u_h1[i], 0);
            phgDofRandomize(u_h2[i], 0);
            phgDofRandomize(u_h3[i], 0);
            phgDofRandomize(u_h4[i], 0);
            phgDofRandomize(u_h5[i], 0);
            phgDofRandomize(u_h6[i], 0);
            phgDofRandomize(u_h7[i], 0);
            phgDofRandomize(u_h8[i], 0);
	}


	if (n == 0) {
	    PEoo = 1e10;
	}
	else {
	    PEoo = 0.0;
	    phgDofSetDataByValue(error, 0.0);
	    for (i = 0; i < n; i++) {
		thres = estimate_error(evals[i], u_h1[i], error);
		if (PEoo < thres)
		    PEoo = thres;
	    }
	}
	phgPrintf("indicator=%0.3le\n", (double)PEoo);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)mem / (1024.0 * 1024.0),
		(double)mem_peak / (1024.0 * 1024.0));
	if (PEoo < tol || mem_peak >= 1024 * (size_t)1024 * mem_max)
	    break;
	//if (n == 0) {
	    //phgRefineAllElements(g, 1);
	//}
	//else {
#if 0
	    /*thres = Pow(1e-2,2) / (double)g->nleaf_global;*/
	    thres = Pow(PEoo * 0.5, 2);
	    ForAllElements(g, e)
		if (*DofElementData(error, e->index) > thres)
		    e->mark = 1;
#else
	//    phgMarkRefine(MARK_DEFAULT, error, Pow(0.8,2), NULL, 0., 1,
	//		  Pow(tol, 2) / g->nleaf_global);
#endif
    //        phgRefineAllElements(g, 1);
	    //phgRefineMarkedElements(g);
	//}
	
    }
        phgMatDestroy(&B);
	phgMatDestroy(&A);
	phgMapDestroy(&map);
    char *sn[2];
    FILE *fp[2];
    sn[0]=(char *)phgAlloc(30*sizeof(char));
    sn[1]=(char *)phgAlloc(30*sizeof(char));
    *sn[0]='\0';
    *sn[1]='\0';
    phgPrintf("\n output eigs00 and eigs01\n");
    sprintf(sn[0],"vtk1/eigs00");
    sprintf(sn[1],"vtk1/eigs01");
    fp[0]=fopen(sn[0],"w");
    fp[1]=fopen(sn[1],"w");
    for(i=0;i<kxN;i++)
    {
      fprintf(fp[0],"%10.9f\n",(double)eigs0[i]);
      fprintf(fp[1],"%10.9f\n",(double)eigs1[i]);
    }
    fclose(fp[0]);
    fclose(fp[1]);
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector1.vtk", nev, u_h1));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector2.vtk", nev, u_h2));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector3.vtk", nev, u_h3));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector4.vtk", nev, u_h4));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector5.vtk", nev, u_h5));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector6.vtk", nev, u_h6));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector7.vtk", nev, u_h7));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk1/D_eigenvector8.vtk", nev, u_h8));

#if 0
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector1.vtk", nev, u_h1));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector2.vtk", nev, u_h2));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector3.vtk", nev, u_h3));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector4.vtk", nev, u_h4));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector5.vtk", nev, u_h5));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector6.vtk", nev, u_h6));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector7.vtk", nev, u_h7));
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "vtk/D_eigenvector8.vtk", nev, u_h8));
#endif

    phgFree(modes);
    phgFree(evals);
    phgFree(eigs0);
    phgFree(eigs1);
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h1 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h2 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h3 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h4 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h5 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h6 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h7 + i);
    }
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h8 + i);
    }
    phgDofFree(&error);
    phgFree(u_h1);
    phgFree(u_h2);
    phgFree(u_h3);
    phgFree(u_h4);
    phgFree(u_h5);
    phgFree(u_h6);
    phgFree(u_h7);
    phgFree(u_h8);
    phgFree(phi);
    phgFree(A1);
    phgFree(A2);
    phgFree(A3);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}




