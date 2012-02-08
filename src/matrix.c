/*** Translated to the C language by N. Kyriazis  20 Aug 2003 ***

 Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,
 tape15,tape16,tape20,tape21)

 Numerical Electromagnetics Code (NEC2)  developed at Lawrence
 Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414
 for problems with the NEC code. For problems with the vax implem-
 entation, contact J. Breakall at 415-422-8196 or E. Domning at 415
 422-5936)
 file created 4/11/80.

                ***********Notice**********
 This computer code material was prepared as an account of work
 sponsored by the United States government.  Neither the United
 States nor the United States Department Of Energy, nor any of
 their employees, nor any of their contractors, subcontractors,
 or their employees, makes any warranty, express or implied, or
 assumes any legal liability or responsibility for the accuracy,
 completeness or usefulness of any information, apparatus, product
 or process disclosed, or represents that its use would not infringe
 privately-owned rights.

*******************************************************************/

#include "nec2c.h"
#include <omp.h>

/* common  /data/ */
extern data_t data;

/* common  /dataj/ */
extern dataj_t dataj;

/* common  /matpar/ */
extern matpar_t matpar;

/* common  /segj/ */
extern segj_t segj;

/* common  /zload/ */
extern zload_t zload;

/* common  /smat/ */
extern smat_t smat;

/* common  /gnd/ */
extern gnd_t gnd;

/* common  /vsorc/ */
extern vsorc_t vsorc;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/*-------------------------------------------------------------------*/

/* cmset sets up the complex structure matrix in the array cm */
void cmset( int nrow, complex long double *cm, long double rkhx, int iexkx )
{
  int mp2, neq, npeq, it, i, j, i1, i2, in2;
  //int iout;
  int im1, im2, ist, ij, ipr, jss, jm1, jm2, jst, k, ka, kk;
  complex long double zaj, deter, *scm = NULL;

  mp2=2* data.mp;
  npeq= data.np+ mp2;
  neq= data.n+2* data.m;
  smat.nop = neq/npeq;

  dataj.rkh= rkhx;
  dataj.iexk= iexkx;
  //iout=2* matpar.npblk* nrow;
  it= matpar.nlast;

  for( i = 0; i < nrow; i++ )
	for( j = 0; j < it; j++ )
	  cm[i+j*nrow]= CPLX_00;

  i1= 1;
  i2= it;
  in2= i2;

  if( in2 > data.np)
	in2= data.np;

  im1= i1- data.np;
  im2= i2- data.np;

  if( im1 < 1)
	im1=1;

  ist=1;
  if( i1 <= data.np)
	ist= data.np- i1+2;

  /* wire source loop */
  if( data.n != 0)
  {
	for( j = 1; j <= data.n; j++ )
	{
	  trio(j);
	  for( i = 0; i < segj.jsno; i++ )
	  {
		ij= segj.jco[i];
		segj.jco[i]=(( ij-1)/ data.np)* mp2+ ij;
	  }

	  if( i1 <= in2)
		cmww( j, i1, in2, cm, nrow, cm, nrow,1);

	  if( im1 <= im2)
		cmws( j, im1, im2, &cm[(ist-1)*nrow], nrow, cm, nrow, 1);

	  /* matrix elements modified by loading */
	  if( zload.nload == 0)
		continue;

	  if( j > data.np)
		continue;

	  ipr= j;
	  if( (ipr < 1) || (ipr > it) )
		continue;

	  zaj= zload.zarray[j-1];

	  for( i = 0; i < segj.jsno; i++ )
	  {
		jss= segj.jco[i];
		cm[(jss-1)+(ipr-1)*nrow] -= ( segj.ax[i]+ segj.cx[i])* zaj;
	  }

	} /* for( j = 1; j <= n; j++ ) */

  } /* if( n != 0) */

  if( data.m != 0)
  {
	/* matrix elements for patch current sources */
	jm1=1- data.mp;
	jm2=0;
	jst=1- mp2;

	for( i = 0; i < smat.nop; i++ )
	{
	  jm1 += data.mp;
	  jm2 += data.mp;
	  jst += npeq;

	  if( i1 <= in2)
		cmsw( jm1, jm2, i1, in2, &cm[(jst-1)], cm, 0, nrow, 1);

	  if( im1 <= im2)
		cmss( jm1, jm2, im1, im2, &cm[(jst-1)+(ist-1)*nrow], nrow, 1);
	}

  } /* if( m != 0) */

  if( matpar.icase == 1)
	return;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, data.np2m * sizeof(complex long double) );

  /* combine elements for symmetry modes */
  for( i = 0; i < it; i++ )
  {
	for( j = 0; j < npeq; j++ )
	{
	  for( k = 0; k < smat.nop; k++ )
	  {
		ka= j+ k*npeq;
		scm[k]= cm[ka+i*nrow];
	  }

	  deter= scm[0];

	  for( kk = 1; kk < smat.nop; kk++ )
		deter += scm[kk];

	  cm[j+i*nrow]= deter;

	  for( k = 1; k < smat.nop; k++ )
	  {
		ka= j+ k*npeq;
		deter= scm[0];

		for( kk = 1; kk < smat.nop; kk++ )
		{
		  deter += scm[kk]* smat.ssx[k+kk*smat.nop];
		  cm[ka+i*nrow]= deter;
		}

	  } /* for( k = 1; k < smat.nop; k++ ) */

	} /* for( j = 0; j < npeq; j++ ) */

  } /* for( i = 0; i < it; i++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* cmss computes matrix elements for surface-surface interactions. */
void cmss( int j1, int j2, int im1, int im2,
	complex long double *cm, int nrow, int itrp )
{
  int i1, i2, icomp, ii1, i, il, ii2, jj1, j, jl, jj2;
  long double t1xi, t1yi, t1zi, t2xi, t2yi, t2zi, xi, yi, zi;
  complex long double g11, g12, g21, g22;

  i1=( im1+1)/2;
  i2=( im2+1)/2;
  icomp= i1*2-3;
  ii1=-2;
  if( icomp+2 < im1)
	ii1=-3;

  /* loop over observation patches */
  il = -1;
  for( i = i1; i <= i2; i++ )
  {
	il++;
	icomp += 2;
	ii1 += 2;
	ii2 = ii1+1;

	t1xi= data.t1x[il]* data.psalp[il];
	t1yi= data.t1y[il]* data.psalp[il];
	t1zi= data.t1z[il]* data.psalp[il];
	t2xi= data.t2x[il]* data.psalp[il];
	t2yi= data.t2y[il]* data.psalp[il];
	t2zi= data.t2z[il]* data.psalp[il];
	xi= data.px[il];
	yi= data.py[il];
	zi= data.pz[il];

	/* loop over source patches */
	jj1=-2;
	for( j = j1; j <= j2; j++ )
	{
	  jl=j-1;
	  jj1 += 2;
	  jj2 = jj1+1;

	  dataj.s= data.pbi[jl];
	  dataj.xj= data.px[jl];
	  dataj.yj= data.py[jl];
	  dataj.zj= data.pz[jl];
	  dataj.t1xj= data.t1x[jl];
	  dataj.t1yj= data.t1y[jl];
	  dataj.t1zj= data.t1z[jl];
	  dataj.t2xj= data.t2x[jl];
	  dataj.t2yj= data.t2y[jl];
	  dataj.t2zj= data.t2z[jl];

	  hintg( xi, yi, zi);

	  g11=-( t2xi* dataj.exk+ t2yi* dataj.eyk+ t2zi* dataj.ezk);
	  g12=-( t2xi* dataj.exs+ t2yi* dataj.eys+ t2zi* dataj.ezs);
	  g21=-( t1xi* dataj.exk+ t1yi* dataj.eyk+ t1zi* dataj.ezk);
	  g22=-( t1xi* dataj.exs+ t1yi* dataj.eys+ t1zi* dataj.ezs);

	  if( i == j )
	  {
		g11 -= .5;
		g22 += .5;
	  }

	  /* normal fill */
	  if( itrp == 0)
	  {
		if( icomp >= im1 )
		{
		  cm[ii1+jj1*nrow]= g11;
		  cm[ii1+jj2*nrow]= g12;
		}

		if( icomp >= im2 )
		  continue;

		cm[ii2+jj1*nrow]= g21;
		cm[ii2+jj2*nrow]= g22;
		continue;

	  } /* if( itrp == 0) */

	  /* transposed fill */
	  if( icomp >= im1 )
	  {
		cm[jj1+ii1*nrow]= g11;
		cm[jj2+ii1*nrow]= g12;
	  }

	  if( icomp >= im2 )
		continue;

	  cm[jj1+ii2*nrow]= g21;
	  cm[jj2+ii2*nrow]= g22;

	} /* for( j = j1; j <= j2; j++ ) */

  } /* for( i = i1; i <= i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* computes matrix elements for e along wires due to patch current */
void cmsw( int j1, int j2, int i1, int i2, complex long double *cm,
	complex long double *cw, int ncw, int nrow, int itrp )
{
  //int neqs; 
  int k, icgo, i, ipch, jl, j, js, il, ip;
  int jsnox; /* -1 offset to "jsno" for array indexing */
  long double xi, yi, zi, cabi, sabi, salpi, fsign=1., pyl, pxl;
  complex long double emel[9];

  //neqs= data.np2m;
  jsnox = segj.jsno-1;

  if( itrp >= 0)
  {
	k=-1;
	icgo=0;

	/* observation loop */
	for( i = i1-1; i < i2; i++ )
	{
	  k++;
	  xi= data.x[i];
	  yi= data.y[i];
	  zi= data.z[i];
	  cabi= data.cab[i];
	  sabi= data.sab[i];
	  salpi= data.salp[i];
	  ipch=0;

	  if( data.icon1[i] >= PCHCON)
	  {
		ipch= data.icon1[i]-PCHCON;
		fsign=-1.;
	  }

	  if( data.icon2[i] >= PCHCON)
	  {
		ipch= data.icon2[i]-PCHCON;
		fsign=1.;
	  }

	  /* source loop */
	  jl = -1;
	  for( j = j1; j <= j2; j++ )
	  {
		jl += 2;
		js = j-1;
		dataj.t1xj= data.t1x[js];
		dataj.t1yj= data.t1y[js];
		dataj.t1zj= data.t1z[js];
		dataj.t2xj= data.t2x[js];
		dataj.t2yj= data.t2y[js];
		dataj.t2zj= data.t2z[js];
		dataj.xj= data.px[js];
		dataj.yj= data.py[js];
		dataj.zj= data.pz[js];
		dataj.s= data.pbi[js];

		/* ground loop */
		for( ip = 1; ip <= gnd.ksymp; ip++ )
		{
		  dataj.ipgnd= ip;

		  if( ((ipch == j) || (icgo != 0)) && (ip != 2) )
		  {
			if( icgo <= 0 )
			{
			  pcint( xi, yi, zi, cabi, sabi, salpi, emel);

			  pyl= PI* data.si[i]* fsign;
			  pxl= sinl( pyl);
			  pyl= cosl( pyl);
			  dataj.exc= emel[8]* fsign;

			  trio(i+1);

			  il= i-ncw;
			  if( i < data.np)
				il += (il/data.np)*2*data.mp;

			  if( itrp == 0 )
				cw[k+il*nrow] +=
				  dataj.exc*( segj.ax[jsnox]+ segj.bx[jsnox]* pxl+ segj.cx[jsnox]* pyl);
			  else
				cw[il+k*nrow] +=
				  dataj.exc*( segj.ax[jsnox]+ segj.bx[jsnox]* pxl+ segj.cx[jsnox]* pyl);

			} /* if( icgo <= 0 ) */

			if( itrp == 0)
			{
			  cm[k+(jl-1)*nrow]= emel[icgo];
			  cm[k+jl*nrow]    = emel[icgo+4];
			}
			else
			{
			  cm[(jl-1)+k*nrow]= emel[icgo];
			  cm[jl+k*nrow]    = emel[icgo+4];
			}

			icgo++;
			if( icgo == 4)
			  icgo=0;

			continue;

		  } /* if( ((ipch == (j+1)) || (icgo != 0)) && (ip != 2) ) */

		  unere( xi, yi, zi);

		  /* normal fill */
		  if( itrp == 0)
		  {
			cm[k+(jl-1)*nrow] += dataj.exk* cabi+ dataj.eyk* sabi+ dataj.ezk* salpi;
			cm[k+jl*nrow]     += dataj.exs* cabi+ dataj.eys* sabi+ dataj.ezs* salpi;
			continue;
		  }

		  /* transposed fill */
		  cm[(jl-1)+k*nrow] += dataj.exk* cabi+ dataj.eyk* sabi+ dataj.ezk* salpi;
		  cm[jl+k*nrow]     += dataj.exs* cabi+ dataj.eys* sabi+ dataj.ezs* salpi;

		} /* for( ip = 1; ip <= gnd.ksymp; ip++ ) */

	  } /* for( j = j1; j <= j2; j++ ) */

	} /* for( i = i1-1; i < i2; i++ ) */

  } /* if( itrp >= 0) */

  return;
}

/*-----------------------------------------------------------------------*/

/* cmws computes matrix elements for wire-surface interactions */
void cmws( int j, int i1, int i2, complex long double *cm,
	int nr, complex long double *cw, int nw, int itrp )
{
  int ipr, i, ipatch, ik, js=0, ij, jx;
  long double xi, yi, zi, tx, ty, tz;
  complex long double etk, ets, etc;
  i=nw;
  j--;
  dataj.s= data.si[j];
  dataj.b= data.bi[j];
  dataj.xj= data.x[j];
  dataj.yj= data.y[j];
  dataj.zj= data.z[j];
  dataj.cabj= data.cab[j];
  dataj.sabj= data.sab[j];
  dataj.salpj= data.salp[j];

  /* observation loop */
  ipr= -1;
  for( i = i1; i <= i2; i++ )
  {
	ipr++;
	ipatch=(i+1)/2;
	ik= i-( i/2)*2;

	if( (ik != 0) || (ipr == 0) )
	{
	  js= ipatch-1;
	  xi= data.px[js];
	  yi= data.py[js];
	  zi= data.pz[js];
	  hsfld( xi, yi, zi,0.);

	  if( ik != 0 )
	  {
		tx= data.t2x[js];
		ty= data.t2y[js];
		tz= data.t2z[js];
	  }
	  else
	  {
		tx= data.t1x[js];
		ty= data.t1y[js];
		tz= data.t1z[js];
	  }

	} /* if( (ik != 0) || (ipr == 0) ) */
	else
	{
	  tx= data.t1x[js];
	  ty= data.t1y[js];
	  tz= data.t1z[js];

	} /* if( (ik != 0) || (ipr == 0) ) */

	etk=-( dataj.exk* tx+ dataj.eyk* ty+ dataj.ezk* tz)* data.psalp[js];
	ets=-( dataj.exs* tx+ dataj.eys* ty+ dataj.ezs* tz)* data.psalp[js];
	etc=-( dataj.exc* tx+ dataj.eyc* ty+ dataj.ezc* tz)* data.psalp[js];

	/* fill matrix elements.  element locations */
	/* determined by connection data. */

	/* normal fill */
	if( itrp == 0)
	{
	  for( ij = 0; ij < segj.jsno; ij++ )
	  {
		jx= segj.jco[ij]-1;
		cm[ipr+jx*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }

	  continue;
	} /* if( itrp == 0) */

	/* transposed fill */
	if( itrp != 2)
	{
	  for( ij = 0; ij < segj.jsno; ij++ )
	  {
		jx= segj.jco[ij]-1;
		cm[jx+ipr*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }

	  continue;
	} /* if( itrp != 2) */

	/* transposed fill - c(ws) and d(ws)prime (=cw) */
	for( ij = 0; ij < segj.jsno; ij++ )
	{
	  jx= segj.jco[ij]-1;
	  if( jx < nr)
		cm[jx+ipr*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  else
	  {
		jx -= nr;
		cw[jx+ipr*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }
	} /* for( ij = 0; ij < segj.jsno; ij++ ) */

  } /* for( i = i1; i <= i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* cmww computes matrix elements for wire-wire interactions */
void cmww( int j, int i1, int i2, complex long double *cm,
	int nr, complex long double *cw, int nw, int itrp)
{
  int ipr, iprx, i, ij, jx;
  long double xi, yi, zi, ai, cabi, sabi, salpi;
  complex long double etk, ets, etc;

  /* set source segment parameters */
  jx = j;
  j--;
  dataj.s= data.si[j];
  dataj.b= data.bi[j];
  dataj.xj= data.x[j];
  dataj.yj= data.y[j];
  dataj.zj= data.z[j];
  dataj.cabj= data.cab[j];
  dataj.sabj= data.sab[j];
  dataj.salpj= data.salp[j];

  /* decide whether ext. t.w. approx. can be used */
  if( dataj.iexk != 0)
  {
	ipr = data.icon1[j];
	if (ipr > PCHCON) dataj.ind1 = 0;
	else if( ipr < 0 )
	{
	  ipr= -ipr;
	  iprx= ipr-1;

	  if( -data.icon1[iprx] != jx )	dataj.ind1 = 2;
	  else
	  {
		xi= fabsl( dataj.cabj* data.cab[iprx]+ dataj.sabj*
			data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
		if( (xi < 0.999999) || (fabsl(data.bi[iprx]/dataj.b-1.) > 1.e-6) )
		  dataj.ind1=2;
		else
		  dataj.ind1=0;

	  } /* if( -data.icon1[iprx] != jx ) */

	} /* if( ipr < 0 ) */
	else
	{
	  iprx = ipr-1;
	  if( ipr == 0 ) dataj.ind1=1;
	  else
	  {
		if( ipr != jx )
		{
		  if( data.icon2[iprx] != jx ) dataj.ind1=2;
		  else
		  {
			xi= fabsl( dataj.cabj* data.cab[iprx]+ dataj.sabj*
				data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
			if( (xi < 0.999999) || (fabsl(data.bi[iprx]/dataj.b-1.) > 1.e-6) )
			  dataj.ind1=2;
			else
			  dataj.ind1=0;

		  } /* if( data.icon2[iprx] != jx ) */

		} /* if( ipr != jx ) */
		else
		  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.e-8)
			dataj.ind1=2;
		  else
			dataj.ind1=0;

	  } /* if( ipr == 0 ) */

	} /* if( ipr < 0 ) */

	ipr = data.icon2[j];
	if (ipr > PCHCON) dataj.ind2 = 2;
	else if( ipr < 0 )
	{
	  ipr= -ipr;
	  iprx = ipr-1;
	  if( -data.icon2[iprx] != jx )
		dataj.ind2=2;
	  else
	  {
		xi= fabsl( dataj.cabj* data.cab[iprx]+ dataj.sabj*
			data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
		if( (xi < 0.999999) || (fabsl(data.bi[iprx]/dataj.b-1.) > 1.e-6) )
		  dataj.ind2=2;
		else
		  dataj.ind2=0;

	  } /* if( -data.icon1[iprx] != jx ) */

	} /* if( ipr < 0 ) */
	else
	{
	  iprx = ipr-1;
	  if( ipr == 0 )
		dataj.ind2=1;
	  else
	  {
		if( ipr != jx )
		{
		  if( data.icon1[iprx] != jx )
			dataj.ind2=2;
		  else
		  {
			xi= fabsl( dataj.cabj* data.cab[iprx]+ dataj.sabj*
				data.sab[iprx]+ dataj.salpj* data.salp[iprx]);
			if( (xi < 0.999999) || (fabsl(data.bi[iprx]/dataj.b-1.) > 1.e-6) )
			  dataj.ind2=2;
			else
			  dataj.ind2=0;

		  } /* if( data.icon2[iprx] != jx ) */

		} /* if( ipr != jx ) */
		else
		  if( dataj.cabj* dataj.cabj+ dataj.sabj* dataj.sabj > 1.e-8)
			dataj.ind2=2;
		  else
			dataj.ind2=0;

	  } /* if( ipr == 0 ) */

	} /* if( ipr < 0 ) */

  } /* if( dataj.iexk != 0) */

  /* observation loop */
  ipr=-1;
  for( i = i1-1; i < i2; i++ )
  {
	ipr++;
	ij= i-j;
	xi= data.x[i];
	yi= data.y[i];
	zi= data.z[i];
	ai= data.bi[i];
	cabi= data.cab[i];
	sabi= data.sab[i];
	salpi= data.salp[i];

	efld( xi, yi, zi, ai, ij);

	etk= dataj.exk* cabi+ dataj.eyk* sabi+ dataj.ezk* salpi;
	ets= dataj.exs* cabi+ dataj.eys* sabi+ dataj.ezs* salpi;
	etc= dataj.exc* cabi+ dataj.eyc* sabi+ dataj.ezc* salpi;

	/* fill matrix elements. element locations */
	/* determined by connection data. */

	/* normal fill */
	if( itrp == 0)
	{
	  for( ij = 0; ij < segj.jsno; ij++ )
	  {
		jx = segj.jco[ij]-1;
		cm[ipr+jx*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }
	  continue;
	}

	/* transposed fill */
	if( itrp != 2)
	{
	  for( ij = 0; ij < segj.jsno; ij++ )
	  {
		jx= segj.jco[ij]-1;
		cm[jx+ipr*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }
	  continue;
	}

	/* trans. fill for c(ww) - test for elements for d(ww)prime.  (=cw) */
	for( ij = 0; ij < segj.jsno; ij++ )
	{
	  jx= segj.jco[ij]-1;
	  if( jx < nr)
		cm[jx+ipr*nr] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  else
	  {
		jx -= nr;
		cw[jx*ipr*nw] += etk* segj.ax[ij]+ ets* segj.bx[ij]+ etc* segj.cx[ij];
	  }

	} /* for( ij = 0; ij < segj.jsno; ij++ ) */

  } /* for( i = i1-1; i < i2; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* etmns fills the array e with the negative of the */
/* electric field incident on the structure. e is the */
/* right hand side of the matrix equation. */
void etmns( long double p1, long double p2, long double p3, long double p4,
	long double p5, long double p6, int ipr, complex long double *e )
{
  int i, is, i1, i2=0, neq;
  long double cth, sth, cph, sph, cet, set, pxl, pyl, pzl, wx;
  long double wy, wz, qx, qy, qz, arg, ds, dsh, rs, r;
  complex long double cx, cy, cz, er, et, ezh, erh, rrv=CPLX_00, rrh=CPLX_00, tt1, tt2;

  neq= data.n+2*data.m;
  vsorc.nqds=0;

  /* applied field of voltage sources for transmitting case */
  if( (ipr == 0) || (ipr == 5) )
  {
	for( i = 0; i < neq; i++ )
	  e[i]=CPLX_00;

	if( vsorc.nsant != 0)
	{
	  for( i = 0; i < vsorc.nsant; i++ )
	  {
		is= vsorc.isant[i]-1;
		e[is]= -vsorc.vsant[i]/( data.si[is]* data.wlam);
	  }
	}

	if( vsorc.nvqd == 0)
	  return;

	for( i = 0; i < vsorc.nvqd; i++ )
	{
	  is= vsorc.ivqd[i];
	  qdsrc( is, vsorc.vqd[i], e);
	}
	return;

  } /* if( (ipr == 0) || (ipr == 5) ) */

  /* incident plane wave, linearly polarized. */
  if( ipr <= 3)
  {
	cth= cosl( p1);
	sth= sinl( p1);
	cph= cosl( p2);
	sph= sinl( p2);
	cet= cosl( p3);
	set= sinl( p3);
	pxl= cth* cph* cet- sph* set;
	pyl= cth* sph* cet+ cph* set;
	pzl=- sth* cet;
	wx=- sth* cph;
	wy=- sth* sph;
	wz=- cth;
	qx= wy* pzl- wz* pyl;
	qy= wz* pxl- wx* pzl;
	qz= wx* pyl- wy* pxl;

	if( gnd.ksymp != 1)
	{
	  if( gnd.iperf != 1)
	  {
		rrv= csqrtl(1.- gnd.zrati* gnd.zrati* sth* sth);
		rrh= gnd.zrati* cth;
		rrh=( rrh- rrv)/( rrh+ rrv);
		rrv= gnd.zrati* rrv;
		rrv=-( cth- rrv)/( cth+ rrv);
	  }
	  else
	  {
		rrv=-CPLX_10;
		rrh=-CPLX_10;
	  } /* if( gnd.iperf != 1) */

	} /* if( gnd.ksymp != 1) */

	if( ipr == 1)
	{
	  if( data.n != 0)
	  {
		for( i = 0; i < data.n; i++ )
		{
		  arg=- TP*( wx* data.x[i]+ wy* data.y[i]+ wz* data.z[i]);
		  e[i]=-( pxl* data.cab[i]+ pyl* data.sab[i]+ pzl*
			  data.salp[i])* cmplx( cosl( arg), sinl( arg));
		}

		if( gnd.ksymp != 1)
		{
		  tt1=( pyl* cph- pxl* sph)*( rrh- rrv);
		  cx= rrv* pxl- tt1* sph;
		  cy= rrv* pyl+ tt1* cph;
		  cz=- rrv* pzl;

		  for( i = 0; i < data.n; i++ )
		  {
			arg=- TP*( wx* data.x[i]+ wy* data.y[i]- wz* data.z[i]);
			e[i]= e[i]-( cx* data.cab[i]+ cy* data.sab[i]+
				cz* data.salp[i])* cmplx(cosl( arg), sinl( arg));
		  }

		} /* if( gnd.ksymp != 1) */

	  } /* if( data.n != 0) */

	  if( data.m == 0)
		return;

	  i= -1;
	  i1= data.n-2;
	  for( is = 0; is < data.m; is++ )
	  {
		i++;
		i1 += 2;
		i2 = i1+1;
		arg=- TP*( wx* data.px[i]+ wy* data.py[i]+ wz* data.pz[i]);
		tt1= cmplx( cosl( arg), sinl( arg))* data.psalp[i]* RETA;
		e[i2]=( qx* data.t1x[i]+ qy* data.t1y[i]+ qz* data.t1z[i])* tt1;
		e[i1]=( qx* data.t2x[i]+ qy* data.t2y[i]+ qz* data.t2z[i])* tt1;
	  }

	  if( gnd.ksymp == 1)
		return;

	  tt1=( qy* cph- qx* sph)*( rrv- rrh);
	  cx=-( rrh* qx- tt1* sph);
	  cy=-( rrh* qy+ tt1* cph);
	  cz= rrh* qz;

	  i= -1;
	  i1= data.n-2;
	  for( is = 0; is < data.m; is++ )
	  {
		i++;
		i1 += 2;
		i2 = i1+1;
		arg=- TP*( wx* data.px[i]+ wy* data.py[i]- wz* data.pz[i]);
		tt1= cmplx( cosl( arg), sinl( arg))* data.psalp[i]* RETA;
		e[i2]= e[i2]+( cx* data.t1x[i]+ cy* data.t1y[i]+ cz* data.t1z[i])* tt1;
		e[i1]= e[i1]+( cx* data.t2x[i]+ cy* data.t2y[i]+ cz* data.t2z[i])* tt1;
	  }
	  return;

	} /* if( ipr == 1) */

	/* incident plane wave, elliptic polarization. */
	tt1=-(CPLX_01)* p6;
	if( ipr == 3)
	  tt1=- tt1;

	if( data.n != 0)
	{
	  cx= pxl+ tt1* qx;
	  cy= pyl+ tt1* qy;
	  cz= pzl+ tt1* qz;

	  for( i = 0; i < data.n; i++ )
	  {
		arg=- TP*( wx* data.x[i]+ wy* data.y[i]+ wz* data.z[i]);
		e[i]=-( cx* data.cab[i]+ cy* data.sab[i]+ cz*
			data.salp[i])* cmplx( cosl( arg), sinl( arg));
	  }

	  if( gnd.ksymp != 1)
	  {
		tt2=( cy* cph- cx* sph)*( rrh- rrv);
		cx= rrv* cx- tt2* sph;
		cy= rrv* cy+ tt2* cph;
		cz=- rrv* cz;

		for( i = 0; i < data.n; i++ )
		{
		  arg=- TP*( wx* data.x[i]+ wy* data.y[i]- wz* data.z[i]);
		  e[i]= e[i]-( cx* data.cab[i]+ cy* data.sab[i]+
			  cz* data.salp[i])* cmplx(cosl( arg), sinl( arg));
		}

	  } /* if( gnd.ksymp != 1) */

	} /* if( n != 0) */

	if( data.m == 0)
	  return;

	cx= qx- tt1* pxl;
	cy= qy- tt1* pyl;
	cz= qz- tt1* pzl;

	i= -1;
	i1= data.n-2;
	for( is = 0; is < data.m; is++ )
	{
	  i++;
	  i1 += 2;
	  i2 = i1+1;
	  arg=- TP*( wx* data.px[i]+ wy* data.py[i]+ wz* data.pz[i]);
	  tt2= cmplx( cosl( arg), sinl( arg))* data.psalp[i]* RETA;
	  e[i2]=( cx* data.t1x[i]+ cy* data.t1y[i]+ cz* data.t1z[i])* tt2;
	  e[i1]=( cx* data.t2x[i]+ cy* data.t2y[i]+ cz* data.t2z[i])* tt2;
	}

	if( gnd.ksymp == 1)
	  return;

	tt1=( cy* cph- cx* sph)*( rrv- rrh);
	cx=-( rrh* cx- tt1* sph);
	cy=-( rrh* cy+ tt1* cph);
	cz= rrh* cz;

	i= -1;
	i1= data.n-2;
	for( is=0; is < data.m; is++ )
	{
	  i++;
	  i1 += 2;
	  i2 = i1+1;
	  arg=- TP*( wx* data.px[i]+ wy* data.py[i]- wz* data.pz[i]);
	  tt1= cmplx( cosl( arg), sinl( arg))* data.psalp[i]* RETA;
	  e[i2]= e[i2]+( cx* data.t1x[i]+ cy* data.t1y[i]+ cz* data.t1z[i])* tt1;
	  e[i1]= e[i1]+( cx* data.t2x[i]+ cy* data.t2y[i]+ cz* data.t2z[i])* tt1;
	}

	return;

  } /* if( ipr <= 3) */

  /* incident field of an elementary current source. */
  wz= cosl( p4);
  wx= wz* cosl( p5);
  wy= wz* sinl( p5);
  wz= sinl( p4);
  ds= p6*59.958;
  dsh= p6/(2.* TP);

  is= 0;
  i1= data.n-2;
  for( i = 0; i < data.npm; i++ )
  {
	if( i >= data.n )
	{
	  i1 += 2;
	  i2 = i1+1;
	  pxl= data.px[is]- p1;
	  pyl= data.py[is]- p2;
	  pzl= data.pz[is]- p3;
	}
	else
	{
	  pxl= data.x[i]- p1;
	  pyl= data.y[i]- p2;
	  pzl= data.z[i]- p3;
	}

	rs= pxl* pxl+ pyl* pyl+ pzl* pzl;
	if( rs < 1.0e-30)
	  continue;

	r= sqrtl( rs);
	pxl= pxl/ r;
	pyl= pyl/ r;
	pzl= pzl/ r;
	cth= pxl* wx+ pyl* wy+ pzl* wz;
	sth= sqrtl(1.- cth* cth);
	qx= pxl- wx* cth;
	qy= pyl- wy* cth;
	qz= pzl- wz* cth;

	arg= sqrtl( qx* qx+ qy* qy+ qz* qz);
	if( arg >= 1.e-30)
	{
	  qx= qx/ arg;
	  qy= qy/ arg;
	  qz= qz/ arg;
	}
	else
	{
	  qx=1.;
	  qy=0.;
	  qz=0.;

	} /* if( arg >= 1.e-30) */

	arg=- TP* r;
	tt1= cmplx( cosl( arg), sinl( arg));

	if( i < data.n )
	{
	  tt2= cmplx(1.0,-1.0/( r* TP))/ rs;
	  er= ds* tt1* tt2* cth;
	  et=.5* ds* tt1*((CPLX_01)* TP/ r+ tt2)* sth;
	  ezh= er* cth- et* sth;
	  erh= er* sth+ et* cth;
	  cx= ezh* wx+ erh* qx;
	  cy= ezh* wy+ erh* qy;
	  cz= ezh* wz+ erh* qz;
	  e[i]=-( cx* data.cab[i]+ cy* data.sab[i]+ cz* data.salp[i]);
	}
	else
	{
	  pxl= wy* qz- wz* qy;
	  pyl= wz* qx- wx* qz;
	  pzl= wx* qy- wy* qx;
	  tt2= dsh* tt1* cmplx(1./ r, TP)/ r* sth* data.psalp[is];
	  cx= tt2* pxl;
	  cy= tt2* pyl;
	  cz= tt2* pzl;
	  e[i2]= cx* data.t1x[is]+ cy* data.t1y[is]+ cz* data.t1z[is];
	  e[i1]= cx* data.t2x[is]+ cy* data.t2y[is]+ cz* data.t2z[is];
	  is++;
	} /* if( i < data.n) */

  } /* for( i = 0; i < npm; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine to factor a matrix into a unit lower triangular matrix */
/* and an upper triangular matrix using the gauss-doolittle algorithm */
/* presented on pages 411-416 of a. ralston--a first course in */
/* numerical analysis.  comments below refer to comments in ralstons */
/* text.    (matrix transposed.) */

void factr( int n, complex long double *a, int *ip, int ndim)
{ int r, rm1, rp1, pj, pr, iflg, k, j, jp1, i;
  long double dmax, elmag;
  complex long double arj, *scm = NULL;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, data.np2m * sizeof(complex long double) );

  /* Un-transpose the matrix for Gauss elimination */
  for( i = 1; i < n; i++ )
	for( j = 0; j < i; j++ )
	{
	  arj = a[i+j*ndim];
	  a[i+j*ndim] = a[j+i*ndim];
	  a[j+i*ndim] = arj;
	}

  iflg=FALSE;
  /* step 1 */
  for( r = 0; r < n; r++ )
  {
	for( k = 0; k < n; k++ )
	  scm[k]= a[k+r*ndim];

	/* steps 2 and 3 */
	rm1= r;
	if( rm1 > 0)
	{
	  for( j = 0; j < rm1; j++ )
	  {
		pj= ip[j]-1;
		arj= scm[pj];
		a[j+r*ndim]= arj;
		scm[pj]= scm[j];
		jp1= j+1;

		for( i = jp1; i < n; i++ )
		  scm[i] -= a[i+j*ndim]* arj;

	  } /* for( j = 0; j < rm1; j++ ) */

	} /* if( rm1 >= 0.) */

	/* step 4 */
	dmax= creal( scm[r]*conjl(scm[r]) );

	rp1= r+1;
	ip[r]= rp1;
	if( rp1 < n)
	{
	  for( i = rp1; i < n; i++ )
	  {
		elmag= creal( scm[i]* conjl(scm[i]) );
		if( elmag >= dmax)
		{
		  dmax= elmag;
		  ip[r]= i+1;
		}
	  }
	} /* if( rp1 < n) */

	if( dmax < 1.e-10)
	  iflg=TRUE;

	pr= ip[r]-1;
	a[r+r*ndim]= scm[pr];
	scm[pr]= scm[r];

	/* step 5 */
	if( rp1 < n)
	{
	  arj=1./ a[r+r*ndim];

	  for( i = rp1; i < n; i++ )
		a[i+r*ndim]= scm[i]* arj;
	}

	if( iflg == TRUE )
	{
	  fprintf( output_fp,
		  "\n  PIVOT(%d)= %16.8LE", r, dmax );
	  iflg=FALSE;
	}

  } /* for( r=0; r < n; r++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* factrs, for symmetric structure, transforms submatricies to form */
/* matricies of the symmetric modes and calls routine to factor */
/* matricies.  if no symmetry, the routine is called to factor the */
/* complete matrix. */
void factrs( int np, int nrow, complex long double *a, int *ip )
{
  int kk, ka;

  smat.nop = nrow/np;
  for( kk = 0; kk < smat.nop; kk++ )
  {
	ka= kk* np;
	factr( np, &a[ka], &ip[ka], nrow );
  }
  return;
}

/*-----------------------------------------------------------------------*/

/* fblock sets parameters for out-of-core */
/* solution for the primary matrix (a) */
void fblock( int nrow, int ncol, int imax, int ipsym )
{
  int i, j, k, ka, kk;
  long double phaz, arg;
  complex long double deter;

  if( nrow*ncol <= imax)
  {
	matpar.npblk= nrow;
	matpar.nlast= nrow;
	matpar.imat= nrow* ncol;

	if( nrow == ncol)
	{
	  matpar.icase=1;
	  return;
	}
	else
	  matpar.icase=2;

  } /* if( nrow*ncol <= imax) */

  smat.nop = ncol/nrow;
  if( smat.nop*nrow != ncol)
  {
	fprintf( output_fp,
		"\n  SYMMETRY ERROR - NROW: %d NCOL: %d", nrow, ncol );
	stop(-1);
  }

  /* set up smat.ssx matrix for rotational symmetry. */
  if( ipsym <= 0)
  {
	phaz = TP/smat.nop;

	for( i = 1; i < smat.nop; i++ )
	{
	  for( j= i; j < smat.nop; j++ )
	  {
		arg= phaz* (long double)i * (long double)j;
		smat.ssx[i+j*smat.nop]= cmplx( cosl( arg), sinl( arg));
		smat.ssx[j+i*smat.nop]= smat.ssx[i+j*smat.nop];
	  }
	}
	return;

  } /* if( ipsym <= 0) */

  /* set up smat.ssx matrix for plane symmetry */
  kk=1;
  smat.ssx[0]=CPLX_10;

  k = 2;
  for( ka = 1; k != smat.nop; ka++ )
	k *= 2;

  for( k = 0; k < ka; k++ )
  {
	for( i = 0; i < kk; i++ )
	{
	  for( j = 0; j < kk; j++ )
	  {
		deter= smat.ssx[i+j*smat.nop];
		smat.ssx[i+(j+kk)*smat.nop]= deter;
		smat.ssx[i+kk+(j+kk)*smat.nop]=- deter;
		smat.ssx[i+kk+j*smat.nop]= deter;
	  }
	}
	kk *= 2;

  } /* for( k = 0; k < ka; k++ ) */

  return;
}

/*-----------------------------------------------------------------------*/


/* subroutine to solve the matrix equation lu*x=b where l is a unit */
/* lower triangular matrix and u is an upper triangular matrix both */
/* of which are stored in a.  the rhs vector b is input and the */
/* solution is returned through vector b.   (matrix transposed. */
void solve( int n, complex long double *a, int *ip,
	complex long double *b, int ndim )
{
  int i, ip1, j, k, pia;
  complex long double sum, *scm = NULL;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, data.np2m * sizeof(complex long double) );

  /* forward substitution */
  for( i = 0; i < n; i++ )
  {
	pia= ip[i]-1;
	scm[i]= b[pia];
	b[pia]= b[i];
	ip1= i+1;

	if( ip1 < n)
	  for( j = ip1; j < n; j++ )
		b[j] -= a[j+i*ndim]* scm[i];
  }

  /* backward substitution */
  for( k = 0; k < n; k++ )
  {
	i= n-k-1;
	sum=CPLX_00;
	ip1= i+1;

	if( ip1 < n)
	  for( j = ip1; j < n; j++ )
		sum += a[i+j*ndim]* b[j];

	b[i]=( scm[i]- sum)/ a[i+i*ndim];
  }

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine solves, for symmetric structures, handles the */
/* transformation of the right hand side vector and solution */
/* of the matrix eq. */
void solves( complex long double *a, int *ip, complex long double *b,
	int neq, int nrh, int np, int n, int mp, int m)
{
  int npeq, nrow, ic, i, kk, ia, ib, j, k;
  long double fnop, fnorm;
  complex long double  sum, *scm = NULL;

  npeq= np+ 2*mp;
  smat.nop = neq/npeq;
  fnop= smat.nop;
  fnorm=1./ fnop;
  nrow= neq;

  /* Allocate to scratch memory */
  mem_alloc( (void *)&scm, data.np2m * sizeof(complex long double) );

  if( smat.nop != 1)
  {
	for( ic = 0; ic < nrh; ic++ )
	{
	  if( (n != 0) && (m != 0) )
	  {
		for( i = 0; i < neq; i++ )
		  scm[i]= b[i+ic*neq];

		kk=2* mp;
		ia= np-1;
		ib= n-1;
		j= np-1;

		for( k = 0; k < smat.nop; k++ )
		{
		  if( k != 0 )
		  {
			for( i = 0; i < np; i++ )
			{
			  ia++;
			  j++;
			  b[j+ic*neq]= scm[ia];
			}

			if( k == (smat.nop-1) )
			  continue;

		  } /* if( k != 0 ) */

		  for( i = 0; i < kk; i++ )
		  {
			ib++;
			j++;
			b[j+ic*neq]= scm[ib];
		  }

		} /* for( k = 0; k < smat.nop; k++ ) */

	  } /* if( (n != 0) && (m != 0) ) */

	  /* transform matrix eq. rhs vector according to symmetry modes */
	  for( i = 0; i < npeq; i++ )
	  {
		for( k = 0; k < smat.nop; k++ )
		{
		  ia= i+ k* npeq;
		  scm[k]= b[ia+ic*neq];
		}

		sum= scm[0];
		for( k = 1; k < smat.nop; k++ )
		  sum += scm[k];

		b[i+ic*neq]= sum* fnorm;

		for( k = 1; k < smat.nop; k++ )
		{
		  ia= i+ k* npeq;
		  sum= scm[0];

		  for( j = 1; j < smat.nop; j++ )
			sum += scm[j]* conjl( smat.ssx[k+j*smat.nop]);

		  b[ia+ic*neq]= sum* fnorm;
		}

	  } /* for( i = 0; i < npeq; i++ ) */

	} /* for( ic = 0; ic < nrh; ic++ ) */

  } /* if( smat.nop != 1) */

  /* solve each mode equation */
  for( kk = 0; kk < smat.nop; kk++ )
  {
	ia= kk* npeq;
	ib= ia;

	for( ic = 0; ic < nrh; ic++ )
	  solve( npeq, &a[ib], &ip[ia], &b[ia+ic*neq], nrow );

  } /* for( kk = 0; kk < smat.nop; kk++ ) */

  if( smat.nop == 1)
  {
	free_ptr( (void *)&scm );
	return;
  }

  /* inverse transform the mode solutions */
  for( ic = 0; ic < nrh; ic++ )
  {
	for( i = 0; i < npeq; i++ )
	{
	  for( k = 0; k < smat.nop; k++ )
	  {
		ia= i+ k* npeq;
		scm[k]= b[ia+ic*neq];
	  }

	  sum= scm[0];
	  for( k = 1; k < smat.nop; k++ )
		sum += scm[k];

	  b[i+ic*neq]= sum;
	  for( k = 1; k < smat.nop; k++ )
	  {
		ia= i+ k* npeq;
		sum= scm[0];

		for( j = 1; j < smat.nop; j++ )
		  sum += scm[j]* smat.ssx[k+j*smat.nop];

		b[ia+ic*neq]= sum;
	  }

	} /* for( i = 0; i < npeq; i++ ) */

	if( (n == 0) || (m == 0) )
	  continue;

	for( i = 0; i < neq; i++ )
	  scm[i]= b[i+ic*neq];

	kk=2* mp;
	ia= np-1;
	ib= n-1;
	j= np-1;

	for( k = 0; k < smat.nop; k++ )
	{
	  if( k != 0 )
	  {
		for( i = 0; i < np; i++ )
		{
		  ia++;
		  j++;
		  b[ia+ic*neq]= scm[j];
		}

		if( k == smat.nop)
		  continue;

	  } /* if( k != 0 ) */

	  for( i = 0; i < kk; i++ )
	  {
		ib++;
		j++;
		b[ib+ic*neq]= scm[j];
	  }

	} /* for( k = 0; k < smat.nop; k++ ) */

  } /* for( ic = 0; ic < nrh; ic++ ) */

  free_ptr( (void *)&scm );

  return;
}

/*-----------------------------------------------------------------------*/

