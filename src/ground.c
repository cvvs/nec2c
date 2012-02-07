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

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/* common  /data/ */
extern data_t data;

/* common  /dataj/ */
extern dataj_t dataj;

/* common  /incom/ */
incom_t incom;;

/* common  /gwav/ */
extern gwav_t gwav;

/* common  /gnd/ */
extern gnd_t gnd;

/*-------------------------------------------------------------------*/

/* segment to obtain the total field due to ground.  the method of */
/* variable interval width romberg integration is used.  there are 9 */
/* field components - the x, y, and z components due to constant, */
/* sine, and cosine current distributions. */
void rom2( long double a, long double b, complex long double *sum, long double dmin )
{
  int i, ns, nt, flag=TRUE;
  int nts = 4, nx = 1, n = 9;
  long double ze, ep, zend, dz=0., dzot=0., tmag1, tmag2, tr, ti;
  long double z, s; /***also global***/
  long double rx = 1.0e-4;
  complex long double g1[9], g2[9], g3[9], g4[9], g5[9];
  complex long double t00, t01[9], t10[9], t02, t11, t20[9];

  z= a;
  ze= b;
  s= b- a;

  if( s < 0.)
  {
	fprintf( output_fp, "\n  ERROR - B LESS THAN A IN ROM2" );
	stop(-1);
  }

  ep= s/(1.e4* data.npm);
  zend= ze- ep;

  for( i = 0; i < n; i++ )
	sum[i]=CPLX_00;

  ns= nx;
  nt=0;
  sflds( z, g1);

  while( TRUE )
  {
	if( flag )
	{
	  dz= s/ ns;
	  if( z+ dz > ze)
	  {
		dz= ze- z;
		if( dz <= ep)
		  return;
	  }

	  dzot= dz*.5;
	  sflds( z+ dzot, g3);
	  sflds( z+ dz, g5);

	} /* if( flag ) */

	tmag1=0.;
	tmag2=0.;

	/* evaluate 3 point romberg result and test convergence. */
	for( i = 0; i < n; i++ )
	{
	  t00=( g1[i]+ g5[i])* dzot;
	  t01[i]=( t00+ dz* g3[i])*.5;
	  t10[i]=(4.* t01[i]- t00)/3.;
	  if( i > 2)
		continue;

	  tr= creal( t01[i]);
	  ti= cimag( t01[i]);
	  tmag1= tmag1+ tr* tr+ ti* ti;
	  tr= creal( t10[i]);
	  ti= cimag( t10[i]);
	  tmag2= tmag2+ tr* tr+ ti* ti;

	} /* for( i = 0; i < n; i++ ) */

	tmag1= sqrtl( tmag1);
	tmag2= sqrtl( tmag2);
	test( tmag1, tmag2, &tr, 0., 0., &ti, dmin);

	if( tr <= rx)
	{
	  for( i = 0; i < n; i++ )
		sum[i] += t10[i];
	  nt += 2;

	  z += dz;
	  if( z > zend)
		return;

	  for( i = 0; i < n; i++ )
		g1[i]= g5[i];

	  if( (nt >= nts) && (ns > nx) )
	  {
		ns= ns/2;
		nt=1;
	  }
	  flag = TRUE;
	  continue;

	} /* if( tr <= rx) */

	sflds( z+ dz*.25, g2);
	sflds( z+ dz*.75, g4);
	tmag1=0.;
	tmag2=0.;

	/* evaluate 5 point romberg result and test convergence. */
	for( i = 0; i < n; i++ )
	{
	  t02=( t01[i]+ dzot*( g2[i]+ g4[i]))*.5;
	  t11=( 4.0 * t02- t01[i] )/3.;
	  t20[i]=(16.* t11- t10[i])/15.;
	  if( i > 2)
		continue;

	  tr= creal( t11);
	  ti= cimag( t11);
	  tmag1= tmag1+ tr* tr+ ti* ti;
	  tr= creal( t20[i]);
	  ti= cimag( t20[i]);
	  tmag2= tmag2+ tr* tr+ ti* ti;

	} /* for( i = 0; i < n; i++ ) */

	tmag1= sqrtl( tmag1);
	tmag2= sqrtl( tmag2);
	test( tmag1, tmag2, &tr, 0.,0., &ti, dmin);

	if( tr > rx)
	{
	  nt=0;
	  if( ns < data.npm )
	  {
		ns= ns*2;
		dz= s/ ns;
		dzot= dz*.5;

		for( i = 0; i < n; i++ )
		{
		  g5[i]= g3[i];
		  g3[i]= g2[i];
		}

		flag=FALSE;
		continue;

	  } /* if( ns < npm) */

	  fprintf( output_fp,
		  "\n  ROM2 -- STEP SIZE LIMITED AT Z = %12.5LE", z );

	} /* if( tr > rx) */

	for( i = 0; i < n; i++ )
	  sum[i]= sum[i]+ t20[i];
	nt= nt+1;

	z= z+ dz;
	if( z > zend)
	  return;

	for( i = 0; i < n; i++ )
	  g1[i]= g5[i];

	flag = TRUE;
	if( (nt < nts) || (ns <= nx) )
	  continue;

	ns= ns/2;
	nt=1;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* sfldx returns the field due to ground for a current element on */
/* the source segment at t relative to the segment center. */
void sflds( long double t, complex long double *e )
{
  long double xt, yt, zt, rhx, rhy, rhs, rho, phx, phy;
  long double cph, sph, zphs, r2s, rk, sfac, thet;
  complex long double  erv, ezv, erh, ezh, eph, er, et, hrv, hzv, hrh;

  xt= dataj.xj+ t* dataj.cabj;
  yt= dataj.yj+ t* dataj.sabj;
  zt= dataj.zj+ t* dataj.salpj;
  rhx= incom.xo- xt;
  rhy= incom.yo- yt;
  rhs= rhx* rhx+ rhy* rhy;
  rho= sqrtl( rhs);

  if( rho <= 0.)
  {
	rhx=1.;
	rhy=0.;
	phx=0.;
	phy=1.;
  }
  else
  {
	rhx= rhx/ rho;
	rhy= rhy/ rho;
	phx=- rhy;
	phy= rhx;
  }

  cph= rhx* incom.xsn+ rhy* incom.ysn;
  sph= rhy* incom.xsn- rhx* incom.ysn;

  if( fabsl( cph) < 1.0e-10)
	cph=0.;
  if( fabsl( sph) < 1.0e-10)
	sph=0.;

  gwav.zph= incom.zo+ zt;
  zphs= gwav.zph* gwav.zph;
  r2s= rhs+ zphs;
  gwav.r2= sqrtl( r2s);
  rk= gwav.r2* TP;
  gwav.xx2= cmplx( cosl( rk),- sinl( rk));

  /* use norton approximation for field due to ground.  current is */
  /* lumped at segment center with current moment for constant, sine, */
  /* or cosine distribution. */
  if( incom.isnor != 1)
  {
	gwav.zmh=1.;
	gwav.r1=1.;
	gwav.xx1=0.;
	gwave( &erv, &ezv, &erh, &ezh, &eph);

	et=-CONST1* gnd.frati* gwav.xx2/( r2s* gwav.r2);
	er=2.* et* cmplx(1.0, rk);
	et= et* cmplx(1.0 - rk* rk, rk);
	hrv=( er+ et)* rho* gwav.zph/ r2s;
	hzv=( zphs* er- rhs* et)/ r2s;
	hrh=( rhs* er- zphs* et)/ r2s;
	erv= erv- hrv;
	ezv= ezv- hzv;
	erh= erh+ hrh;
	ezh= ezh+ hrv;
	eph= eph+ et;
	erv= erv* dataj.salpj;
	ezv= ezv* dataj.salpj;
	erh= erh* incom.sn* cph;
	ezh= ezh* incom.sn* cph;
	eph= eph* incom.sn* sph;
	erh= erv+ erh;
	e[0]=( erh* rhx+ eph* phx)* dataj.s;
	e[1]=( erh* rhy+ eph* phy)* dataj.s;
	e[2]=( ezv+ ezh)* dataj.s;
	e[3]=0.;
	e[4]=0.;
	e[5]=0.;
	sfac= PI* dataj.s;
	sfac= sinl( sfac)/ sfac;
	e[6]= e[0]* sfac;
	e[7]= e[1]* sfac;
	e[8]= e[2]* sfac;

	return;
  } /* if( smat.isnor != 1) */

  /* interpolate in sommerfeld field tables */
  if( rho >= 1.0e-12)
	thet= atanl( gwav.zph/ rho);
  else
	thet= POT;

  /* combine vertical and horizontal components and convert */
  /* to x,y,z components. multiply by exp(-jkr)/r. */
  intrp( gwav.r2, thet, &erv, &ezv, &erh, &eph );
  gwav.xx2= gwav.xx2/ gwav.r2;
  sfac= incom.sn* cph;
  erh= gwav.xx2*( dataj.salpj* erv+ sfac* erh);
  ezh= gwav.xx2*( dataj.salpj* ezv- sfac* erv);
  /* x,y,z fields for constant current */
  eph= incom.sn* sph* gwav.xx2* eph;
  e[0]= erh* rhx+ eph* phx;
  e[1]= erh* rhy+ eph* phy;
  e[2]= ezh;
  /* x,y,z fields for sine current */
  rk= TP* t;
  sfac= sinl( rk);
  e[3]= e[0]* sfac;
  e[4]= e[1]* sfac;
  /* x,y,z fields for cosine current */
  e[5]= e[2]* sfac;
  sfac= cosl( rk);
  e[6]= e[0]* sfac;
  e[7]= e[1]* sfac;
  e[8]= e[2]* sfac;

  return;
}

/*-----------------------------------------------------------------------*/

