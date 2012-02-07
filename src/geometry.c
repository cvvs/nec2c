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

/* common  /data/ */
extern data_t data;

/* common  /segj/ */
extern segj_t segj;

/* pointers to input/output files */
extern FILE *input_fp, *output_fp, *plot_fp;

/* common  /plot/ */
extern plot_t plot;

/*-------------------------------------------------------------------*/

/* arc generates segment geometry data for an arc of ns segments */
void arc( int itg, int ns, long double rada,
	long double ang1, long double ang2, long double rad )
{
  int ist, i, mreq;
  long double ang, dang, xs1, xs2, zs1, zs2;

  ist= data.n;
  data.n += ns;
  data.np= data.n;
  data.mp= data.m;
  data.ipsym=0;

  if( ns < 1)
	return;

  if( fabsl( ang2- ang1) < 360.00001)
  {
	/* Reallocate tags buffer */
	mem_realloc( (void *)&data.itag, (data.n+data.m) * sizeof(int) );

	/* Reallocate wire buffers */
	mreq = data.n * sizeof(long double);
	mem_realloc( (void *)&data.x1, mreq );
	mem_realloc( (void *)&data.y1, mreq );
	mem_realloc( (void *)&data.z1, mreq );
	mem_realloc( (void *)&data.x2, mreq );
	mem_realloc( (void *)&data.y2, mreq );
	mem_realloc( (void *)&data.z2, mreq );
	mem_realloc( (void *)&data.bi, mreq );

	ang= ang1* TA;
	dang=( ang2- ang1)* TA/ ns;
	xs1= rada* cosl( ang);
	zs1= rada* sinl( ang);

	for( i = ist; i < data.n; i++ )
	{
	  ang += dang;
	  xs2= rada* cosl( ang);
	  zs2= rada* sinl( ang);
	  data.x1[i]= xs1;

	  data.y1[i]=0.;
	  data.z1[i]= zs1;
	  data.x2[i]= xs2;
	  data.y2[i]=0.;
	  data.z2[i]= zs2;
	  xs1= xs2;
	  zs1= zs2;
	  data.bi[i]= rad;
	  data.itag[i]= itg;

	} /* for( i = ist; i < data.n; i++ ) */

  } /* if( fabsl( ang2- ang1) < 360.00001) */
  else
  {
	fprintf( output_fp, "\n  ERROR -- ARC ANGLE EXCEEDS 360 DEGREES");
	stop(-1);
  }

  return;
}

/*-----------------------------------------------------------------------*/

/* connect sets up segment connection data in arrays icon1 and */
/* icon2 by searching for segment ends that are in contact. */
void conect( int ignd )
{
  int i, iz, ic, j, jx, ix, ixx, iseg, iend, jend, jump, ipf;
  //int nsflg;
  long double sep=0., xi1, yi1, zi1, xi2, yi2, zi2;
  long double slen, xa, ya, za, xs, ys, zs;

  segj.maxcon = 1;

  if( ignd != 0)
  {
	fprintf( output_fp, "\n\n     GROUND PLANE SPECIFIED." );

	if( ignd > 0)
	  fprintf( output_fp,
		  "\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL"
		  " BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n" );

	if( data.ipsym == 2)
	{
	  data.np=2* data.np;
	  data.mp=2* data.mp;
	}

	if( abs( data.ipsym) > 2 )
	{
	  data.np= data.n;
	  data.mp= data.m;
	}

	/*** possibly should be error condition?? **/
	if( data.np > data.n)
	{
	  fprintf( output_fp,
		  "\n ERROR: NP > N IN CONECT()" );
	  stop(-1);
	}

	if( (data.np == data.n) && (data.mp == data.m) )
	  data.ipsym=0;

  } /* if( ignd != 0) */

  if( data.n != 0)
  {
	/* Allocate memory to connections */
	mem_realloc( (void *)&data.icon1, (data.n+data.m) * sizeof(int) );
	mem_realloc( (void *)&data.icon2, (data.n+data.m) * sizeof(int) );

	for( i = 0; i < data.n; i++ )
	{
	  data.icon1[i] = data.icon2[i] = 0;
	  iz = i+1;
	  xi1= data.x1[i];
	  yi1= data.y1[i];
	  zi1= data.z1[i];
	  xi2= data.x2[i];
	  yi2= data.y2[i];
	  zi2= data.z2[i];
	  slen= sqrtl( (xi2- xi1)*(xi2- xi1) + (yi2- yi1) *
		  (yi2- yi1) + (zi2- zi1)*(zi2- zi1) ) * SMIN;

	  /* determine connection data for end 1 of segment. */
	  jump = FALSE;
	  if( ignd > 0)
	  {
		if( zi1 <= -slen)
		{
		  fprintf( output_fp,
			  "\n  GEOMETRY DATA ERROR -- SEGMENT"
			  " %d EXTENDS BELOW GROUND", iz );
		  stop(-1);
		}

		if( zi1 <= slen)
		{
		  data.icon1[i]= iz;
		  data.z1[i]=0.;
		  jump = TRUE;

		} /* if( zi1 <= slen) */

	  } /* if( ignd > 0) */

	  if( ! jump )
	  {
		ic= i;
		for( j = 1; j < data.n; j++)
		{
		  ic++;
		  if( ic >= data.n)
			ic=0;

		  sep= fabsl( xi1- data.x1[ic])+ fabsl(yi1- data.y1[ic])+ fabsl(zi1- data.z1[ic]);
		  if( sep <= slen)
		  {
			data.icon1[i]= -(ic+1);
			break;
		  }

		  sep= fabsl( xi1- data.x2[ic])+ fabsl(yi1- data.y2[ic])+ fabsl(zi1- data.z2[ic]);
		  if( sep <= slen)
		  {
			data.icon1[i]= (ic+1);
			break;
		  }

		} /* for( j = 1; j < data.n; j++) */

	  } /* if( ! jump ) */

	  /* determine connection data for end 2 of segment. */
	  if( (ignd > 0) || jump )
	  {
		if( zi2 <= -slen)
		{
		  fprintf( output_fp,
			  "\n  GEOMETRY DATA ERROR -- SEGMENT"
			  " %d EXTENDS BELOW GROUND", iz );
		  stop(-1);
		}

		if( zi2 <= slen)
		{
		  if( data.icon1[i] == iz )
		  {
			fprintf( output_fp,
				"\n  GEOMETRY DATA ERROR -- SEGMENT"
				" %d LIES IN GROUND PLANE", iz );
			stop(-1);
		  }

		  data.icon2[i]= iz;
		  data.z2[i]=0.;
		  continue;

		} /* if( zi2 <= slen) */

	  } /* if( ignd > 0) */

	  ic= i;
	  for( j = 1; j < data.n; j++ )
	  {
		ic++;
		if( ic >= data.n)
		  ic=0;

		sep= fabsl(xi2- data.x1[ic])+ fabsl(yi2- data.y1[ic])+ fabsl(zi2- data.z1[ic]);
		if( sep <= slen)
		{
		  data.icon2[i]= (ic+1);
		  break;
		}

		sep= fabsl(xi2- data.x2[ic])+ fabsl(yi2- data.y2[ic])+ fabsl(zi2- data.z2[ic]);
		if( sep <= slen)
		{
		  data.icon2[i]= -(ic+1);
		  break;
		}

	  } /* for( j = 1; j < data.n; j++ ) */

	} /* for( i = 0; i < data.n; i++ ) */

	/* find wire-surface connections for new patches */
	if( data.m != 0)
	{
	  ix = -1;
	  i = 0;
	  while( ++i <= data.m )
	  {
		ix++;
		xs= data.px[ix];
		ys= data.py[ix];
		zs= data.pz[ix];

		for( iseg = 0; iseg < data.n; iseg++ )
		{
		  xi1= data.x1[iseg];
		  yi1= data.y1[iseg];
		  zi1= data.z1[iseg];
		  xi2= data.x2[iseg];
		  yi2= data.y2[iseg];
		  zi2= data.z2[iseg];

		  /* for first end of segment */
		  slen=( fabsl(xi2- xi1)+ fabsl(yi2- yi1)+ fabsl(zi2- zi1))* SMIN;
		  sep= fabsl(xi1- xs)+ fabsl(yi1- ys)+ fabsl(zi1- zs);

		  /* connection - divide patch into 4 patches at present array loc. */
		  if( sep <= slen)
		  {
			data.icon1[iseg]=PCHCON+ i;
			ic=0;
			subph( i, ic );
			break;
		  }

		  sep= fabsl(xi2- xs)+ fabsl(yi2- ys)+ fabsl(zi2- zs);
		  if( sep <= slen)
		  {
			data.icon2[iseg]=PCHCON+ i;
			ic=0;
			subph( i, ic );
			break;
		  }

		} /* for( iseg = 0; iseg < data.n; iseg++ ) */

	  } /* while( ++i <= data.m ) */

	} /* if( data.m != 0) */

  } /* if( data.n != 0) */

  fprintf( output_fp, "\n\n"
	  "     TOTAL SEGMENTS USED: %d   SEGMENTS IN A"
	  " SYMMETRIC CELL: %d   SYMMETRY FLAG: %d",
	  data.n, data.np, data.ipsym );

  if( data.m > 0)
	fprintf( output_fp,	"\n"
		"       TOTAL PATCHES USED: %d   PATCHES"
		" IN A SYMMETRIC CELL: %d",  data.m, data.mp );

  iseg=( data.n+ data.m)/( data.np+ data.mp);
  if( iseg != 1)
  {
	/*** may be error condition?? ***/
	if( data.ipsym == 0 )
	{
	  fprintf( output_fp,
		  "\n  ERROR: IPSYM=0 IN CONECT()" );
	  stop(-1);
	}

	if( data.ipsym < 0 )
	  fprintf( output_fp,
		  "\n  STRUCTURE HAS %d FOLD ROTATIONAL SYMMETRY\n", iseg );
	else
	{
	  ic= iseg/2;
	  if( iseg == 8)
		ic=3;
	  fprintf( output_fp,
		  "\n  STRUCTURE HAS %d PLANES OF SYMMETRY\n", ic );
	} /* if( data.ipsym < 0 ) */

  } /* if( iseg == 1) */

  if( data.n == 0)
	return;

  /* Allocate to connection buffers */
  mem_realloc( (void *)&segj.jco, segj.maxcon * sizeof(int) );

  /* adjust connected seg. ends to exactly coincide.  print junctions */
  /* of 3 or more seg.  also find old seg. connecting to new seg. */
  iseg = 0;
  ipf = FALSE;
  for( j = 0; j < data.n; j++ )
  {
	jx = j+1;
	iend=-1;
	jend=-1;
	ix= data.icon1[j];
	ic=1;
	segj.jco[0]= -jx;
	xa= data.x1[j];
	ya= data.y1[j];
	za= data.z1[j];

	while( TRUE )
	{
	  if( (ix != 0) && (ix != (j+1)) && (ix <= PCHCON) )
	  {
		//nsflg=0;

		do
		{
		  if( ix == 0 )
		  {
			fprintf( output_fp,
				"\n  CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT: %d", ix );
			stop(-1);
		  }

		  if( ix < 0 )
			ix= -ix;
		  else
			jend= -jend;

		  jump = FALSE;

		  if( ix == jx )
			break;

		  if( ix < jx )
		  {
			jump = TRUE;
			break;
		  }

		  /* Record max. no. of connections */
		  ic++;
		  if( ic >= segj.maxcon )
		  {
			segj.maxcon = ic+1;
			mem_realloc( (void *)&segj.jco, segj.maxcon * sizeof(int) );
		  }
		  segj.jco[ic-1]= ix* jend;

		  /*if( ix > 0){
			//nsflg=1;
		  }*/
		  ixx = ix-1;
		  if( jend != 1)
		  {
			xa= xa+ data.x1[ixx];
			ya= ya+ data.y1[ixx];
			za= za+ data.z1[ixx];
			ix= data.icon1[ixx];
			continue;
		  }

		  xa= xa+ data.x2[ixx];
		  ya= ya+ data.y2[ixx];
		  za= za+ data.z2[ixx];
		  ix= data.icon2[ixx];

		} /* do */
		while( ix != 0 );

		if( jump && (iend == 1) )
		  break;
		else
		  if( jump )
		  {
			iend=1;
			jend=1;
			ix= data.icon2[j];
			ic=1;
			segj.jco[0]= jx;
			xa= data.x2[j];
			ya= data.y2[j];
			za= data.z2[j];
			continue;
		  }

		sep= (long double)ic;
		xa= xa/ sep;
		ya= ya/ sep;
		za= za/ sep;

		for( i = 0; i < ic; i++ )
		{
		  ix= segj.jco[i];
		  if( ix <= 0)
		  {
			ix=- ix;
			ixx = ix-1;
			data.x1[ixx]= xa;
			data.y1[ixx]= ya;
			data.z1[ixx]= za;
			continue;
		  }

		  ixx = ix-1;
		  data.x2[ixx]= xa;
		  data.y2[ixx]= ya;
		  data.z2[ixx]= za;

		} /* for( i = 0; i < ic; i++ ) */

		if( ic >= 3)
		{
		  if( ! ipf )
		  {
			fprintf( output_fp, "\n\n"
				"    ---------- MULTIPLE WIRE JUNCTIONS ----------\n"
				"    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)" );
			ipf = TRUE;
		  }

		  iseg++;
		  fprintf( output_fp, "\n   %5d      ", iseg );

		  for( i = 1; i <= ic; i++ )
		  {
			fprintf( output_fp, "%5d", segj.jco[i-1] );
			if( !(i % 20) )
			  fprintf( output_fp, "\n              " );
		  }

		} /* if( ic >= 3) */

	  } /*if( (ix != 0) && (ix != j) && (ix <= PCHCON) ) */

	  if( iend == 1)
		break;

	  iend=1;
	  jend=1;
	  ix= data.icon2[j];
	  ic=1;
	  segj.jco[0]= jx;
	  xa= data.x2[j];
	  ya= data.y2[j];
	  za= data.z2[j];

	} /* while( TRUE ) */

  } /* for( j = 0; j < data.n; j++ ) */

  mem_realloc( (void *)&segj.ax, segj.maxcon * sizeof(long double) );
  mem_realloc( (void *)&segj.bx, segj.maxcon * sizeof(long double) );
  mem_realloc( (void *)&segj.cx, segj.maxcon * sizeof(long double) );

  return;
}

/*-----------------------------------------------------------------------*/

/* datagn is the main routine for input of geometry data. */
void datagn( void )
{
  char gm[3];
  char ifx[2] = {'*', 'X'}, ify[2]={'*','Y'}, ifz[2]={'*','Z'};
  char ipt[4] = { 'P', 'R', 'T', 'Q' };

  /* input card mnemonic list */
  /* "XT" stands for "exit", added for testing */
#define GM_NUM  12
  char *atst[GM_NUM] =
  {
	"GW", "GX", "GR", "GS", "GE", "GM", \
	  "SP", "SM", "GA", "SC", "GH", "GF"
  };

  int nwire, isct, iphd, i1, i2, itg, iy, iz, mreq;
  int ix, i, ns, gm_num; /* geometry card id as a number */
  long double rad, xs1, xs2, ys1, ys2, zs1, zs2, x4=0, y4=0, z4=0;
  long double x3=0, y3=0, z3=0, xw1, xw2, yw1, yw2, zw1, zw2;
  long double dummy;

  data.ipsym=0;
  nwire=0;
  data.n=0;
  data.np=0;
  data.m=0;
  data.mp=0;
  isct=0;
  iphd = FALSE;

  /* read geometry data card and branch to */
  /* section for operation requested */
  do
  {
	readgm( gm, &itg, &ns, &xw1, &yw1, &zw1, &xw2, &yw2, &zw2, &rad);

	/* identify card id mnemonic */
	for( gm_num = 0; gm_num < GM_NUM; gm_num++ )
	  if( strncmp( gm, atst[gm_num], 2) == 0 )
		break;

	if( iphd == FALSE )
	{
	  fprintf( output_fp, "\n\n\n"
		  "                               "
		  "-------- STRUCTURE SPECIFICATION --------\n"
		  "                                     "
		  "COORDINATES MUST BE INPUT IN\n"
		  "                                     "
		  "METERS OR BE SCALED TO METERS\n"
		  "                                     "
		  "BEFORE STRUCTURE INPUT IS ENDED\n" );

	  fprintf( output_fp, "\n"
		  "  WIRE                                           "
		  "                                      SEG FIRST  LAST  TAG\n"
		  "   No:        X1         Y1         Z1         X2      "
		  "   Y2         Z2       RADIUS   No:   SEG   SEG  No:" );

	  iphd=1;
	}

	if( gm_num != 10 )
	  isct=0;

	switch( gm_num )
	{
	  case 0: /* "gw" card, generate segment data for straight wire. */

		nwire++;
		i1= data.n+1;
		i2= data.n+ ns;

		fprintf( output_fp, "\n"
			" %5d  %10.4LF %10.4LF %10.4LF %10.4LF"
			" %10.4LF %10.4LF %10.4LF %5d %5d %5d %4d",
			nwire, xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, i1, i2, itg );

		if( rad != 0)
		{
		  xs1=1.;
		  ys1=1.;
		}
		else
		{
		  readgm( gm, &ix, &iy, &xs1, &ys1, &zs1,
			  &dummy, &dummy, &dummy, &dummy);

		  if( strcmp(gm, "GC" ) != 0 )
		  {
			fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
			stop(-1);
		  }

		  fprintf( output_fp,
			  "\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: %9.5LF\n"
			  "                                 "
			  "RADIUS FROM: %9.5LF TO: %9.5LF", xs1, ys1, zs1 );

		  if( (ys1 == 0) || (zs1 == 0) )
		  {
			fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
			stop(-1);
		  }

		  rad= ys1;
		  ys1= powl( (zs1/ys1), (1./(ns-1.)) );
		}

		wire( xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1, ns, itg);

		continue;

		/* reflect structure along x,y, or z */
		/* axes or rotate to form cylinder.  */
	  case 1: /* "gx" card */

		iy= ns/10;
		iz= ns- iy*10;
		ix= iy/10;
		iy= iy- ix*10;

		if( ix != 0)
		  ix=1;
		if( iy != 0)
		  iy=1;
		if( iz != 0)
		  iz=1;

		fprintf( output_fp,
			"\n  STRUCTURE REFLECTED ALONG THE AXES %c %c %c"
			" - TAGS INCREMENTED BY %d\n",
			ifx[ix], ify[iy], ifz[iz], itg );

		reflc( ix, iy, iz, itg, ns);

		continue;

	  case 2: /* "gr" card */

		fprintf( output_fp,
			"\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
			" - LABELS INCREMENTED BY %d\n", ns, itg );

		ix=-1;
		iz = 0;
		reflc( ix, iy, iz, itg, ns);

		continue;

	  case 3: /* "gs" card, scale structure dimensions by factor xw1. */

		if( data.n > 0)
		{
		  for( i = 0; i < data.n; i++ )
		  {
			data.x1[i]= data.x1[i]* xw1;
			data.y1[i]= data.y1[i]* xw1;
			data.z1[i]= data.z1[i]* xw1;
			data.x2[i]= data.x2[i]* xw1;
			data.y2[i]= data.y2[i]* xw1;
			data.z2[i]= data.z2[i]* xw1;
			data.bi[i]= data.bi[i]* xw1;
		  }
		} /* if( data.n >= n2) */

		if( data.m > 0)
		{
		  yw1= xw1* xw1;
		  for( i = 0; i < data.m; i++ )
		  {
			data.px[i]= data.px[i]* xw1;
			data.py[i]= data.py[i]* xw1;
			data.pz[i]= data.pz[i]* xw1;
			data.pbi[i]= data.pbi[i]* yw1;
		  }
		} /* if( data.m >= m2) */

		fprintf( output_fp,
			"\n     STRUCTURE SCALED BY FACTOR: %10.5LF", xw1 );

		continue;

	  case 4: /* "ge" card, terminate structure geometry input. */

		if( ns != 0)
		{
		  plot.iplp1=1;
		  plot.iplp2=1;
		}

		conect( itg);

		if( data.n != 0)
		{
		  /* Allocate wire buffers */
		  mreq = data.n * sizeof(long double);
		  mem_realloc( (void *)&data.si, mreq );
		  mem_realloc( (void *)&data.sab, mreq );
		  mem_realloc( (void *)&data.cab, mreq );
		  mem_realloc( (void *)&data.salp, mreq );
		  mem_realloc( (void *)&data.x, mreq );
		  mem_realloc( (void *)&data.y, mreq );
		  mem_realloc( (void *)&data.z, mreq );

		  fprintf( output_fp, "\n\n\n"
			  "                              "
			  " ---------- SEGMENTATION DATA ----------\n"
			  "                                       "
			  " COORDINATES IN METERS\n"
			  "                           "
			  " I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n" );

		  fprintf( output_fp, "\n"
			  "   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION"
			  " ANGLES    WIRE    CONNECTION DATA   TAG\n"
			  "   No:       X         Y         Z      LENGTH     ALPHA     "
			  " BETA    RADIUS    I-     I    I+   No:" );

		  for( i = 0; i < data.n; i++ )
		  {
			xw1= data.x2[i]- data.x1[i];
			yw1= data.y2[i]- data.y1[i];
			zw1= data.z2[i]- data.z1[i];
			data.x[i]=( data.x1[i]+ data.x2[i])/2.;
			data.y[i]=( data.y1[i]+ data.y2[i])/2.;
			data.z[i]=( data.z1[i]+ data.z2[i])/2.;
			xw2= xw1* xw1+ yw1* yw1+ zw1* zw1;
			yw2= sqrtl( xw2);
			yw2=( xw2/ yw2+ yw2)*.5;
			data.si[i]= yw2;
			data.cab[i]= xw1/ yw2;
			data.sab[i]= yw1/ yw2;
			xw2= zw1/ yw2;

			if( xw2 > 1.)
			  xw2=1.;
			if( xw2 < -1.)
			  xw2=-1.;

			data.salp[i]= xw2;
			xw2= asinl( xw2)* TD;
			yw2= atan2l( yw1, xw1)* TD;

			fprintf( output_fp, "\n"
				" %5d %9.4LF %9.4LF %9.4LF %9.4LF"
				" %9.4LF %9.4LF %9.4LF %5d %5d %5d %5d",
				i+1, data.x[i], data.y[i], data.z[i], data.si[i], xw2, yw2,
				data.bi[i], data.icon1[i], i+1, data.icon2[i], data.itag[i] );

			if( plot.iplp1 == 1)
			  fprintf( plot_fp, "%12.4LE %12.4LE %12.4LE "
				  "%12.4LE %12.4LE %12.4LE %12.4LE %5d %5d %5d\n",
				  data.x[i],data.y[i],data.z[i],data.si[i],xw2,yw2,
				  data.bi[i],data.icon1[i],i+1,data.icon2[i] );

			if( (data.si[i] <= 1.e-20) || (data.bi[i] <= 0.) )
			{
			  fprintf( output_fp, "\n SEGMENT DATA ERROR" );
			  stop(-1);
			}

		  } /* for( i = 0; i < data.n; i++ ) */

		} /* if( data.n != 0) */

		if( data.m != 0)
		{
		  fprintf( output_fp, "\n\n\n"
			  "                                   "
			  " --------- SURFACE PATCH DATA ---------\n"
			  "                                            "
			  " COORDINATES IN METERS\n\n"
			  " PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      "
			  " PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n"
			  "  No:       X          Y          Z          X        Y        Z      "
			  " AREA         X1       Y1       Z1        X2       Y2      Z2" );

		  for( i = 0; i < data.m; i++ )
		  {
			xw1=( data.t1y[i]* data.t2z[i]- data.t1z[i]* data.t2y[i])* data.psalp[i];
			yw1=( data.t1z[i]* data.t2x[i]- data.t1x[i]* data.t2z[i])* data.psalp[i];
			zw1=( data.t1x[i]* data.t2y[i]- data.t1y[i]* data.t2x[i])* data.psalp[i];

			fprintf( output_fp, "\n"
				" %4d %10.5LF %10.5LF %10.5LF  %8.4LF %8.4LF %8.4LF"
				" %10.5LF  %8.4LF %8.4LF %8.4LF  %8.4LF %8.4LF %8.4LF",
				i+1, data.px[i], data.py[i], data.pz[i], xw1, yw1, zw1, data.pbi[i],
				data.t1x[i], data.t1y[i], data.t1z[i], data.t2x[i], data.t2y[i], data.t2z[i] );

		  } /* for( i = 0; i < data.m; i++ ) */

		} /* if( data.m == 0) */

		data.npm  = data.n+data.m;
		data.np2m = data.n+2*data.m;
		data.np3m = data.n+3*data.m;

		return;

		/* "gm" card, move structure or reproduce */
		/* original structure in new positions.   */
	  case 5:

		fprintf( output_fp,
			"\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n"
			"   %3d %5d %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
			itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

		xw1= xw1* TA;
		yw1= yw1* TA;
		zw1= zw1* TA;

		move( xw1, yw1, zw1, xw2, yw2, zw2, (int)( rad+.5), ns, itg);
		continue;

	  case 6: /* "sp" card, generate single new patch */

		i1= data.m+1;
		ns++;

		if( itg != 0)
		{
		  fprintf( output_fp, "\n  PATCH DATA ERROR" );
		  stop(-1);
		}

		fprintf( output_fp, "\n"
			" %5d%c %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
			i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

		if( (ns == 2) || (ns == 4) )
		  isct=1;

		if( ns > 1)
		{
		  readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

		  if( (ns == 2) || (itg > 0) )
		  {
			x4= xw1+ x3- xw2;
			y4= yw1+ y3- yw2;
			z4= zw1+ z3- zw2;
		  }

		  fprintf( output_fp, "\n"
			  "      %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
			  x3, y3, z3, x4, y4, z4 );

		  if( strcmp(gm, "SC") != 0 )
		  {
			fprintf( output_fp, "\n  PATCH DATA ERROR" );
			stop(-1);
		  }

		} /* if( ns > 1) */
		else
		{
		  xw2= xw2* TA;
		  yw2= yw2* TA;
		}

		patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

		continue;

	  case 7: /* "sm" card, generate multiple-patch surface */

		i1= data.m+1;
		fprintf( output_fp, "\n"
			" %5d%c %10.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF"
			"     SURFACE - %d BY %d PATCHES",
			i1, ipt[1], xw1, yw1, zw1, xw2, yw2, zw2, itg, ns );

		if( (itg < 1) || (ns < 1) )
		{
		  fprintf( output_fp, "\n  PATCH DATA ERROR" );
		  stop(-1);
		}

		readgm( gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);

		if( (ns == 2) || (itg > 0) )
		{
		  x4= xw1+ x3- xw2;
		  y4= yw1+ y3- yw2;
		  z4= zw1+ z3- zw2;
		}

		fprintf( output_fp, "\n"
			"      %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
			x3, y3, z3, x4, y4, z4 );

		if( strcmp(gm, "SC" ) != 0 )
		{
		  fprintf( output_fp, "\n  PATCH DATA ERROR" );
		  stop(-1);
		}

		patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

		continue;

	  case 8: /* "ga" card, generate segment data for wire arc */

		nwire++;
		i1= data.n+1;
		i2= data.n+ ns;

		fprintf( output_fp, "\n"
			" %5d  ARC RADIUS: %9.5LF  FROM: %8.3LF TO: %8.3LF DEGREES"
			"       %11.5LF %5d %5d %5d %4d",
			nwire, xw1, yw1, zw1, xw2, ns, i1, i2, itg );

		arc( itg, ns, xw1, yw1, zw1, xw2);

		continue;

	  case 9: /* "sc" card */

		if( isct == 0)
		{
		  fprintf( output_fp, "\n  PATCH DATA ERROR" );
		  stop(-1);
		}

		i1= data.m+1;
		ns++;

		if( (itg != 0) || ((ns != 2) && (ns != 4)) )
		{
		  fprintf( output_fp, "\n  PATCH DATA ERROR" );
		  stop(-1);
		}

		xs1= x4;
		ys1= y4;
		zs1= z4;
		xs2= x3;
		ys2= y3;
		zs2= z3;
		x3= xw1;
		y3= yw1;
		z3= zw1;

		if( ns == 4)
		{
		  x4= xw2;
		  y4= yw2;
		  z4= zw2;
		}

		xw1= xs1;
		yw1= ys1;
		zw1= zs1;
		xw2= xs2;
		yw2= ys2;
		zw2= zs2;

		if( ns != 4)
		{
		  x4= xw1+ x3- xw2;
		  y4= yw1+ y3- yw2;
		  z4= zw1+ z3- zw2;
		}

		fprintf( output_fp, "\n"
			" %5d%c %10.5LF %11.5LF %11.5LF %11.5LF %11.5LF %11.5LF",
			i1, ipt[ns-1], xw1, yw1, zw1, xw2, yw2, zw2 );

		fprintf( output_fp, "\n"
			"      %11.5LF %11.5LF %11.5LF  %11.5LF %11.5LF %11.5LF",
			x3, y3, z3, x4, y4, z4 );

		patch( itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);

		continue;

	  case 10: /* "gh" card, generate helix */

		nwire++;
		i1= data.n+1;
		i2= data.n+ ns;

		fprintf( output_fp, "\n"
			" %5d HELIX STRUCTURE - SPACING OF TURNS: %8.3LF AXIAL"
			" LENGTH: %8.3LF  %8.3LF %5d %5d %5d %4d\n      "
			" RADIUS X1:%8.3LF Y1:%8.3LF X2:%8.3LF Y2:%8.3LF ",
			nwire, xw1, yw1, rad, ns, i1, i2, itg, zw1, xw2, yw2, zw2 );

		helix( xw1, yw1, zw1, xw2, yw2, zw2, rad, ns, itg);

		continue;

	  case 11: /* "gf" card, not supported */
		abort_on_error(-5);

	  default: /* error message */

		fprintf( output_fp, "\n  GEOMETRY DATA CARD ERROR" );
		fprintf( output_fp, "\n"
			" %2s %3d %5d %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF %10.5LF",
			gm, itg, ns, xw1, yw1, zw1, xw2, yw2, zw2, rad );

		stop(-1);

	} /* switch( gm_num ) */

  } /* do */
  while( TRUE );

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine helix generates segment geometry */
/* data for a helix of ns segments */
void helix( long double s, long double hl, long double a1, long double b1,
	long double a2, long double b2, long double rad, int ns, int itg )
{
  int ist, i, mreq;
  //long double turns;
  long double zinc, copy, sangle, hdia, turn, pitch, hmaj, hmin;

  ist= data.n;
  data.n += ns;
  data.np= data.n;
  data.mp= data.m;
  data.ipsym=0;

  if( ns < 1)
	return;

  //turns= fabsl( hl/ s);
  zinc= fabsl( hl/ ns);

  /* Reallocate tags buffer */
  mem_realloc( (void *)&data.itag, (data.n+data.m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = data.n * sizeof(long double);
  mem_realloc( (void *)&data.x1, mreq );
  mem_realloc( (void *)&data.y1, mreq );
  mem_realloc( (void *)&data.z1, mreq );
  mem_realloc( (void *)&data.x2, mreq );
  mem_realloc( (void *)&data.y2, mreq );
  mem_realloc( (void *)&data.z2, mreq );
  mem_realloc( (void *)&data.bi, mreq );

  data.z1[ist]=0.;
  for( i = ist; i < data.n; i++ )
  {
	data.bi[i]= rad;
	data.itag[i]= itg;

	if( i != ist )
	  data.z1[i]= data.z1[i-1]+ zinc;

	data.z2[i]= data.z1[i]+ zinc;

	if( a2 == a1)
	{
	  if( b1 == 0.)
		b1= a1;

	  data.x1[i]= a1* cosl(2.* PI* data.z1[i]/ s);
	  data.y1[i]= b1* sinl(2.* PI* data.z1[i]/ s);
	  data.x2[i]= a1* cosl(2.* PI* data.z2[i]/ s);
	  data.y2[i]= b1* sinl(2.* PI* data.z2[i]/ s);
	}
	else
	{
	  if( b2 == 0.)
		b2= a2;

	  data.x1[i]=( a1+( a2- a1)* data.z1[i]/ fabsl( hl))* cosl(2.* PI* data.z1[i]/ s);
	  data.y1[i]=( b1+( b2- b1)* data.z1[i]/ fabsl( hl))* sinl(2.* PI* data.z1[i]/ s);
	  data.x2[i]=( a1+( a2- a1)* data.z2[i]/ fabsl( hl))* cosl(2.* PI* data.z2[i]/ s);
	  data.y2[i]=( b1+( b2- b1)* data.z2[i]/ fabsl( hl))* sinl(2.* PI* data.z2[i]/ s);

	} /* if( a2 == a1) */

	if( hl > 0.)
	  continue;

	copy= data.x1[i];
	data.x1[i]= data.y1[i];
	data.y1[i]= copy;
	copy= data.x2[i];
	data.x2[i]= data.y2[i];
	data.y2[i]= copy;

  } /* for( i = ist; i < data.n; i++ ) */

  if( a2 != a1)
  {
	sangle= atanl( a2/( fabsl( hl)+( fabsl( hl)* a1)/( a2- a1)));
	fprintf( output_fp,
		"\n       THE CONE ANGLE OF THE SPIRAL IS %10.4LF", sangle );
	return;
  }

  if( a1 == b1)
  {
	hdia=2.* a1;
	turn= hdia* PI;
	pitch= atanl( s/( PI* hdia));
	turn= turn/ cosl( pitch);
	pitch=180.* pitch/ PI;
  }
  else
  {
	if( a1 >= b1)
	{
	  hmaj=2.* a1;
	  hmin=2.* b1;
	}
	else
	{
	  hmaj=2.* b1;
	  hmin=2.* a1;
	}

	hdia= sqrtl(( hmaj*hmaj+ hmin*hmin)/2* hmaj);
	turn=2.* PI* hdia;
	pitch=(180./ PI)* atanl( s/( PI* hdia));

  } /* if( a1 == b1) */

  fprintf( output_fp, "\n"
	  "       THE PITCH ANGLE IS: %.4LF    THE LENGTH OF WIRE/TURN IS: %.4LF",
	  pitch, turn );

  return;
}

/*-----------------------------------------------------------------------*/

/* isegno returns the segment number of the mth segment having the */
/* tag number itagi.  if itagi=0 segment number m is returned. */
int isegno( int itagi, int mx)
{
  int icnt, i, iseg;

  if( mx <= 0)
  {
	fprintf( output_fp,
		"\n  CHECK DATA, PARAMETER SPECIFYING SEGMENT"
		" POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO" );
	stop(-1);
  }

  icnt=0;
  if( itagi == 0)
  {
	iseg = mx;
	return( iseg );
  }

  if( data.n > 0)
  {
	for( i = 0; i < data.n; i++ )
	{
	  if( data.itag[i] != itagi )
		continue;

	  icnt++;
	  if( icnt == mx)
	  {
		iseg= i+1;
		return( iseg );
	  }

	} /* for( i = 0; i < data.n; i++ ) */

  } /* if( data.n > 0) */

  fprintf( output_fp, "\n\n"
	  "  NO SEGMENT HAS AN ITAG OF %d",  itagi );
  stop(-1);

  return(0);
}

/*-----------------------------------------------------------------------*/

/* subroutine move moves the structure with respect to its */
/* coordinate system or reproduces structure in new positions. */
/* structure is rotated about x,y,z axes by rox,roy,roz */
/* respectively, then shifted by xs,ys,zs */
void move( long double rox, long double roy, long double roz, long double xs,
	long double ys, long double zs, int its, int nrpt, int itgi )
{
  int nrp, ix, i1, k, ir, i, ii, mreq;
  long double sps, cps, sth, cth, sph, cph, xx, xy;
  long double xz, yx, yy, yz, zx, zy, zz, xi, yi, zi;

  if( fabsl( rox)+ fabsl( roy) > 1.0e-10)
	data.ipsym= data.ipsym*3;

  sps= sinl( rox);
  cps= cosl( rox);
  sth= sinl( roy);
  cth= cosl( roy);
  sph= sinl( roz);
  cph= cosl( roz);
  xx= cph* cth;
  xy= cph* sth* sps- sph* cps;
  xz= cph* sth* cps+ sph* sps;
  yx= sph* cth;
  yy= sph* sth* sps+ cph* cps;
  yz= sph* sth* cps- cph* sps;
  zx=- sth;
  zy= cth* sps;
  zz= cth* cps;

  if( nrpt == 0)
	nrp=1;
  else
	nrp= nrpt;

  ix=1;
  if( data.n > 0)
  {
	i1= isegno( its, 1);
	if( i1 < 1)
	  i1= 1;

	ix= i1;
	if( nrpt == 0)
	  k= i1-1;
	else
	{
	  k= data.n;
	  /* Reallocate tags buffer */
	  mreq = data.n+data.m + (data.n+1-i1)*nrpt;
	  mem_realloc( (void *)&data.itag, mreq * sizeof(int) );

	  /* Reallocate wire buffers */
	  mreq = (data.n+(data.n+1-i1)*nrpt) * sizeof(long double);
	  mem_realloc( (void *)&data.x1, mreq );
	  mem_realloc( (void *)&data.y1, mreq );
	  mem_realloc( (void *)&data.z1, mreq );
	  mem_realloc( (void *)&data.x2, mreq );
	  mem_realloc( (void *)&data.y2, mreq );
	  mem_realloc( (void *)&data.z2, mreq );
	  mem_realloc( (void *)&data.bi, mreq );
	}

	for( ir = 0; ir < nrp; ir++ )
	{
	  for( i = i1-1; i < data.n; i++ )
	  {
		xi= data.x1[i];
		yi= data.y1[i];
		zi= data.z1[i];
		data.x1[k]= xi* xx+ yi* xy+ zi* xz+ xs;
		data.y1[k]= xi* yx+ yi* yy+ zi* yz+ ys;
		data.z1[k]= xi* zx+ yi* zy+ zi* zz+ zs;
		xi= data.x2[i];
		yi= data.y2[i];
		zi= data.z2[i];
		data.x2[k]= xi* xx+ yi* xy+ zi* xz+ xs;
		data.y2[k]= xi* yx+ yi* yy+ zi* yz+ ys;
		data.z2[k]= xi* zx+ yi* zy+ zi* zz+ zs;
		data.bi[k]= data.bi[i];
		data.itag[k]= data.itag[i];
		if( data.itag[i] != 0)
		  data.itag[k]= data.itag[i]+ itgi;

		k++;

	  } /* for( i = i1; i < data.n; i++ ) */

	  i1= data.n+1;
	  data.n= k;

	} /* for( ir = 0; ir < nrp; ir++ ) */

  } /* if( data.n >= n2) */

  if( data.m > 0)
  {
	i1 = 0;
	if( nrpt == 0)
	  k= 0;
	else
	  k = data.m;

	/* Reallocate patch buffers */
	mreq = data.m * (1+nrpt) * sizeof(long double);
	mem_realloc( (void *)&data.px, mreq );
	mem_realloc( (void *)&data.py, mreq );
	mem_realloc( (void *)&data.pz, mreq );
	mem_realloc( (void *)&data.t1x, mreq );
	mem_realloc( (void *)&data.t1y, mreq );
	mem_realloc( (void *)&data.t1z, mreq );
	mem_realloc( (void *)&data.t2x, mreq );
	mem_realloc( (void *)&data.t2y, mreq );
	mem_realloc( (void *)&data.t2z, mreq );
	mem_realloc( (void *)&data.pbi, mreq );
	mem_realloc( (void *)&data.psalp, mreq );

	for( ii = 0; ii < nrp; ii++ )
	{
	  for( i = i1; i < data.m; i++ )
	  {
		xi= data.px[i];
		yi= data.py[i];
		zi= data.pz[i];
		data.px[k]= xi* xx+ yi* xy+ zi* xz+ xs;
		data.py[k]= xi* yx+ yi* yy+ zi* yz+ ys;
		data.pz[k]= xi* zx+ yi* zy+ zi* zz+ zs;
		xi= data.t1x[i];
		yi= data.t1y[i];
		zi= data.t1z[i];
		data.t1x[k]= xi* xx+ yi* xy+ zi* xz;
		data.t1y[k]= xi* yx+ yi* yy+ zi* yz;
		data.t1z[k]= xi* zx+ yi* zy+ zi* zz;
		xi= data.t2x[i];
		yi= data.t2y[i];
		zi= data.t2z[i];
		data.t2x[k]= xi* xx+ yi* xy+ zi* xz;
		data.t2y[k]= xi* yx+ yi* yy+ zi* yz;
		data.t2z[k]= xi* zx+ yi* zy+ zi* zz;
		data.psalp[k]= data.psalp[i];
		data.pbi[k]= data.pbi[i];
		k++;

	  } /* for( i = i1; i < data.m; i++ ) */

	  i1= data.m;
	  data.m = k;

	} /* for( ii = 0; ii < nrp; ii++ ) */

  } /* if( data.m >= m2) */

  if( (nrpt == 0) && (ix == 1) )
	return;

  data.np= data.n;
  data.mp= data.m;
  data.ipsym=0;

  return;
}

/*-----------------------------------------------------------------------*/

/* patch generates and modifies patch geometry data */
void patch( int nx, int ny,
	long double ax1, long double ay1, long double az1,
	long double ax2, long double ay2, long double az2,
	long double ax3, long double ay3, long double az3,
	long double ax4, long double ay4, long double az4 )
{
  int mi, ntp, iy, ix, mreq;
  long double s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., xst=0.;
  long double znv, xnv, ynv, xa, xn2, yn2, zn2, salpn, xs, ys, zs, xt, yt, zt;

  /* new patches.  for nx=0, ny=1,2,3,4 patch is (respectively) */
  /* arbitrary, rectagular, triangular, or quadrilateral. */
  /* for nx and ny  > 0 a rectangular surface is produced with */
  /* nx by ny rectangular patches. */

  data.m++;
  mi= data.m-1;

  /* Reallocate patch buffers */
  mreq = data.m * sizeof(long double);
  mem_realloc( (void *)&data.px, mreq );
  mem_realloc( (void *)&data.py, mreq );
  mem_realloc( (void *)&data.pz, mreq );
  mem_realloc( (void *)&data.t1x, mreq );
  mem_realloc( (void *)&data.t1y, mreq );
  mem_realloc( (void *)&data.t1z, mreq );
  mem_realloc( (void *)&data.t2x, mreq );
  mem_realloc( (void *)&data.t2y, mreq );
  mem_realloc( (void *)&data.t2z, mreq );
  mem_realloc( (void *)&data.pbi, mreq );
  mem_realloc( (void *)&data.psalp, mreq );

  if( nx > 0)
	ntp=2;
  else
	ntp= ny;

  if( ntp <= 1)
  {
	data.px[mi]= ax1;
	data.py[mi]= ay1;
	data.pz[mi]= az1;
	data.pbi[mi]= az2;
	znv= cosl( ax2);
	xnv= znv* cosl( ay2);
	ynv= znv* sinl( ay2);
	znv= sinl( ax2);
	xa= sqrtl( xnv* xnv+ ynv* ynv);

	if( xa >= 1.0e-6)
	{
	  data.t1x[mi]=- ynv/ xa;
	  data.t1y[mi]= xnv/ xa;
	  data.t1z[mi]=0.;
	}
	else
	{
	  data.t1x[mi]=1.;
	  data.t1y[mi]=0.;
	  data.t1z[mi]=0.;
	}

  } /* if( ntp <= 1) */
  else
  {
	s1x= ax2- ax1;
	s1y= ay2- ay1;
	s1z= az2- az1;
	s2x= ax3- ax2;
	s2y= ay3- ay2;
	s2z= az3- az2;

	if( nx != 0)
	{
	  s1x= s1x/ nx;
	  s1y= s1y/ nx;
	  s1z= s1z/ nx;
	  s2x= s2x/ ny;
	  s2y= s2y/ ny;
	  s2z= s2z/ ny;
	}

	xnv= s1y* s2z- s1z* s2y;
	ynv= s1z* s2x- s1x* s2z;
	znv= s1x* s2y- s1y* s2x;
	xa= sqrtl( xnv* xnv+ ynv* ynv+ znv* znv);
	xnv= xnv/ xa;
	ynv= ynv/ xa;
	znv= znv/ xa;
	xst= sqrtl( s1x* s1x+ s1y* s1y+ s1z* s1z);
	data.t1x[mi]= s1x/ xst;
	data.t1y[mi]= s1y/ xst;
	data.t1z[mi]= s1z/ xst;

	if( ntp <= 2)
	{
	  data.px[mi]= ax1+.5*( s1x+ s2x);
	  data.py[mi]= ay1+.5*( s1y+ s2y);
	  data.pz[mi]= az1+.5*( s1z+ s2z);
	  data.pbi[mi]= xa;
	}
	else
	{
	  if( ntp != 4)
	  {
		data.px[mi]=( ax1+ ax2+ ax3)/3.;
		data.py[mi]=( ay1+ ay2+ ay3)/3.;
		data.pz[mi]=( az1+ az2+ az3)/3.;
		data.pbi[mi]=.5* xa;
	  }
	  else
	  {
		s1x= ax3- ax1;
		s1y= ay3- ay1;
		s1z= az3- az1;
		s2x= ax4- ax1;
		s2y= ay4- ay1;
		s2z= az4- az1;
		xn2= s1y* s2z- s1z* s2y;
		yn2= s1z* s2x- s1x* s2z;
		zn2= s1x* s2y- s1y* s2x;
		xst= sqrtl( xn2* xn2+ yn2* yn2+ zn2* zn2);
		salpn=1./(3.*( xa+ xst));
		data.px[mi]=( xa*( ax1+ ax2+ ax3)+ xst*( ax1+ ax3+ ax4))* salpn;
		data.py[mi]=( xa*( ay1+ ay2+ ay3)+ xst*( ay1+ ay3+ ay4))* salpn;
		data.pz[mi]=( xa*( az1+ az2+ az3)+ xst*( az1+ az3+ az4))* salpn;
		data.pbi[mi]=.5*( xa+ xst);
		s1x=( xnv* xn2+ ynv* yn2+ znv* zn2)/ xst;

		if( s1x <= 0.9998)
		{
		  fprintf( output_fp,
			  "\n  ERROR -- CORNERS OF QUADRILATERAL"
			  " PATCH DO NOT LIE IN A PLANE" );
		  stop(-1);
		}

	  } /* if( ntp != 4) */

	} /* if( ntp <= 2) */

  } /* if( ntp <= 1) */

  data.t2x[mi]= ynv* data.t1z[mi]- znv* data.t1y[mi];
  data.t2y[mi]= znv* data.t1x[mi]- xnv* data.t1z[mi];
  data.t2z[mi]= xnv* data.t1y[mi]- ynv* data.t1x[mi];
  data.psalp[mi]=1.;

  if( nx != 0)
  {
	data.m += nx*ny-1;

	/* Reallocate patch buffers */
	mreq = data.m * sizeof(long double);
	mem_realloc( (void *)&data.px, mreq );
	mem_realloc( (void *)&data.py, mreq );
	mem_realloc( (void *)&data.pz, mreq );
	mem_realloc( (void *)&data.t1x, mreq );
	mem_realloc( (void *)&data.t1y, mreq );
	mem_realloc( (void *)&data.t1z, mreq );
	mem_realloc( (void *)&data.t2x, mreq );
	mem_realloc( (void *)&data.t2y, mreq );
	mem_realloc( (void *)&data.t2z, mreq );
	mem_realloc( (void *)&data.pbi, mreq );
	mem_realloc( (void *)&data.psalp, mreq );

	xn2= data.px[mi]- s1x- s2x;
	yn2= data.py[mi]- s1y- s2y;
	zn2= data.pz[mi]- s1z- s2z;
	xs= data.t1x[mi];
	ys= data.t1y[mi];
	zs= data.t1z[mi];
	xt= data.t2x[mi];
	yt= data.t2y[mi];
	zt= data.t2z[mi];

	for( iy = 0; iy < ny; iy++ )
	{
	  xn2 += s2x;
	  yn2 += s2y;
	  zn2 += s2z;

	  for( ix = 1; ix <= nx; ix++ )
	  {
		xst= (long double)ix;
		data.px[mi]= xn2+ xst* s1x;
		data.py[mi]= yn2+ xst* s1y;
		data.pz[mi]= zn2+ xst* s1z;
		data.pbi[mi]= xa;
		data.psalp[mi]=1.;
		data.t1x[mi]= xs;
		data.t1y[mi]= ys;
		data.t1z[mi]= zs;
		data.t2x[mi]= xt;
		data.t2y[mi]= yt;
		data.t2z[mi]= zt;
		mi++;
	  } /* for( ix = 0; ix < nx; ix++ ) */

	} /* for( iy = 0; iy < ny; iy++ ) */

  } /* if( nx != 0) */

  data.ipsym=0;
  data.np= data.n;
  data.mp= data.m;

  return;
}

/*-----------------------------------------------------------------------*/

/*** this function was an 'entry point' (part of) 'patch()' ***/
void subph( int nx, int ny )
{
  int mia, ix, iy, mi, mreq;
  long double xs, ys, zs, xa, xst, s1x, s1y, s1z, s2x, s2y, s2z, saln, xt, yt;

  /* Reallocate patch buffers */
  if( ny == 0 )
	data.m += 3;
  else
	data.m += 4;

  mreq = data.m * sizeof(long double);
  mem_realloc( (void *)&data.px, mreq );
  mem_realloc( (void *)&data.py, mreq );
  mem_realloc( (void *)&data.pz, mreq );
  mem_realloc( (void *)&data.t1x, mreq );
  mem_realloc( (void *)&data.t1y, mreq );
  mem_realloc( (void *)&data.t1z, mreq );
  mem_realloc( (void *)&data.t2x, mreq );
  mem_realloc( (void *)&data.t2y, mreq );
  mem_realloc( (void *)&data.t2z, mreq );
  mem_realloc( (void *)&data.pbi, mreq );
  mem_realloc( (void *)&data.psalp, mreq );
  mem_realloc( (void *)&data.icon1, (data.n+data.m) * sizeof(int) );
  mem_realloc( (void *)&data.icon2, (data.n+data.m) * sizeof(int) );


  /* Shift patches to make room for new ones */
  if( (ny == 0) && (nx != data.m) )
  {
	for( iy = data.m-1; iy > nx+2; iy-- )
	{
	  ix = iy-3;
	  data.px[iy]= data.px[ix];
	  data.py[iy]= data.py[ix];
	  data.pz[iy]= data.pz[ix];
	  data.pbi[iy]= data.pbi[ix];
	  data.psalp[iy]= data.psalp[ix];
	  data.t1x[iy]= data.t1x[ix];
	  data.t1y[iy]= data.t1y[ix];
	  data.t1z[iy]= data.t1z[ix];
	  data.t2x[iy]= data.t2x[ix];
	  data.t2y[iy]= data.t2y[ix];
	  data.t2z[iy]= data.t2z[ix];
	}

  } /* if( (ny == 0) || (nx != m) ) */

  /* divide patch for connection */
  mi= nx-1;
  xs= data.px[mi];
  ys= data.py[mi];
  zs= data.pz[mi];
  xa= data.pbi[mi]/4.;
  xst= sqrtl( xa)/2.;
  s1x= data.t1x[mi];
  s1y= data.t1y[mi];
  s1z= data.t1z[mi];
  s2x= data.t2x[mi];
  s2y= data.t2y[mi];
  s2z= data.t2z[mi];
  saln= data.psalp[mi];
  xt= xst;
  yt= xst;

  if( ny == 0)
	mia= mi;
  else
  {
	data.mp++;
	mia= data.m-1;
  }

  for( ix = 1; ix <= 4; ix++ )
  {
	data.px[mia]= xs+ xt* s1x+ yt* s2x;
	data.py[mia]= ys+ xt* s1y+ yt* s2y;
	data.pz[mia]= zs+ xt* s1z+ yt* s2z;
	data.pbi[mia]= xa;
	data.t1x[mia]= s1x;
	data.t1y[mia]= s1y;
	data.t1z[mia]= s1z;
	data.t2x[mia]= s2x;
	data.t2y[mia]= s2y;
	data.t2z[mia]= s2z;
	data.psalp[mia]= saln;

	if( ix == 2)
	  yt=- yt;

	if( (ix == 1) || (ix == 3) )
	  xt=- xt;

	mia++;
  }

  if( nx <= data.mp)
	data.mp += 3;

  if( ny > 0 )
	data.pz[mi]=10000.;

  return;
}

/*-----------------------------------------------------------------------*/

void readgm( char *gm, int *i1, int *i2, long double *x1, long double *y1,
	long double *z1, long double *x2, long double *y2, long double *z2, long double *rad )
{
  char line_buf[134];
  int nlin, i, line_idx;
  int nint = 2, nflt = 7;
  int iarr[2] = { 0, 0 };
  long double rarr[7] = { 0., 0., 0., 0., 0., 0., 0. };


  /* read a line from input file */
  load_line( line_buf, input_fp );

  /* get line length */
  nlin= strlen( line_buf );

  /* abort if card's mnemonic too short or missing */
  if( nlin < 2 )
  {
	fprintf( output_fp,
		"\n  GEOMETRY DATA CARD ERROR:"
		"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
	stop(-1);
  }

  /* extract card's mnemonic code */
  strncpy( gm, line_buf, 2 );
  gm[2] = '\0';

  /* Exit if "XT" command read (for testing) */
  if( strcmp( gm, "XT" ) == 0 )
  {
	fprintf( stderr,
		"\nnec2c: Exiting after an \"XT\" command in readgm()\n" );
	fprintf( output_fp,
		"\n\n  nec2c: Exiting after an \"XT\" command in readgm()" );
	stop(0);
  }

  /* Return if only mnemonic on card */
  if( nlin == 2 )
  {
	*i1 = *i2 = 0;
	*x1 = *y1 = *z1 = *x2 = *y2 = *z2 = *rad = 0.;
	return;
  }

  /* read integers from line */
  line_idx = 1;
  for( i = 0; i < nint; i++ )
  {
	/* Find first numerical character */
	while( ((line_buf[++line_idx] <  '0')  ||
		  (line_buf[  line_idx] >  '9')) &&
		(line_buf[  line_idx] != '+')  &&
		(line_buf[  line_idx] != '-') )
	  if( (line_buf[line_idx] == '\0') )
	  {
		*i1= iarr[0];
		*i2= iarr[1];
		*x1= rarr[0];
		*y1= rarr[1];
		*z1= rarr[2];
		*x2= rarr[3];
		*y2= rarr[4];
		*z2= rarr[5];
		*rad= rarr[6];
		return;
	  }

	/* read an integer from line */
	iarr[i] = atoi( &line_buf[line_idx] );

	/* traverse numerical field to next ' ' or ',' or '\0' */
	line_idx--;
	while(
		(line_buf[++line_idx] != ' ') &&
		(line_buf[  line_idx] != '	') &&
		(line_buf[  line_idx] != ',') &&
		(line_buf[  line_idx] != '\0') )
	{
	  /* test for non-numerical characters */
	  if( ((line_buf[line_idx] <  '0')  ||
			(line_buf[line_idx] >  '9')) &&
		  (line_buf[line_idx] != '+')  &&
		  (line_buf[line_idx] != '-') )
	  {
		fprintf( output_fp,
			"\n  GEOMETRY DATA CARD \"%s\" ERROR:"
			"\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
			gm, line_buf[line_idx], (line_idx+1)  );
		stop(-1);
	  }

	} /* while( (line_buff[++line_idx] ... */

	/* Return on end of line */
	if( line_buf[line_idx] == '\0' )
	{
	  *i1= iarr[0];
	  *i2= iarr[1];
	  *x1= rarr[0];
	  *y1= rarr[1];
	  *z1= rarr[2];
	  *x2= rarr[3];
	  *y2= rarr[4];
	  *z2= rarr[5];
	  *rad= rarr[6];
	  return;
	}

  } /* for( i = 0; i < nint; i++ ) */

  /* read long doubles from line */
  for( i = 0; i < nflt; i++ )
  {
	/* Find first numerical character */
	while( ((line_buf[++line_idx] <  '0')  ||
		  (line_buf[  line_idx] >  '9')) &&
		(line_buf[  line_idx] != '+')  &&
		(line_buf[  line_idx] != '-')  &&
		(line_buf[  line_idx] != '.') )
	  if( (line_buf[line_idx] == '\0') )
	  {
		*i1= iarr[0];
		*i2= iarr[1];
		*x1= rarr[0];
		*y1= rarr[1];
		*z1= rarr[2];
		*x2= rarr[3];
		*y2= rarr[4];
		*z2= rarr[5];
		*rad= rarr[6];
		return;
	  }

	/* read a long double from line */
	rarr[i] = atof( &line_buf[line_idx] );

	/* traverse numerical field to next ' ' or ',' or '\0' */
	line_idx--;
	while(
		(line_buf[++line_idx] != ' ')  &&
		(line_buf[  line_idx] != '	') &&
		(line_buf[  line_idx] != ',')  &&
		(line_buf[  line_idx] != '\0') )
	{
	  /* test for non-numerical characters */
	  if( ((line_buf[line_idx] <  '0')  ||
			(line_buf[line_idx] >  '9')) &&
		  (line_buf[line_idx] != '.')  &&
		  (line_buf[line_idx] != '+')  &&
		  (line_buf[line_idx] != '-')  &&
		  (line_buf[line_idx] != 'E')  &&
		  (line_buf[line_idx] != 'e') )
	  {
		fprintf( output_fp,
			"\n  GEOMETRY DATA CARD \"%s\" ERROR:"
			"\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d.\n",
			gm, line_buf[line_idx], (line_idx+1) );
		stop(-1);
	  }

	} /* while( (line_buff[++line_idx] ... */

	/* Return on end of line */
	if( line_buf[line_idx] == '\0' )
	{
	  *i1= iarr[0];
	  *i2= iarr[1];
	  *x1= rarr[0];
	  *y1= rarr[1];
	  *z1= rarr[2];
	  *x2= rarr[3];
	  *y2= rarr[4];
	  *z2= rarr[5];
	  *rad= rarr[6];
	  return;
	}

  } /* for( i = 0; i < nflt; i++ ) */

  *i1  = iarr[0];
  *i2  = iarr[1];
  *x1  = rarr[0];
  *y1  = rarr[1];
  *z1  = rarr[2];
  *x2  = rarr[3];
  *y2  = rarr[4];
  *z2  = rarr[5];
  *rad = rarr[6];

  return;
}

/*-----------------------------------------------------------------------*/

/* reflc reflects partial structure along x,y, or z axes or rotates */
/* structure to complete a symmetric structure. */
void reflc( int ix, int iy, int iz, int itx, int nop )
{
  int iti, i, nx, itagi, k, mreq;
  long double e1, e2, fnop, sam, cs, ss, xk, yk;

  data.np= data.n;
  data.mp= data.m;
  data.ipsym=0;
  iti= itx;

  if( ix >= 0)
  {
	if( nop == 0)
	  return;

	data.ipsym=1;

	/* reflect along z axis */
	if( iz != 0)
	{
	  data.ipsym=2;

	  if( data.n > 0 )
	  {
		/* Reallocate tags buffer */
		mem_realloc( (void *)&data.itag, (2*data.n+data.m) * sizeof(int) );

		/* Reallocate wire buffers */
		mreq = 2*data.n * sizeof(long double);
		mem_realloc( (void *)&data.x1, mreq );
		mem_realloc( (void *)&data.y1, mreq );
		mem_realloc( (void *)&data.z1, mreq );
		mem_realloc( (void *)&data.x2, mreq );
		mem_realloc( (void *)&data.y2, mreq );
		mem_realloc( (void *)&data.z2, mreq );
		mem_realloc( (void *)&data.bi, mreq );

		for( i = 0; i < data.n; i++ )
		{
		  nx= i+ data.n;
		  e1= data.z1[i];
		  e2= data.z2[i];

		  if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
		  {
			fprintf( output_fp,
				"\n  GEOMETRY DATA ERROR--SEGMENT %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
			stop(-1);
		  }

		  data.x1[nx]= data.x1[i];
		  data.y1[nx]= data.y1[i];
		  data.z1[nx]=- e1;
		  data.x2[nx]= data.x2[i];
		  data.y2[nx]= data.y2[i];
		  data.z2[nx]=- e2;
		  itagi= data.itag[i];

		  if( itagi == 0)
			data.itag[nx]=0;
		  if( itagi != 0)
			data.itag[nx]= itagi+ iti;

		  data.bi[nx]= data.bi[i];

		} /* for( i = 0; i < data.n; i++ ) */

		data.n= data.n*2;
		iti= iti*2;

	  } /* if( data.n > 0) */

	  if( data.m > 0 )
	  {
		/* Reallocate patch buffers */
		mreq = 2*data.m * sizeof(long double);
		mem_realloc( (void *)&data.px, mreq );
		mem_realloc( (void *)&data.py, mreq );
		mem_realloc( (void *)&data.pz, mreq );
		mem_realloc( (void *)&data.t1x, mreq );
		mem_realloc( (void *)&data.t1y, mreq );
		mem_realloc( (void *)&data.t1z, mreq );
		mem_realloc( (void *)&data.t2x, mreq );
		mem_realloc( (void *)&data.t2y, mreq );
		mem_realloc( (void *)&data.t2z, mreq );
		mem_realloc( (void *)&data.pbi, mreq );
		mem_realloc( (void *)&data.psalp, mreq );

		for( i = 0; i < data.m; i++ )
		{
		  nx = i+data.m;
		  if( fabsl(data.pz[i]) <= 1.0e-10)
		  {
			fprintf( output_fp,
				"\n  GEOMETRY DATA ERROR--PATCH %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
			stop(-1);
		  }

		  data.px[nx]= data.px[i];
		  data.py[nx]= data.py[i];
		  data.pz[nx]=- data.pz[i];
		  data.t1x[nx]= data.t1x[i];
		  data.t1y[nx]= data.t1y[i];
		  data.t1z[nx]=- data.t1z[i];
		  data.t2x[nx]= data.t2x[i];
		  data.t2y[nx]= data.t2y[i];
		  data.t2z[nx]=- data.t2z[i];
		  data.psalp[nx]=- data.psalp[i];
		  data.pbi[nx]= data.pbi[i];
		}

		data.m= data.m*2;

	  } /* if( data.m >= m2) */

	} /* if( iz != 0) */

	/* reflect along y axis */
	if( iy != 0)
	{
	  if( data.n > 0)
	  {
		/* Reallocate tags buffer */
		mem_realloc( (void *)&data.itag, (2*data.n+data.m) * sizeof(int) );/*????*/

		/* Reallocate wire buffers */
		mreq = 2*data.n * sizeof(long double);
		mem_realloc( (void *)&data.x1, mreq );
		mem_realloc( (void *)&data.y1, mreq );
		mem_realloc( (void *)&data.z1, mreq );
		mem_realloc( (void *)&data.x2, mreq );
		mem_realloc( (void *)&data.y2, mreq );
		mem_realloc( (void *)&data.z2, mreq );
		mem_realloc( (void *)&data.bi, mreq );

		for( i = 0; i < data.n; i++ )
		{
		  nx= i+ data.n;
		  e1= data.y1[i];
		  e2= data.y2[i];

		  if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
		  {
			fprintf( output_fp,
				"\n  GEOMETRY DATA ERROR--SEGMENT %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
			stop(-1);
		  }

		  data.x1[nx]= data.x1[i];
		  data.y1[nx]=- e1;
		  data.z1[nx]= data.z1[i];
		  data.x2[nx]= data.x2[i];
		  data.y2[nx]=- e2;
		  data.z2[nx]= data.z2[i];
		  itagi= data.itag[i];

		  if( itagi == 0)
			data.itag[nx]=0;
		  if( itagi != 0)
			data.itag[nx]= itagi+ iti;

		  data.bi[nx]= data.bi[i];

		} /* for( i = n2-1; i < data.n; i++ ) */

		data.n= data.n*2;
		iti= iti*2;

	  } /* if( data.n >= n2) */

	  if( data.m > 0 )
	  {
		/* Reallocate patch buffers */
		mreq = 2*data.m * sizeof(long double);
		mem_realloc( (void *)&data.px, mreq );
		mem_realloc( (void *)&data.py, mreq );
		mem_realloc( (void *)&data.pz, mreq );
		mem_realloc( (void *)&data.t1x, mreq );
		mem_realloc( (void *)&data.t1y, mreq );
		mem_realloc( (void *)&data.t1z, mreq );
		mem_realloc( (void *)&data.t2x, mreq );
		mem_realloc( (void *)&data.t2y, mreq );
		mem_realloc( (void *)&data.t2z, mreq );
		mem_realloc( (void *)&data.pbi, mreq );
		mem_realloc( (void *)&data.psalp, mreq );

		for( i = 0; i < data.m; i++ )
		{
		  nx= i+data.m;
		  if( fabsl( data.py[i]) <= 1.0e-10)
		  {
			fprintf( output_fp,
				"\n  GEOMETRY DATA ERROR--PATCH %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
			stop(-1);
		  }

		  data.px[nx]= data.px[i];
		  data.py[nx]=- data.py[i];
		  data.pz[nx]= data.pz[i];
		  data.t1x[nx]= data.t1x[i];
		  data.t1y[nx]=- data.t1y[i];
		  data.t1z[nx]= data.t1z[i];
		  data.t2x[nx]= data.t2x[i];
		  data.t2y[nx]=- data.t2y[i];
		  data.t2z[nx]= data.t2z[i];
		  data.psalp[nx]=- data.psalp[i];
		  data.pbi[nx]= data.pbi[i];

		} /* for( i = m2; i <= data.m; i++ ) */

		data.m= data.m*2;

	  } /* if( data.m >= m2) */

	} /* if( iy != 0) */

	/* reflect along x axis */
	if( ix == 0 )
	  return;

	if( data.n > 0 )
	{
	  /* Reallocate tags buffer */
	  mem_realloc( (void *)&data.itag, (2*data.n+data.m) * sizeof(int) );/*????*/

	  /* Reallocate wire buffers */
	  mreq = 2*data.n * sizeof(long double);
	  mem_realloc( (void *)&data.x1, mreq );
	  mem_realloc( (void *)&data.y1, mreq );
	  mem_realloc( (void *)&data.z1, mreq );
	  mem_realloc( (void *)&data.x2, mreq );
	  mem_realloc( (void *)&data.y2, mreq );
	  mem_realloc( (void *)&data.z2, mreq );
	  mem_realloc( (void *)&data.bi, mreq );

	  for( i = 0; i < data.n; i++ )
	  {
		nx= i+ data.n;
		e1= data.x1[i];
		e2= data.x2[i];

		if( (fabsl(e1)+fabsl(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
		{
		  fprintf( output_fp,
			  "\n  GEOMETRY DATA ERROR--SEGMENT %d"
			  " LIES IN PLANE OF SYMMETRY", i+1 );
		  stop(-1);
		}

		data.x1[nx]=- e1;
		data.y1[nx]= data.y1[i];
		data.z1[nx]= data.z1[i];
		data.x2[nx]=- e2;
		data.y2[nx]= data.y2[i];
		data.z2[nx]= data.z2[i];
		itagi= data.itag[i];

		if( itagi == 0)
		  data.itag[nx]=0;
		if( itagi != 0)
		  data.itag[nx]= itagi+ iti;

		data.bi[nx]= data.bi[i];
	  }

	  data.n= data.n*2;

	} /* if( data.n > 0) */

	if( data.m == 0 )
	  return;

	/* Reallocate patch buffers */
	mreq = 2*data.m * sizeof(long double);
	mem_realloc( (void *)&data.px, mreq );
	mem_realloc( (void *)&data.py, mreq );
	mem_realloc( (void *)&data.pz, mreq );
	mem_realloc( (void *)&data.t1x, mreq );
	mem_realloc( (void *)&data.t1y, mreq );
	mem_realloc( (void *)&data.t1z, mreq );
	mem_realloc( (void *)&data.t2x, mreq );
	mem_realloc( (void *)&data.t2y, mreq );
	mem_realloc( (void *)&data.t2z, mreq );
	mem_realloc( (void *)&data.pbi, mreq );
	mem_realloc( (void *)&data.psalp, mreq );

	for( i = 0; i < data.m; i++ )
	{
	  nx= i+data.m;
	  if( fabsl( data.px[i]) <= 1.0e-10)
	  {
		fprintf( output_fp,
			"\n  GEOMETRY DATA ERROR--PATCH %d"
			" LIES IN PLANE OF SYMMETRY", i+1 );
		stop(-1);
	  }

	  data.px[nx]=- data.px[i];
	  data.py[nx]= data.py[i];
	  data.pz[nx]= data.pz[i];
	  data.t1x[nx]=- data.t1x[i];
	  data.t1y[nx]= data.t1y[i];
	  data.t1z[nx]= data.t1z[i];
	  data.t2x[nx]=- data.t2x[i];
	  data.t2y[nx]= data.t2y[i];
	  data.t2z[nx]= data.t2z[i];
	  data.psalp[nx]=- data.psalp[i];
	  data.pbi[nx]= data.pbi[i];
	}

	data.m= data.m*2;
	return;

  } /* if( ix >= 0) */

  /* reproduce structure with rotation to form cylindrical structure */
  fnop= (long double)nop;
  data.ipsym=-1;
  sam=TP/ fnop;
  cs= cosl( sam);
  ss= sinl( sam);

  if( data.n > 0)
  {
	data.n *= nop;
	nx= data.np;

	/* Reallocate tags buffer */
	mem_realloc( (void *)&data.itag, (data.n+data.m) * sizeof(int) );/*????*/

	/* Reallocate wire buffers */
	mreq = data.n * sizeof(long double);
	mem_realloc( (void *)&data.x1, mreq );
	mem_realloc( (void *)&data.y1, mreq );
	mem_realloc( (void *)&data.z1, mreq );
	mem_realloc( (void *)&data.x2, mreq );
	mem_realloc( (void *)&data.y2, mreq );
	mem_realloc( (void *)&data.z2, mreq );
	mem_realloc( (void *)&data.bi, mreq );

	for( i = nx; i < data.n; i++ )
	{
	  k= i- data.np;
	  xk= data.x1[k];
	  yk= data.y1[k];
	  data.x1[i]= xk* cs- yk* ss;
	  data.y1[i]= xk* ss+ yk* cs;
	  data.z1[i]= data.z1[k];
	  xk= data.x2[k];
	  yk= data.y2[k];
	  data.x2[i]= xk* cs- yk* ss;
	  data.y2[i]= xk* ss+ yk* cs;
	  data.z2[i]= data.z2[k];
	  data.bi[i]= data.bi[k];
	  itagi= data.itag[k];

	  if( itagi == 0)
		data.itag[i]=0;
	  if( itagi != 0)
		data.itag[i]= itagi+ iti;
	}

  } /* if( data.n >= n2) */

  if( data.m == 0 )
	return;

  data.m *= nop;
  nx= data.mp;

  /* Reallocate patch buffers */
  mreq = data.m * sizeof(long double);
  mem_realloc( (void *)&data.px, mreq  );
  mem_realloc( (void *)&data.py, mreq  );
  mem_realloc( (void *)&data.pz, mreq );
  mem_realloc( (void *)&data.t1x, mreq );
  mem_realloc( (void *)&data.t1y, mreq );
  mem_realloc( (void *)&data.t1z, mreq );
  mem_realloc( (void *)&data.t2x, mreq );
  mem_realloc( (void *)&data.t2y, mreq );
  mem_realloc( (void *)&data.t2z, mreq );
  mem_realloc( (void *)&data.pbi, mreq );
  mem_realloc( (void *)&data.psalp, mreq );

  for( i = nx; i < data.m; i++ )
  {
	k = i-data.mp;
	xk= data.px[k];
	yk= data.py[k];
	data.px[i]= xk* cs- yk* ss;
	data.py[i]= xk* ss+ yk* cs;
	data.pz[i]= data.pz[k];
	xk= data.t1x[k];
	yk= data.t1y[k];
	data.t1x[i]= xk* cs- yk* ss;
	data.t1y[i]= xk* ss+ yk* cs;
	data.t1z[i]= data.t1z[k];
	xk= data.t2x[k];
	yk= data.t2y[k];
	data.t2x[i]= xk* cs- yk* ss;
	data.t2y[i]= xk* ss+ yk* cs;
	data.t2z[i]= data.t2z[k];
	data.psalp[i]= data.psalp[k];
	data.pbi[i]= data.pbi[k];

  } /* for( i = nx; i < data.m; i++ ) */

  return;
}

/*-----------------------------------------------------------------------*/

/* subroutine wire generates segment geometry */
/* data for a straight wire of ns segments. */
void wire( long double xw1, long double yw1, long double zw1,
	long double xw2, long double yw2, long double zw2, long double rad,
	long double rdel, long double rrad, int ns, int itg )
{
  int ist, i, mreq;
  long double xd, yd, zd, delz, rd, fns, radz;
  long double xs1, ys1, zs1, xs2, ys2, zs2;

  ist= data.n;
  data.n= data.n+ ns;
  data.np= data.n;
  data.mp= data.m;
  data.ipsym=0;

  if( ns < 1)
	return;

  /* Reallocate tags buffer */
  mem_realloc( (void *)&data.itag, (data.n+data.m) * sizeof(int) );/*????*/

  /* Reallocate wire buffers */
  mreq = data.n * sizeof(long double);
  mem_realloc( (void *)&data.x1, mreq );
  mem_realloc( (void *)&data.y1, mreq );
  mem_realloc( (void *)&data.z1, mreq );
  mem_realloc( (void *)&data.x2, mreq );
  mem_realloc( (void *)&data.y2, mreq );
  mem_realloc( (void *)&data.z2, mreq );
  mem_realloc( (void *)&data.bi, mreq );

  xd= xw2- xw1;
  yd= yw2- yw1;
  zd= zw2- zw1;

  if( fabsl( rdel-1.) >= 1.0e-6)
  {
	delz= sqrtl( xd* xd+ yd* yd+ zd* zd);
	xd= xd/ delz;
	yd= yd/ delz;
	zd= zd/ delz;
	delz= delz*(1.- rdel)/(1.- powl(rdel, ns) );
	rd= rdel;
  }
  else
  {
	fns= ns;
	xd= xd/ fns;
	yd= yd/ fns;
	zd= zd/ fns;
	delz=1.;
	rd=1.;
  }

  radz= rad;
  xs1= xw1;
  ys1= yw1;
  zs1= zw1;

  for( i = ist; i < data.n; i++ )
  {
	data.itag[i]= itg;
	xs2= xs1+ xd* delz;
	ys2= ys1+ yd* delz;
	zs2= zs1+ zd* delz;
	data.x1[i]= xs1;
	data.y1[i]= ys1;
	data.z1[i]= zs1;
	data.x2[i]= xs2;
	data.y2[i]= ys2;
	data.z2[i]= zs2;
	data.bi[i]= radz;
	delz= delz* rd;
	radz= radz* rrad;
	xs1= xs2;
	ys1= ys2;
	zs1= zs2;
  }

  data.x2[data.n-1]= xw2;
  data.y2[data.n-1]= yw2;
  data.z2[data.n-1]= zw2;

  return;
}

/*-----------------------------------------------------------------------*/

