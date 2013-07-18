/* =====================================================================

This program illustrates the use of the ScaLAPACK routines PDDTTRF
and PDDTTRS to factor and solve a (diagonally dominant) tridiagonal
system of linear equations T*x = b.          

Using Matlab notation, we define T = diag(D)+diag(DL,-1)+diag(DU,1)
and set 

D  = [ 1.8180 1.6602 1.3420 1.2897 1.3412 1.5341 1.7271 1.3093 ]
DL = [ 0.8385 0.5681 0.3704 0.7027 0.5466 0.4449 0.6946 ]
DU = [ 0.6946 0.4449 0.5466 0.7027 0.3704 0.5681 0.8385 ]

therefore, if

b = [ 1 2 3 4 5 6 7 8 ]

then

x = [ 0.2716 0.7289 1.2635 1.6287 2.0370 3.0330 0.4112 5.8920 ]

Important notes:
1. D, DL, DU and b have been hardcoded, which means that the code is 
   not malleable and works only for two processors.
2. DL(1) and DU(N) are not referenced (see the documentation of  
   PDDTTRF and PDDTTRS) and this has been reflected in the 
   distribution of DL and DU below.
3. PDDTTRS overwrites b with x

PE =  0, X(:) =  0.2716E+00  0.7289E+00  0.1264E+01  0.1629E+01
PE =  1, X(:) =  0.2037E+01  0.3033E+01  0.4112E+00  0.5892E+01

===================================================================== */

#include <stdio.h>
#include <math.h>
#include "mpi.h"

extern void Cblacs_get();
extern void Cblacs_pinfo();
extern void Cblacs_exit();
extern void Cblacs_gridexit();

main()
{

   int     context, desca[9], descb[9], ib, info, ja, laf, lda, ldb,
           lwork, mb, mype, n, nb, npcol, npe, nprow, nrhs;
   double  b[4], d[4], dl[4], du[4];
   double  *af, *work;
   char    trans='N';

   /* Set array dimensions and blocking */

   n = 8;            /* dimension of the problem */
   lda = 4;          /* leading dimension of A */
   ldb = 4;          /* leading dimension of B */
   nrhs = 1;         /* number of right-hand sides  */
   npcol = 2;        /* number of processor columns */
   ja = 1;           /* offset for A */
   ib = ja;          /* offset for B */
   mb = 4;           /* blocking */
   nb = 4;           /* blocking */

   laf = 12*npcol+3*nb;
   af = (double *)malloc(laf*sizeof(double));
   lwork = 10*npcol+4*nrhs;
   work = (double *)malloc(lwork*sizeof(double));

   /* Start BLACS */

   Cblacs_pinfo( &mype, &npe );
   Cblacs_get( 0, 0, &context );
   Cblacs_gridinit( &context, "R", 1, npe );

   if      ( mype == 0 ){ 
           /* PE = 0 gets D(1:4), DL(1:4), DU(1:4) and B(1:4) */
           d[0] = 1.8180; d[1] = 1.6602; d[2] = 1.3420; d[3] = 1.2897;
           dl[0] = 0.0000; dl[1] = 0.8385; dl[2] = 0.5681; dl[3] = 0.3704;
           du[0] = 0.6946; du[1] = 0.4449; du[2] = 0.5466; du[3] = 0.7027;
           b[0] = 1.0; b[1] = 2.0; b[2] = 3.0; b[3] = 4.0;
   }
   else if ( mype == 1 ){
           /* PE = 1 gets D(5:8), DL(5:8), DU(5:8) and B(5:8) */
           d[0] = 1.3412; d[1] = 1.5341; d[2] = 1.7271; d[3] = 1.3093;
           dl[0] = 0.7027; dl[1] = 0.5466; dl[2] = 0.4449; dl[3] = 0.6946;
           du[0] = 0.3704; du[1] = 0.5681; du[2] = 0.8385; du[3] = 0.0000;
           b[0] = 5.0; b[1] = 6.0; b[2] = 7.0; b[3] = 8.0;
   }

   /* Array descriptor for A (D, DL and DU) */

   desca[0] = 501; desca[1] = context; desca[2] = n; desca[3] = nb; 
   desca[4] = 0; desca[5] = lda; desca[6] = 0;

   /* Array descriptor for B */

   descb[0] = 502; descb[1] = context; descb[2] = n; descb[3] = nb; 
   descb[4] = 0; descb[5] = ldb; descb[6] = 0;

   /* Factorization */

   pddttrf_( &n, dl, d, du, &ja, desca, af, &laf, work, &lwork, &info );

   /* Solution */

   pddttrs_( &trans, &n, &nrhs, dl, d, du, &ja, desca, b, &ib, descb,
            af, &laf, work, &lwork, &info );

   printf( "MYPE=%i: x[:] = %7.4f %7.4f %7.4f %7.4f\n",
           mype, b[0], b[1], b[2], b[3]);

   Cblacs_gridexit( context );
   Cblacs_exit( 0 );

}
