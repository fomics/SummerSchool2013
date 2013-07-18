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

   /*   ***  RETRIEVE THE DEFAULT SYSTEM CONTEXT FROM BLACS

SYNOPSIS
                     void Cblacs_get( int context, int request, int* value);

DESCRIPTION
       blacs_get reports information about the BLACS state. The request may be
       tied to a particular context or handle, or may request system wide
       information, depending upon the nature of the request. blacs_get can be
       used to return the system default context needed to construct a first
       BLACS grid and its associated context. In this implementation the
       default context corresponds to the MPI default context, of
       MPI_COMM_WORLD.

       The routines take the following arguments:

       context   Integer. (input)
                 The context/handle which is the target of the enquiry. For
                 context-independent requests, this is unused.

       request   Integer. (input)
                 Which BLACS parameters should be returned.

                 request = 0 : Report the default system context.

                 request = 1 : Report the BLACS message ID range.

                 request = 2 : Report the debug level, BLACS was  compiled
                 with.

                 request = 10 : Report the system context from which the
                 specified context was constructed from.

                 request = 11 : Report the number of rings, multiring topology
                 is presently using.

                 request = 12 : Report the number of branches, the general

    */

   Cblacs_get( 0, 0, *** );  /* Get the context (communicator) */

   /* ***   INITIALIZE THE VIRTUAL TOPOLOGY AS A 1D ARRAY OF SIZE 1 X NPE

SYNOPSIS
                     int Cblacs_gridinit( int* context, char * order, int
                     np_row, int np_col);


DESCRIPTION
       blacs_gridinit maps available processes onto a BLACS grid. The BLACS
       grid created is identified by a context handle that is subsequently
       used to identify that particular process grid among the many that may
       be generated. A BLACS grid must be created before using other BLACS
       calls, except where noted. blacs_gridinit is a collective or globally
       blocking operation amongst the participating processes.

       The routines take the following arguments:

       context   Integer. (input/output)
                 On input the handle of the system context to use (a default
                 system context can be obtained from blacs_get. On return, the
                 handle to the newly created grid.

       order     character*1. (input)
                 Specifies the nature of the mapping of processes onto the new
                 BLACS grid. A value that is not explicitly specified as
                 below, will default to row major ordering.

                 order = "R" : Use row major ordering

                 order = "C" : Use column major ordering

       np_row    Integer (input)
                 Specifies the number of rows in the process grid.

       np_col    Integer (input)
                 Specifies the number of columns in the process grid.

    */

   /* 1 x npe process grid with Row major ordering */
   Cblacs_gridinit( &context, ***, ***, *** ); 

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

   /* *** FACTORIZE THE MATRIX 

SYNOPSIS

  void pddttrf( int n, double *dl, double *d, double *du, int ja, int *desca, 
                double *af, int laf, double *work, int lwork, int info );

DESCRIPTION
       PDDTTRF computes a LU factorization of an N-by-N real tridiagonal
       diagonally dominant-like distributed matrix A(1:N, JA:JA+N-1).
       Reordering is used to increase parallelism in the factorization.  This
       reordering results in factors that are DIFFERENT from those produced by
       equivalent sequential codes. These factors cannot be used directly by
       users; however, they can be used in
       subsequent calls to PDDTTRS to solve linear systems.

       The factorization has the form

               P A(1:N, JA:JA+N-1) P^T = L U

       where U is a tridiagonal upper triangular matrix and L is tridiagonal
       lower triangular, and P is a permutation matrix.

   */

   pddttrf_( &n, dl, ***, ***, &ja, desca, ***, &laf, work, &lwork, &info );

   /* Solution */

   /* *** SOLVE THE SYSTEM USING THE ABOVE FACTORIZATION

SYNOPSIS

  void pddttrs( char trans, int n, int nrhs, double *dl, double *d, double *du, int ja, int *desca, double *b, int ib, 
		int *descb, double *af, int laf, double *work, int lwork, int info );

DESCRIPTION
       PDDTTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X =
       B(IB:IB+N-1, 1:NRHS)                                   or
                 A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)

       where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
       stored in A(1:N,JA:JA+N-1) and AF by PDDTTRF.
       A(1:N, JA:JA+N-1) is an N-by-N real
       tridiagonal diagonally dominant-like distributed
       matrix.

       Routine PDDTTRF MUST be called first.

   */

   pddttrs_( &trans, ***, &nrhs, ***, ***, ***, ***, desca, ***, &ib, ***, 
            af, &laf, work, &lwork, &info );

   printf( "MYPE=%i: x[:] = %7.4f %7.4f %7.4f %7.4f\n",
           mype, b[0], b[1], b[2], b[3]);

   /*    *** EXIT THE GRID

SYNOPSIS
                     void Cblacs_gridexit( int context);

DESCRIPTION
       blacs_gridexit frees resources associated with the specified context.
       After the call to blacs_gridexit the context is invalid. The numeric
       value of the context may be recycled if subsequent process grids are
       created with blacs_gridinit or blacs_gridmap.

       The routine takes the following arguments:

       context   Integer. (input)
                 The handle of the context to be freed.

   */

   Cblacs_gridexit( *** );
   Cblacs_exit( 0 );

}
