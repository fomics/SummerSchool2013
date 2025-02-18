static char help[] = "Solves a linear system in parallel with KSP.  The matrix\n\
uses simple bilinear elements on the unit square.  To test the parallel\n\
matrix assembly, the matrix is intentionally laid out across processors\n\
differently from the way it is assembled.  Input arguments are:\n\
  -m <size> : problem size\n\n";

/*T
   Concepts: KSP^basic parallel example
   Concepts: Matrices^inserting elements by blocks
   Processors: n
T*/

/* 
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include <petscksp.h>

/* Declare user-defined routines */
extern PetscErrorCode FormElementStiffness(PetscReal,PetscScalar*);
extern PetscErrorCode FormElementRhs(PetscReal,PetscReal,PetscReal,PetscScalar*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            u,b,ustar; /* approx solution, RHS, exact solution */
  Mat            A;           /* linear system matrix */
  KSP            ksp;         /* Krylov subspace method context */
  PetscInt       N;           /* dimension of system (global) */
  PetscInt       M;           /* number of elements (global) */
  PetscMPIInt    rank;        /* processor rank */
  PetscMPIInt    size;        /* size of communicator */
  PetscScalar    Ke[16];      /* element matrix */
  PetscScalar    r[4];        /* element vector */
  PetscReal      h;           /* mesh width */
  PetscReal      norm;        /* norm of solution error */
  PetscReal      x,y;
  PetscScalar    val;
  PetscErrorCode ierr;
  PetscInt       idx[4],count,*rows,i,m = 5,start,end,its;

  PetscLogDouble tsetup,tsetup1,tsetup2,tsolve,tsolve1,tsolve2;
#if defined (PETSC_USE_LOG)
  PetscLogEvent  USER_EVENT;
#endif

  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
  N = (m+1)*(m+1);
  M = m*m;
  h = 1.0/m;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         Compute the matrix and right-hand-side vector that define
         the linear system, Au = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscGetTime(&tsetup1);

  /* 
     Create stiffness matrix of size NxN
     Careful! MatCreate expects a pointer to the object
  */
  ierr = MatCreate(PETSC_COMM_WORLD,***);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,***,***);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,9,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,9,PETSC_NULL,5,PETSC_NULL);CHKERRQ(ierr); /* More than necessary */
  start = rank*(M/size) + ((M%size) < rank ? (M%size) : rank);
  end   = start + M/size + ((M%size) > rank); 

  /*
     Assemble matrix
  */
  PetscLogEventRegister("Matrix Set",KSP_CLASSID,&USER_EVENT);
  ierr = FormElementStiffness(h*h,Ke);  /* Ke is 4x4 matrix to insert */
  for (i=start; i<end; i++) {
     /* location of lower left corner of element */
     x = h*(i % m); y = h*(i/m); 
     /* node numbers for the four corners of element */
     idx[0] = (m+1)*(i/m) + (i % m);
     idx[1] = idx[0]+1; idx[2] = idx[1] + m + 1; idx[3] = idx[2] - 1;
     ierr = MatSetValues(A,4,idx,***,***,***,ADD_VALUES);CHKERRQ(ierr);
  }
  /* Perform the (final) assembly of A */
  ierr = MatAssemblyBegin(***,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(***,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /*
     Create right-hand-side and solution vectors
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr); 
  ierr = VecSetSizes(u,PETSC_DECIDE,N);CHKERRQ(ierr); 
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)u,"Approx. Solution");CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)b,"Right hand side");CHKERRQ(ierr);
  ierr = VecDuplicate(b,&ustar);CHKERRQ(ierr);
  ierr = VecSet(u,0.0);CHKERRQ(ierr);
  ierr = VecSet(b,0.0);CHKERRQ(ierr);

  /* 
     Assemble right-hand-side vector
  */
  for (i=start; i<end; i++) {
     /* location of lower left corner of element */
     x = h*(i % m); y = h*(i/m); 
     /* node numbers for the four corners of element */
     idx[0] = (m+1)*(i/m) + (i % m);
     idx[1] = idx[0]+1; idx[2] = idx[1] + m + 1; idx[3] = idx[2] - 1;
     ierr = FormElementRhs(x,y,h*h,r);CHKERRQ(ierr);
     ierr = VecSetValues(b,4,idx,r,ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* 
     Modify matrix and right-hand-side for Dirichlet boundary conditions
  */
  ierr = PetscMalloc(4*m*sizeof(PetscInt),&rows);CHKERRQ(ierr);
  for (i=0; i<m+1; i++) {
    rows[i] = i; /* bottom */
    rows[3*m - 1 +i] = m*(m+1) + i; /* top */
  }
  count = m+1; /* left side */
  for (i=m+1; i<m*(m+1); i+= m+1) {
    rows[count++] = i;
  }
  count = 2*m; /* left side */
  for (i=2*m+1; i<m*(m+1); i+= m+1) {
    rows[count++] = i;
  }
  for (i=0; i<4*m; i++) {
     x = h*(rows[i] % (m+1)); y = h*(rows[i]/(m+1)); 
     val = y; /* value at boundary is linear interpolation in y direction*/
     ierr = VecSetValues(u,1,&rows[i],&val,INSERT_VALUES);CHKERRQ(ierr);
     /* insert one value (val) at the given row index */
     ierr = VecSetValues(b,***,***,***,INSERT_VALUES);CHKERRQ(ierr);
  }    
  ierr = MatZeroRows(A,4*m,rows,1.0,0,0);CHKERRQ(ierr);
  ierr = PetscFree(rows);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(u);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr); 
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = KSPCreate(PETSC_COMM_WORLD,***);CHKERRQ(ierr);  /* Create ksp operator, as in lecture */
  /* Bind the matrix A to the operator, use A to derive preconditioner */
  ierr = KSPSetOperators(***,***,***,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  PetscGetTime(&tsetup2);

  tsetup = tsetup2 - tsetup1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscGetTime(&tsolve1);
  ierr = KSPSolve(***,***,***);CHKERRQ(ierr); /* Solve: A u = b */
  PetscGetTime(&tsolve2);
  tsolve = tsolve2 - tsolve1;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* Check error */
  ierr = VecGetOwnershipRange(ustar,&start,&end);CHKERRQ(ierr);
  for (i=start; i<end; i++) {
     x = h*(i % (m+1)); y = h*(i/(m+1)); 
     val = y;
     ierr = VecSetValues(ustar,1,&i,&val,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = VecAssemblyBegin(ustar);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(ustar);CHKERRQ(ierr);
  ierr = VecAXPY(u,***,***);CHKERRQ(ierr); /* Subtract ustar from u */
  ierr = VecNorm(***,***,&norm);CHKERRQ(ierr); /* Find the 2-norm of u */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G Iterations %D\n",norm*h,its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "Time for setup %e Time for solve %e \n", tsetup, tsolve);CHKERRQ(ierr);

  /* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr); ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&ustar);CHKERRQ(ierr); ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();
  return 0;
}

/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormElementStiffness"
   /* element stiffness for Laplacian */
PetscErrorCode FormElementStiffness(PetscReal H,PetscScalar *Ke)
{
  PetscFunctionBegin;
  Ke[0]  = H/6.0;    Ke[1]  = -.125*H; Ke[2]  = H/12.0;   Ke[3]  = -.125*H;
  Ke[4]  = -.125*H;  Ke[5]  = H/6.0;   Ke[6]  = -.125*H;  Ke[7]  = H/12.0;
  Ke[8]  = H/12.0;   Ke[9]  = -.125*H; Ke[10] = H/6.0;    Ke[11] = -.125*H;
  Ke[12] = -.125*H;  Ke[13] = H/12.0;  Ke[14] = -.125*H;  Ke[15] = H/6.0;
  PetscFunctionReturn(0);
}
/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormElementRhs"
PetscErrorCode FormElementRhs(PetscReal x,PetscReal y,PetscReal H,PetscScalar *r)
{
  PetscFunctionBegin;
  r[0] = 0.; r[1] = 0.; r[2] = 0.; r[3] = 0.0; 
  PetscFunctionReturn(0);
}
