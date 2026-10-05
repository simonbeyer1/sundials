/* -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2026, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Unit test for adjoint sensitivity analysis when recomputing the forward
 * solution from a check point takes more steps than the first pass. The
 * frequency of y' = cos(w t) is raised after CVodeF, so the recomputation in
 * CVodeB needs more steps than were stored per check point. CVodeB must return
 * CV_FWD_FAIL instead of writing past the end of the interpolation data.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRcos(x) (cos((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRcos(x) (cosf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRcos(x) (cosl((x)))
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype w               = *((sunrealtype*)user_data);
  N_VGetArrayPointer(ydot)[0] = SUNRcos(w * t);
  return 0;
}

static int ode_rhsB(sunrealtype t, N_Vector y, N_Vector yB, N_Vector yBdot,
                    void* user_data)
{
  N_VGetArrayPointer(yBdot)[0] = ZERO;
  return 0;
}

static int run_test(int interp)
{
  SUNContext sunctx   = NULL;
  N_Vector y          = NULL;
  N_Vector yB         = NULL;
  SUNMatrix A         = NULL;
  SUNMatrix AB        = NULL;
  SUNLinearSolver LS  = NULL;
  SUNLinearSolver LSB = NULL;
  void* cvode_mem     = NULL;

  sunrealtype w    = ONE;
  sunrealtype tf   = SUN_RCONST(10.0);
  sunrealtype tret = ZERO;
  long int steps   = 10;
  int ncheck       = 0;
  int which        = 0;
  int flag         = 0;
  int result       = 1;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ZERO, y);

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { goto cleanup; }
  if (CVodeInit(cvode_mem, ode_rhs, ZERO, y)) { goto cleanup; }
  if (CVodeSetUserData(cvode_mem, &w)) { goto cleanup; }
  if (CVodeSStolerances(cvode_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }

  A  = SUNDenseMatrix(1, 1, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!A || !LS) { goto cleanup; }
  if (CVodeSetLinearSolver(cvode_mem, LS, A)) { goto cleanup; }

  if (CVodeAdjInit(cvode_mem, steps, interp)) { goto cleanup; }

  flag = CVodeF(cvode_mem, tf, y, &tret, CV_NORMAL, &ncheck);
  if (flag < 0)
  {
    fprintf(stderr, "CVodeF returned %i\n", flag);
    goto cleanup;
  }

  /* The recomputation from a check point now takes more steps */
  w = SUN_RCONST(20.0);

  yB = N_VNew_Serial(1, sunctx);
  if (!yB) { goto cleanup; }
  N_VConst(ONE, yB);

  if (CVodeCreateB(cvode_mem, CV_BDF, &which)) { goto cleanup; }
  if (CVodeInitB(cvode_mem, which, ode_rhsB, tf, yB)) { goto cleanup; }
  if (CVodeSStolerancesB(cvode_mem, which, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }

  AB  = SUNDenseMatrix(1, 1, sunctx);
  LSB = SUNLinSol_Dense(yB, AB, sunctx);
  if (!AB || !LSB) { goto cleanup; }
  if (CVodeSetLinearSolverB(cvode_mem, which, LSB, AB)) { goto cleanup; }

  flag = CVodeB(cvode_mem, ZERO, CV_NORMAL);
  printf("interp = %i, check points = %i, CVodeB returned %i\n", interp, ncheck,
         flag);

  if (flag == CV_FWD_FAIL) { result = 0; }
  else { printf("ERROR: Expected CV_FWD_FAIL (%i)\n", CV_FWD_FAIL); }

cleanup:
  N_VDestroy(y);
  N_VDestroy(yB);
  SUNMatDestroy(A);
  SUNMatDestroy(AB);
  SUNLinSolFree(LS);
  SUNLinSolFree(LSB);
  CVodeFree(&cvode_mem);
  SUNContext_Free(&sunctx);
  return result;
}

int main(int argc, char* argv[])
{
  int fails = 0;

  fails += run_test(CV_HERMITE);
  fails += run_test(CV_POLYNOMIAL);

  if (fails) { printf("FAIL: %i test(s) failed\n", fails); }
  else { printf("SUCCESS\n"); }

  return fails;
}
