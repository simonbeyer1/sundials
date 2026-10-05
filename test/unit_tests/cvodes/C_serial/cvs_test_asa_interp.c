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
 * Unit test for the interpolation of the forward solution in adjoint
 * sensitivity analysis. The forward problem y' = -y, y(0) = 1 is integrated to
 * tf = 2 and to tf = -2 with Hermite and polynomial interpolation and a range
 * of steps per check point, so that some check point intervals hold fewer
 * steps than the interpolation order. The backward right-hand side compares
 * the interpolated forward solution with exp(-t).
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cvodes/cvodes.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#else
#define GSYM "g"
#endif

#if defined(SUNDIALS_DOUBLE_PRECISION)
#define SUNRexp(x) (exp((x)))
#elif defined(SUNDIALS_SINGLE_PRECISION)
#define SUNRexp(x) (expf((x)))
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define SUNRexp(x) (expl((x)))
#endif

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* The polynomial interpolation reaches up to about 6e-3 in single precision,
   while the errors this test checks for are larger than 9e-2 */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define ERR_TOL SUN_RCONST(1.0e-2)
#else
#define ERR_TOL SUN_RCONST(1.0e-3)
#endif

static int ode_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  N_VGetArrayPointer(ydot)[0] = -N_VGetArrayPointer(y)[0];
  return 0;
}

/* Record the largest relative error of the interpolated forward solution */
static int ode_rhsB(sunrealtype t, N_Vector y, N_Vector yB, N_Vector yBdot,
                    void* user_data)
{
  sunrealtype* max_err = (sunrealtype*)user_data;
  sunrealtype y_true   = SUNRexp(-t);
  sunrealtype err      = SUNRabs(N_VGetArrayPointer(y)[0] - y_true) / y_true;

  if (err > *max_err) { *max_err = err; }
  N_VGetArrayPointer(yBdot)[0] = N_VGetArrayPointer(yB)[0];
  return 0;
}

/* Returns the largest relative interpolation error or a negative value on
   failure */
static sunrealtype run_test(sunrealtype tf, long int steps, int interp)
{
  SUNContext sunctx   = NULL;
  N_Vector y          = NULL;
  N_Vector yB         = NULL;
  SUNMatrix A         = NULL;
  SUNMatrix AB        = NULL;
  SUNLinearSolver LS  = NULL;
  SUNLinearSolver LSB = NULL;
  void* cvode_mem     = NULL;

  sunrealtype rtol    = SUNMAX(SUN_RCONST(1.0e-10), 1000 * SUN_UNIT_ROUNDOFF);
  sunrealtype atol    = SUNMAX(SUN_RCONST(1.0e-14), SUN_UNIT_ROUNDOFF);
  sunrealtype tret    = ZERO;
  sunrealtype max_err = ZERO;
  sunrealtype result  = -ONE;
  int ncheck          = 0;
  int which           = 0;
  int flag            = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return -ONE; }

  y = N_VNew_Serial(1, sunctx);
  if (!y) { goto cleanup; }
  N_VConst(ONE, y);

  cvode_mem = CVodeCreate(CV_BDF, sunctx);
  if (!cvode_mem) { goto cleanup; }
  if (CVodeInit(cvode_mem, ode_rhs, ZERO, y)) { goto cleanup; }
  if (CVodeSStolerances(cvode_mem, rtol, atol)) { goto cleanup; }

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

  yB = N_VNew_Serial(1, sunctx);
  if (!yB) { goto cleanup; }
  N_VConst(ONE, yB);

  if (CVodeCreateB(cvode_mem, CV_BDF, &which)) { goto cleanup; }
  if (CVodeInitB(cvode_mem, which, ode_rhsB, tf, yB)) { goto cleanup; }
  if (CVodeSStolerancesB(cvode_mem, which, rtol, atol)) { goto cleanup; }
  if (CVodeSetUserDataB(cvode_mem, which, &max_err)) { goto cleanup; }

  AB  = SUNDenseMatrix(1, 1, sunctx);
  LSB = SUNLinSol_Dense(yB, AB, sunctx);
  if (!AB || !LSB) { goto cleanup; }
  if (CVodeSetLinearSolverB(cvode_mem, which, LSB, AB)) { goto cleanup; }

  flag = CVodeB(cvode_mem, ZERO, CV_NORMAL);
  if (flag < 0)
  {
    fprintf(stderr, "CVodeB returned %i\n", flag);
    goto cleanup;
  }

  result = max_err;

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
  const sunrealtype tf[2]  = {SUN_RCONST(2.0), SUN_RCONST(-2.0)};
  const int interp[2]      = {CV_HERMITE, CV_POLYNOMIAL};
  const sunrealtype tol    = ERR_TOL;
  const long int max_steps = 40;

  sunrealtype err = ZERO;
  long int steps  = 0;
  int fails       = 0;
  int i, j;

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
    {
      for (steps = 2; steps <= max_steps + 1; steps++)
      {
        /* The last run uses a single check point interval */
        err = run_test(tf[i], steps > max_steps ? 1000 : steps, interp[j]);
        if (err < ZERO || err > tol)
        {
          printf("FAIL: tf = %" GSYM
                 ", interp = %i, steps = %ld, error = %" GSYM "\n",
                 tf[i], interp[j], steps > max_steps ? 1000 : steps, err);
          fails++;
        }
      }
    }
  }

  if (fails) { printf("FAIL: %i test(s) failed\n", fails); }
  else { printf("SUCCESS\n"); }

  return fails;
}
