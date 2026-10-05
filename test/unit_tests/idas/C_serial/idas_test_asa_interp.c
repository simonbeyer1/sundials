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
 * steps than the interpolation order. The backward residual compares the
 * interpolated forward solution with exp(-t).
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "idas/idas.h"
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

static int dae_res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr,
                   void* user_data)
{
  N_VGetArrayPointer(rr)[0] = N_VGetArrayPointer(yp)[0] +
                              N_VGetArrayPointer(yy)[0];
  return 0;
}

/* Record the largest relative error of the interpolated forward solution */
static int dae_resB(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector yyB,
                    N_Vector ypB, N_Vector rrB, void* user_data)
{
  sunrealtype* max_err = (sunrealtype*)user_data;
  sunrealtype y_true   = SUNRexp(-t);
  sunrealtype err      = SUNRabs(N_VGetArrayPointer(yy)[0] - y_true) / y_true;

  if (err > *max_err) { *max_err = err; }
  N_VGetArrayPointer(rrB)[0] = N_VGetArrayPointer(ypB)[0] -
                               N_VGetArrayPointer(yyB)[0];
  return 0;
}

/* Returns the largest relative interpolation error or a negative value on
   failure */
static sunrealtype run_test(sunrealtype tf, long int steps, int interp)
{
  SUNContext sunctx   = NULL;
  N_Vector y          = NULL;
  N_Vector yp         = NULL;
  N_Vector yB         = NULL;
  N_Vector ypB        = NULL;
  SUNMatrix A         = NULL;
  SUNMatrix AB        = NULL;
  SUNLinearSolver LS  = NULL;
  SUNLinearSolver LSB = NULL;
  void* ida_mem       = NULL;

  sunrealtype rtol    = SUNMAX(SUN_RCONST(1.0e-10), 1000 * SUN_UNIT_ROUNDOFF);
  sunrealtype atol    = SUNMAX(SUN_RCONST(1.0e-14), SUN_UNIT_ROUNDOFF);
  sunrealtype tret    = ZERO;
  sunrealtype max_err = ZERO;
  sunrealtype result  = -ONE;
  int ncheck          = 0;
  int which           = 0;
  int flag            = 0;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return -ONE; }

  y  = N_VNew_Serial(1, sunctx);
  yp = N_VNew_Serial(1, sunctx);
  if (!y || !yp) { goto cleanup; }
  N_VConst(ONE, y);
  N_VConst(-ONE, yp);

  ida_mem = IDACreate(sunctx);
  if (!ida_mem) { goto cleanup; }
  if (IDAInit(ida_mem, dae_res, ZERO, y, yp)) { goto cleanup; }
  if (IDASStolerances(ida_mem, rtol, atol)) { goto cleanup; }

  A  = SUNDenseMatrix(1, 1, sunctx);
  LS = SUNLinSol_Dense(y, A, sunctx);
  if (!A || !LS) { goto cleanup; }
  if (IDASetLinearSolver(ida_mem, LS, A)) { goto cleanup; }

  if (IDAAdjInit(ida_mem, steps, interp)) { goto cleanup; }

  flag = IDASolveF(ida_mem, tf, &tret, y, yp, IDA_NORMAL, &ncheck);
  if (flag < 0)
  {
    fprintf(stderr, "IDASolveF returned %i\n", flag);
    goto cleanup;
  }

  yB  = N_VNew_Serial(1, sunctx);
  ypB = N_VNew_Serial(1, sunctx);
  if (!yB || !ypB) { goto cleanup; }
  N_VConst(ONE, yB);
  N_VConst(ONE, ypB);

  if (IDACreateB(ida_mem, &which)) { goto cleanup; }
  if (IDAInitB(ida_mem, which, dae_resB, tf, yB, ypB)) { goto cleanup; }
  if (IDASStolerancesB(ida_mem, which, rtol, atol)) { goto cleanup; }
  if (IDASetUserDataB(ida_mem, which, &max_err)) { goto cleanup; }

  AB  = SUNDenseMatrix(1, 1, sunctx);
  LSB = SUNLinSol_Dense(yB, AB, sunctx);
  if (!AB || !LSB) { goto cleanup; }
  if (IDASetLinearSolverB(ida_mem, which, LSB, AB)) { goto cleanup; }

  flag = IDASolveB(ida_mem, ZERO, IDA_NORMAL);
  if (flag < 0)
  {
    fprintf(stderr, "IDASolveB returned %i\n", flag);
    goto cleanup;
  }

  result = max_err;

cleanup:
  N_VDestroy(y);
  N_VDestroy(yp);
  N_VDestroy(yB);
  N_VDestroy(ypB);
  SUNMatDestroy(A);
  SUNMatDestroy(AB);
  SUNLinSolFree(LS);
  SUNLinSolFree(LSB);
  IDAFree(&ida_mem);
  SUNContext_Free(&sunctx);
  return result;
}

int main(int argc, char* argv[])
{
  const sunrealtype tf[2]  = {SUN_RCONST(2.0), SUN_RCONST(-2.0)};
  const int interp[2]      = {IDA_HERMITE, IDA_POLYNOMIAL};
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
