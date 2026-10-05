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
 * frequency of y' = cos(w t) is raised after IDASolveF, so the recomputation
 * in IDASolveB needs more steps than were stored per check point. IDASolveB
 * must return IDA_FWD_FAIL instead of writing past the end of the
 * interpolation data.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "idas/idas.h"
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

static int dae_res(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector rr,
                   void* user_data)
{
  sunrealtype w             = *((sunrealtype*)user_data);
  N_VGetArrayPointer(rr)[0] = N_VGetArrayPointer(yp)[0] - SUNRcos(w * t);
  return 0;
}

static int dae_resB(sunrealtype t, N_Vector yy, N_Vector yp, N_Vector yyB,
                    N_Vector ypB, N_Vector rrB, void* user_data)
{
  N_VGetArrayPointer(rrB)[0] = N_VGetArrayPointer(ypB)[0];
  return 0;
}

static int run_test(int interp)
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

  sunrealtype w    = ONE;
  sunrealtype tf   = SUN_RCONST(10.0);
  sunrealtype tret = ZERO;
  long int steps   = 10;
  int ncheck       = 0;
  int which        = 0;
  int flag         = 0;
  int result       = 1;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx)) { return 1; }

  y  = N_VNew_Serial(1, sunctx);
  yp = N_VNew_Serial(1, sunctx);
  if (!y || !yp) { goto cleanup; }
  N_VConst(ZERO, y);
  N_VConst(ONE, yp);

  ida_mem = IDACreate(sunctx);
  if (!ida_mem) { goto cleanup; }
  if (IDAInit(ida_mem, dae_res, ZERO, y, yp)) { goto cleanup; }
  if (IDASetUserData(ida_mem, &w)) { goto cleanup; }
  if (IDASStolerances(ida_mem, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }

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

  /* The recomputation from a check point now takes more steps */
  w = SUN_RCONST(20.0);

  yB  = N_VNew_Serial(1, sunctx);
  ypB = N_VNew_Serial(1, sunctx);
  if (!yB || !ypB) { goto cleanup; }
  N_VConst(ONE, yB);
  N_VConst(ZERO, ypB);

  if (IDACreateB(ida_mem, &which)) { goto cleanup; }
  if (IDAInitB(ida_mem, which, dae_resB, tf, yB, ypB)) { goto cleanup; }
  if (IDASStolerancesB(ida_mem, which, SUN_RCONST(1.0e-6), SUN_RCONST(1.0e-8)))
  {
    goto cleanup;
  }

  AB  = SUNDenseMatrix(1, 1, sunctx);
  LSB = SUNLinSol_Dense(yB, AB, sunctx);
  if (!AB || !LSB) { goto cleanup; }
  if (IDASetLinearSolverB(ida_mem, which, LSB, AB)) { goto cleanup; }

  flag = IDASolveB(ida_mem, ZERO, IDA_NORMAL);
  printf("interp = %i, check points = %i, IDASolveB returned %i\n", interp,
         ncheck, flag);

  if (flag == IDA_FWD_FAIL) { result = 0; }
  else { printf("ERROR: Expected IDA_FWD_FAIL (%i)\n", IDA_FWD_FAIL); }

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
  int fails = 0;

  fails += run_test(IDA_HERMITE);
  fails += run_test(IDA_POLYNOMIAL);

  if (fails) { printf("FAIL: %i test(s) failed\n", fails); }
  else { printf("SUCCESS\n"); }

  return fails;
}
