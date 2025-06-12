//------------------------------------------------------------------------------
// LDL/Include/ldl.h: include file for the LDL package
//------------------------------------------------------------------------------

// LDL, Copyright (c) 2005-2024 by Timothy A. Davis. All Rights Reserved.
// SPDX-License-Identifier: LGPL-2.1+

//------------------------------------------------------------------------------

#ifndef LDL_H
#define LDL_H

// #include "SuiteSparse_config.h"

#ifdef LDL_LONG
#define LDL_int int64_t
#define LDL_ID  "%" PRId64

#define LDL_symbolic ldl_l_symbolic
#define LDL_numeric ldl_l_numeric
#define LDL_lsolve ldl_l_lsolve
#define LDL_dsolve ldl_l_dsolve
#define LDL_ltsolve ldl_l_ltsolve
#define LDL_perm ldl_l_perm
#define LDL_permt ldl_l_permt
#define LDL_valid_perm ldl_l_valid_perm
#define LDL_valid_matrix ldl_l_valid_matrix

#else
#ifdef MATLAB_COMP
    #define LDL_int long long int
#else
    #define LDL_int int
#endif
#define LDL_ID "%d"

#define LDL_symbolic ldl_symbolic
#define LDL_numeric ldl_numeric
#define LDL_lsolve ldl_lsolve
#define LDL_dsolve ldl_dsolve
#define LDL_ltsolve ldl_ltsolve
#define LDL_perm ldl_perm
#define LDL_permt ldl_permt
#define LDL_valid_perm ldl_valid_perm
#define LDL_valid_matrix ldl_valid_matrix

#endif

/*-------Custom overwrite to include compiling in single precision float format----------*/
#ifndef USE_SINGLE_PRECISION
    #define CSfloat double
#else
    #define CSfloat float
#endif

/* make it easy for C++ programs to include LDL */
#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* === LDL_int version ====================================================== */
/* ========================================================================== */

void ldl_symbolic (LDL_int n, LDL_int Ap [ ], LDL_int Ai [ ], LDL_int Lp [ ],
    LDL_int Parent [ ], LDL_int Lnz [ ], LDL_int Flag [ ], LDL_int P [ ],
    LDL_int Pinv [ ]) ;

LDL_int ldl_numeric (LDL_int n, LDL_int Ap [ ], LDL_int Ai [ ], CSfloat Ax [ ],
    LDL_int Lp [ ], LDL_int Parent [ ], LDL_int Lnz [ ], LDL_int Li [ ],
    CSfloat Lx [ ], CSfloat D [ ], CSfloat Y [ ], LDL_int Pattern [ ],
    LDL_int Flag [ ], LDL_int P [ ], LDL_int Pinv [ ]) ;

void ldl_lsolve (LDL_int n, CSfloat X [ ], LDL_int Lp [ ], LDL_int Li [ ],
    CSfloat Lx [ ]) ;

void ldl_dsolve (LDL_int n, CSfloat X [ ], CSfloat D [ ]) ;

void ldl_ltsolve (LDL_int n, CSfloat X [ ], LDL_int Lp [ ], LDL_int Li [ ],
    CSfloat Lx [ ]) ;

void ldl_perm  (LDL_int n, CSfloat X [ ], CSfloat B [ ], LDL_int P [ ]) ;
void ldl_permt (LDL_int n, CSfloat X [ ], CSfloat B [ ], LDL_int P [ ]) ;

LDL_int ldl_valid_perm (LDL_int n, LDL_int P [ ], LDL_int Flag [ ]) ;
LDL_int ldl_valid_matrix ( LDL_int n, LDL_int Ap [ ], LDL_int Ai [ ]) ;

#ifdef __cplusplus
}
#endif

#endif
