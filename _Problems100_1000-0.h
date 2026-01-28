/*=============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Simplex (No MPI)
Module: _Problems100_1000-0.h (LP problems of dimensions 100...1000 without random inequalities)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems were obtained using LPP-Generator https://github.com/leonid-sokolinsky/LPP-Generator
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Rnd-LP
===============================================================================*/
#pragma once

//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
#define PP_EPS_RELATIVE_ERROR			1E-11			// Used if defined PP_CHECK_MAX_OBJ_VALUE 

//============================== Problem Parameters ============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
//#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
//#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*1000)	// Accuracy of belonging to hyperplane
//==============================================================================

/*============================== cone100-1000-1 LP problem =========================*/
#define PP_PROBLEM_NAME	"cone100-1000-1"
#define PP_M	2100		// Number of equations (number of rows in *.mtx)
#define PP_N	2200		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 20000 
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	1E-10	// Accuracy of belonging to hyperplane
#define PP_MIN_COS				0.8		// Minimum allowable cosine of angle between launch vector and direction vector
#define PP_LAUNCH_VECTOR_LENGTH	1E+8	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 9.3181261 (MPI)
// Number of iterations: 919
// Computed objective value: 20000
// Maximal objective value:  20000
// Relative error = 0
// Distance to polytope: 2.1884716e-12
//------------------------------------------------------------------------------

/*============================== rnd100-0 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd100-0"
#define PP_KK	100		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	101		// Number of equations (number of rows in *.mtx)
#define PP_N	201		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 1009900 // =200*(n-1)*(2+n)/2+100
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------ ifdef _DEBUG --------------------------------
#define PP_ITER_COUNT			10				// Each PP_ITER_COUNT-th iteration to be outputted inside PC_bsf_MapF(*)
#define PP_PROJECTION_COUNT		1000000			// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)
//------------------------------------------------------------------------------
// Elapsed time: 13.555384
// Number of iterations: 5
// Computed objective value: 1009900
// Maximal objective value:  1009900
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== rnd100-10-1 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd100-10-1"
#define PP_KK	100		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	110		// Number of equations (number of rows in *.mtx)
#define PP_N	210		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 770948.875353609793819487 
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
#define PP_EPS_INVERSE			1E-10	// Accuracy for comparison with zero when calculating inverse matrix
#define PP_EPS_ON_HYPERPLANE	1E-10	// Accuracy of belonging to hyperplane
#define PP_MIN_COS				0.8		// Minimum allowable cosine of angle between launch vector and direction vector
#define PP_LAUNCH_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 64
// Computed objective value: 770948.875353609793819487
// Maximal objective value:  770948.875353609793819487
// Relative error = 0
// Distance to polytope: 1.9700285e-16
//------------------------------------------------------------------------------

/*============================== rnd100-100 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd100-100"
#define PP_KK	100		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	200		// Number of equations (number of rows in *.mtx)
#define PP_N	300		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 751657.607965730829164386 
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	1E-10	// Accuracy of belonging to hyperplane
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 90
// Computed objective value: 751657.607965730829164386
// Maximal objective value:  751657.607965730829164386
// Relative error = 0
// Distance to polytope: 2.8421709e-13
//------------------------------------------------------------------------------

/*============================== rnd150-0 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd150-0"
#define PP_KK	150		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	151		// Number of equations (number of rows in *.mtx)
#define PP_N	301		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 2264900
//-----------------------------------------------------------------------------

/*============================== rnd200-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd200-0"
#define PP_KK	200		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	201		// Number of equations (number of rows in *.mtx)
#define PP_N	401		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 4019900
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------ ifdef _DEBUG --------------------------------
#define PP_ITER_COUNT			10				// Each PP_ITER_COUNT-th iteration to be outputted inside PC_bsf_MapF(*)
#define PP_PROJECTION_COUNT		1000000			// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)
//-----------------------------------------------------------------------------

/*============================== rnd250-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd250-0"
#define PP_KK	250		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	251		// Number of equations (number of rows in *.mtx)
#define PP_N	501		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 6274900
//-----------------------------------------------------------------------------

/*============================== rnd400-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd400-0"
#define PP_KK	400		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	401		// Number of equations (number of rows in *.mtx)
#define PP_N	801		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 16039900
//-----------------------------------------------------------------------------

/*============================== rnd600-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd600-0"
#define PP_KK	600		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	601		// Number of equations (number of rows in *.mtx)
#define PP_N	1201	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 36059900
//-----------------------------------------------------------------------------

/*============================== rnd800-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd800-0"
#define PP_KK	800		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	801		// Number of equations (number of rows in *.mtx)
#define PP_N	1601	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 64079900
//-----------------------------------------------------------------------------

/*============================== tcube1K LP problem =========================*
#define PP_PROBLEM_NAME	"tcube1K" // Truncated hypercube
#define PP_KK	1000		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	1001		// Number of equations (number of rows in *.mtx)
#define PP_N	2001	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 100099900
//-----------------------------------------------------------------------------
// Elapsed time: 21
// Number of iterations: 6
// Computed objective value: 100099899.999809265136719
// Maximal objective value:  100099900
// Relative error = 1.91e-12
// Distance to polytope: 0
//-----------------------------------------------------------------------------

/*============================== tcube1K5 LP problem =========================*
#define PP_PROBLEM_NAME	"tcube1K5" // Truncated hypercube
#define PP_KK	1500	// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	1501	// Number of equations (number of rows in *.mtx)
#define PP_N	3001	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 225149900
//-----------------------------------------------------------------------------

/*============================== tcube2K LP problem =========================*
#define PP_PROBLEM_NAME	"tcube2K" // Truncated hypercube
#define PP_KK	2000	// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	2001	// Number of equations (number of rows in *.mtx)
#define PP_N	4001	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 400199900
//-----------------------------------------------------------------------------

/*=============================================================================*/