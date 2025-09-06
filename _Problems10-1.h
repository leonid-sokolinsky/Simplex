/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Simplex (No MPI)
Module: _Problems10-1.h (LP problems of dimension 10 with 1 random inequality: LP-Rnd-Problems Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems were obtained using LPP-Generator https://github.com/leonid-sokolinsky/LPP-Generator
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Rnd-LP
================================================================================*/
#pragma once

//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
#define PP_MAXPROJECTION
//------------------------------------------------------------------------------

//============================== Problem Parameters =============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000				// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)

//-----------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11			// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-3			// Used if defined PP_CHECK_MAX_OBJ_VALUE 
//-----------------------------------------------------------------------------
#define PP_M	11		// Number of equations (number of rows in *.mtx)
#define PP_N	21		// Number of variables (number of cols in *.mtx)
//-----------------------------------------------------------------------------

/*============================== rnd10-0 LP problem ==============================*/
#define PP_PROBLEM_NAME	"rnd10-0"
#define PP_MAX_OBJ_VALUE 10900			
//-----------------------------------------------------------------------------
// Elapsed time: 0.0010709
// Number of iterations: 0
// Computed objective value: 10900
// Maximal objective value:  10900
// Relative error = 0
//-----------------------------------------------------------------------------

/*============================== rnd10-1-1 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-1"
#define PP_MAX_OBJ_VALUE 9551.382889057777
//-----------------------------------------------------------------------------

/*============================== rnd10-1-2 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-2"
#define PP_MAX_OBJ_VALUE 9123.80496737513
//-----------------------------------------------------------------------------

/*============================== rnd10-1-3 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-3"
#define PP_MAX_OBJ_VALUE 9960.789316531023
//-----------------------------------------------------------------------------

/*============================== rnd10-1-4 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-4"
#define PP_MAX_OBJ_VALUE 9440.134567875428
//-----------------------------------------------------------------------------

/*============================== rnd10-1-5 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-5"
#define PP_MAX_OBJ_VALUE 10248.35536348364
//-----------------------------------------------------------------------------

/*============================== rnd10-1-6 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd10-1-6"
#define PP_MAX_OBJ_VALUE 10022.93578020061
//-----------------------------------------------------------------------------

/*===============================================================================*/