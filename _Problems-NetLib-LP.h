/*==========================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Simplex (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/NetLib-LP
============================================================================*/
#pragma once

#define PP_MPS_FORMAT

//=========================== Problem Parameters ===============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
#define PP_MAX_PSEUDOPROJECTING_ITER	100000000		// Maximum acceptable number of iterations in Pseudoprojection on flat

/*============================== adlittle LP problem ===========================*/
// Number of equations: 15
// Subspace dimension: 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		-225494.96316238038228101176621492
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8		// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
/// Elapsed time: 0
// Number of iterations: 4
// Computed objective value: -225494.963164675747975707
// Maximal objective value:  -225494.963162380387075245
// Relative error = 1.02e-11
// Distance to polytope: 5.8608178e-09
//------------------------------------------------------------------------------

/*============================== afiro LP problem ==============================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27			// Number of constraints in mps-file
#define PP_N 32			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11		// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE			1E-10		// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 464.753142857050022485055
// Maximal objective value:  464.753142857142847788054
// Relative error = 2e-13
// Distance to polytope: 1.2893561e-11
//------------------------------------------------------------------------------

/*============================== beaconfd LP problem ===========================*
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE			1E-7	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 16
// Computed objective value: -33592.4858072347997222096
// Maximal objective value:  -33592.4858071999988169409
// Relative error = 1.04e-12
// Distance to polytope: 2.7633211e-09
//------------------------------------------------------------------------------

/*============================== blend LP problem ==============================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74			// Number of constraints in mps-file
#define PP_N 83			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9		// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 321
// Computed objective value: 30.8121498458278431087365
// Maximal objective value:  30.8121498458282196963864
// Relative error = 1.22e-14
// Distance to polytope: 5.0086275e-14
//------------------------------------------------------------------------------

/*============================== fit1d LP problem ==============================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24		// Number of constraints 
#define PP_N 1026	// Number of variables
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
//#define PP_EPS_ZERO				1E-4	// Accuracy for comparison with zero
//#define PP_EPS_ON_HYPERPLANE		1E-6	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 273
// Number of iterations: 305
// Computed objective value: 9146.37811542919916973915
// Maximal objective value:  9146.37809242092771455646
// Relative error = 2.52e-09
// Distance to polytope: 3.8550375e-08
//------------------------------------------------------------------------------

/*============================== grow7 LP problem ============================*
// Number of equations: 140
// Subspace dimension: 161
#define PP_PROBLEM_NAME		"grow7"
#define PP_M 140	// Number of equations (after conversion to standard form)
#define PP_N 301	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 47787811.814711502616766956242865	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE			1E-5	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH		1E+6	// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 1
// Number of iterations: 11
// Computed objective value: 47787811.8147159442305565
// Maximal objective value:  47787811.8147115036845207
// Relative error = 9.29e-14
// Distance to polytope: 7.2316328e-06
//------------------------------------------------------------------------------

/*============================== israel LP problem =========================*
// Number of equations: 0
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 896644.821862765122205019
// Maximal objective value:  896644.821863045683130622
// Relative error = 3.13e-13
// Distance to polytope: 2.1396484e-10
//------------------------------------------------------------------------------

/*============================== kb2 LP problem ================================*
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-6	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 16
// Computed objective value: 1749.90011625635452219285
// Maximal objective value:  1749.90012990620562050026
// Relative error = 7.8e-09
// Distance to polytope: 3.789716e-10
//------------------------------------------------------------------------------

/*============================== recipe LP problem =============================*
// Number of equations: 67
// Subspace dimension: 92 
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91	// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_BASIC_VECTORS_ONLY
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 14
// Computed objective value: 266.616000000034091499401
// Maximal objective value:  266.615999999999985448085
// Relative error = 1.28e-13
// Distance to polytope: 1.0569e-12
//------------------------------------------------------------------------------

/*============================== sc105 LP problem ==============================*
// Number of equations: 45
// Subspace dimension: 58
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in mps-file
#define PP_N 103	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-10	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 32
// Computed objective value: 52.2020612117204976243556
// Maximal objective value:  52.2020612117072460023337
// Relative error = 2.54e-13
// Distance to polytope: 2.4098883e-11
//------------------------------------------------------------------------------

/*============================== sc50a LP problem ==============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//-------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-10	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//----------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 64.5750770586596019029457
// Maximal objective value:  64.5750770585645028631916
// Relative error = 1.47e-12
// Distance to polytope: 1.9414848e-11
//----------------------------------------------------------------------------

/*============================== sc50b LP problem ============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-11	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//--------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 70
// Maximal objective value:  70
// Relative error = 0
// Distance to polytope: 5.4835683e-14
//--------------------------------------------------------------------------

/*============================== share2b LP problem ============================*
// Number of equations: 13
// Subspace dimension: 66
#define PP_PROBLEM_NAME		"share2b"
#define PP_M 96	// Number of constraints in *.mps
#define PP_N 79	// Number of variables in *.mps
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 8
// Computed objective value: 415.732256237170588519803
// Maximal objective value:  415.732240741419502683129
// Relative error = 3.73e-08
// Distance to polytope: 1.3048846e-07
//------------------------------------------------------------------------------

/*============================== scagr7 LP problem =============================*
// Number of equations : 84
// Subspace dimension : 56
#define PP_PROBLEM_NAME	"scagr7"
#define PP_M 129		// Number of constraints in mps-file
#define PP_N 140		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 2331389.824330984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to hyperplane
//-------------------------- Compilation Modes ---------------------------------
#define PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 30
// Computed objective value: 2331389.82436896860599518
// Maximal objective value:  2331389.8243309841491282
// Relative error = 1.63e-11
// Distance to polytope: 3.8679237e-08
//------------------------------------------------------------------------------

/*============================== stocfor1 LP problem ============================*
// Number of equations: 63
// Subspace dimension: 48
#define PP_PROBLEM_NAME		"stocfor1"	
#define PP_M 117	// Number of constraints in mps-file
#define PP_N 111	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 41131.976219436406065682760731514 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		42	// Precision for point to be in halfspace
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 57
// Computed objective value: 41131.9762194380891742185
// Maximal objective value:  41131.9762194364084280096
// Relative error = 4.09e-14
// Distance to polytope: 2.4025777e-12
//------------------------------------------------------------------------------

//==============================================================================*/