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
#define PP_GRADIENT

//------------------------------------------------------------------------------
// If the jump length is very tiny then you should to decrease the PP_EPS_ZERO!
//------------------------------------------------------------------------------

/*============================== adlittle LP problem ===========================*
// Number of equations: 15
// Subspace dimension: 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		-225494.96316238038228101176621492
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------- lp_adlittle_v - VeRSAl.mtx -------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 59
// Computed objective value: -225494.963162381172878668
// Maximal objective value:  -225494.963162380387075245
// Relative error = 3.48e-15
// Distance to polytope: 4.0489652e-13
//------------------------------------------------------------------------------

/*============================== afiro LP problem ==============================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27			// Number of constraints in mps-file
#define PP_N 32			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 4
// Computed objective value: 464.753142857142904631473
// Maximal objective value:  464.753142857142847788054
// Relative error = 1.22e-16
// Distance to polytope: 7.9945703e-14
//------------------------------------------------------------------------------

/*============================== agg LP problem ================================*
// Number of equations : 36
// Subspace dimension : 127
#define PP_PROBLEM_NAME		"agg"
#define PP_M 488	// Number of constraints in mps-file
#define PP_N 163	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 35991767.286576506712640824319636
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-8	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+1	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 713
// Computed objective value: 35991767.286576509475708
// Maximal objective value:  35991767.286576509475708
// Relative error = 0
// Distance to polytope: 2.565356e-09
//------------------------------------------------------------------------------

/*============================== agg2 LP problem ===============================*
// Number of equations : 60
// Subspace dimension : 242
#define PP_PROBLEM_NAME		"agg2"
#define PP_M 516	// Number of constraints in mps-file
#define PP_N 302	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		20239252.355977109024317661926133
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-8	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+1	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 3      Simplex noMPI
// Number of iterations: 133
// Computed objective value: 20239252.3559771180152893
// Maximal objective value:  20239252.3559771105647087
// Relative error = 3.68e-16
// Distance to polytope: 2.06631e-11
//------------------------------------------------------------------------------

/*============================== beaconfd LP problem ===========================*
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+4	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 23
// Computed objective value: -33592.4858072000060928985
// Maximal objective value:  -33592.4858071999988169409
// Relative error = 2.17e-16
// Distance to polytope: 2.246e-13
//------------------------------------------------------------------------------

/*============================== blend LP problem ==============================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74			// Number of constraints in mps-file
#define PP_N 83			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-8	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+7	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 70
// Computed objective value: 30.8121498458278288978818
// Maximal objective value:  30.8121498458282196963864
// Relative error = 1.27e-14
// Distance to polytope: 2.1884247e-13
//------------------------------------------------------------------------------

/*============================== e226 LP problem ===============================*
// Number of equations: 33
// Subspace dimension : 249
#define PP_PROBLEM_NAME		"e226"
#define PP_M 223	// Number of constraints in mps-file
#define PP_N 282	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 18.751929066370549102605687681285
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-7	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+7	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
#define PP_RND_SEED 0
//-------------------------- lp_e226_v - VeRSAl.mtx ----------------------------
// Elapsed time: 5
// Number of iterations: 246
// Computed objective value: 18.7519290663703515065208
// Maximal objective value:  18.7519290663705504584868
// Relative error = 1.06e-14
// Distance to polytope: 3.5186531e-13
//------------------------------------------------------------------------------

/*============================== fit1d LP problem ==============================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24		// Number of constraints 
#define PP_N 1026	// Number of variables
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//------------------------------------------------------------------------------
// If the jump length is very tiny then you should to decrease the PP_EPS_ZERO!
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-7	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+5	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
#define PP_LOAD_BASIS
//-------------------------- lp_fit1d_v - zero.mtx -----------------------------
// Elapsed time: 529    Simplex noMPI
// Number of iterations: 478
// Computed objective value: 9146.3780924202819733182
// Maximal objective value:  9146.37809242092771455646
// Relative error = 7.06e-14
// Distance to polytope: 2.3614169e-12
//------------------------------------------------------------------------------

/*============================== grow7 LP problem ==============================*
// Number of equations: 140
// Subspace dimension: 161
#define PP_PROBLEM_NAME		"grow7"
#define PP_M 140	// Number of equations (after conversion to standard form)
#define PP_N 301	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 47787811.814711502616766956242865	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-6	// Accuracy for comparison with zero
#define PP_EPS_INVERSE			1E-10	// Accuracy for comparison with zero when calculating inverse matrix	
#define PP_EPS_ON_HYPERPLANE	1E-7	// Accuracy of belonging to hyperplane
#define PP_EPS_RELATIVE_ERROR	1E-8	// Acceptable error for optimum of objective function
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+1	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_RND_SEED 6
//#define PP_LOAD_BASIS
//------------------------------------------------------------------------------
// ?
//------------------------------------------------------------------------------

/*============================== grow15 LP problem =============================*
// Number of equations: 300
// Subspace dimension: 345
#define PP_PROBLEM_NAME		"grow15"
#define PP_M 300	// Number of equations (after conversion to standard form)
#define PP_N 645	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 106870941.29357533671604040930313	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	1E-2	// Accuracy of belonging to hyperplane
#define PP_EPS_RELATIVE_ERROR	1E-8	// Acceptable error for optimum of objective function
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1		// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- lp_grow7_v - zero.mtx -----------------------------
// ?
//------------------------------------------------------------------------------

/*============================== israel LP problem =============================*
// Number of equations: 0
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
#define PP_EPS_RELATIVE_ERROR 1E-8 // Acceptable error for optimum of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10	// Accuracy for comparison with zero
//-------------------------- lp_israel_v - VeRSAl.mtx --------------------------
// Elapsed time: 0
// Number of iterations: 160
// Computed objective value: 896644.82186304300557822
// Maximal objective value:  896644.821863045683130622
// Relative error = 2.99e-15
// Distance to polytope: 1.525458e-12
//------------------------------------------------------------------------------

/*============================== kb2 LP problem ================================*
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
#define PP_EPS_RELATIVE_ERROR 1E-8 // Acceptable error for optimum of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+5	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
#undef PP_NORMALIZATION
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 78
// Computed objective value: 1749.90012990620766686334
// Maximal objective value:  1749.90012990620562050026
// Relative error = 1.17e-15
// Distance to polytope: 9.5655205e-14
//------------------------------------------------------------------------------

/*============================== lotfi LP problem ==============================*
// Number of equations: 95
// Subspace dimension: 213
#define PP_PROBLEM_NAME		"lotfi"
#define PP_M 153	// Number of equations (after conversion to standard form)
#define PP_N 308	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 25.26470606188
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-6	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+7	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- lp_lotfi_v - VeRSAl.mtx ---------------------------
// Elapsed time: 7.9372972      Simplex MPI
// Number of iterations: 131
// Computed objective value: 25.2647060618771028828178
// Maximal objective value:  25.264706061879998344466
// Relative error = 1.15e-13
// Distance to polytope: 6.4115615e-11
//------------------------------------------------------------------------------

/*============================== recipe LP problem =============================*
// Number of equations: 67
// Subspace dimension: 92 
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91	// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//-------------------------- p_recipe_v - VeRSAl.mtx ---------------------------
// Elapsed time: 0
// Number of iterations: 17
// Computed objective value: 266.616000000000269665179
// Maximal objective value:  266.615999999999985448085
// Relative error = 1.07e-15
// Distance to polytope: 1.0255801e-15
//------------------------------------------------------------------------------

/*============================== sc105 LP problem ==============================*
// Number of equations: 45
// Subspace dimension: 58
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in mps-file
#define PP_N 103	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+7	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 219
// Computed objective value: 52.2020612117072460023337
// Maximal objective value:  52.2020612117072460023337
// Relative error = 0
// Distance to polytope: 3.2818563e-14
//------------------------------------------------------------------------------

/*============================== sc50a LP problem ==============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO			1E-10	// Accuracy for comparison with zero
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 7
// Computed objective value: 64.5750770585645028631916
// Maximal objective value:  64.5750770585645028631916
// Relative error = 0
// Distance to polytope: 1.6409282e-14
//------------------------------------------------------------------------------

/*============================== sc50b LP problem ==============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+7	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 5
// Computed objective value: 70
// Maximal objective value:  70
// Relative error = 0
// Distance to polytope: 3.9583342e-14
//------------------------------------------------------------------------------

/*============================== scagr7 LP problem =============================*
// Number of equations: 84
// Subspace dimension : 56
#define PP_PROBLEM_NAME	"scagr7"
#define PP_M 129		// Number of constraints in mps-file
#define PP_N 140		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 2331389.824330984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+2	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 30
// Computed objective value: 2331389.82433097949251533
// Maximal objective value:  2331389.8243309841491282
// Relative error = 2e-15
// Distance to polytope: 1.0989728e-12
//------------------------------------------------------------------------------

/*============================== share2b LP problem ============================*
// Number of equations: 13
// Subspace dimension: 66
#define PP_PROBLEM_NAME		"share2b"
#define PP_M 96	// Number of constraints in *.mps
#define PP_N 79	// Number of variables in *.mps
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+4	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 24
// Computed objective value: 415.732240741420582708088
// Maximal objective value:  415.732240741419502683129
// Relative error = 2.6e-15
// Distance to polytope: 1.156e-13
//------------------------------------------------------------------------------

/*============================== stocfor1 LP problem ============================*/
// Number of equations: 63
// Subspace dimension: 48
#define PP_PROBLEM_NAME		"stocfor1"	
#define PP_M 117	// Number of constraints in mps-file
#define PP_N 111	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 41131.976219436406065682760731514 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-10	// Accuracy for comparison with zero
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR			1E+4	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- lp_stocfor1_v - VeRSAl.mtx ------------------------
// Elapsed time: 0      Simplex noMPI
// Number of iterations: 22
// Computed objective value: 41131.9762194368740892969
// Maximal objective value:  41131.9762194364084280096
// Relative error = 1.13e-14
// Distance to polytope: 1.843625e-12
//------------------------------------------------------------------------------

//==============================================================================*/