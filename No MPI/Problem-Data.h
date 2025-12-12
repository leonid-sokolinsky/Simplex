/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Simplex (No MPI)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
================================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Algorithm-independent data ========================
static int PD_m;					// Total number of constraints
static int PD_n;					// Space dimension
static PT_matrix_T PD_A;			// Matrix of constraint coefficients
static PT_bitscale_T PD_isEquation;	// Constraint is equation
static PT_column_T PD_b;			// Column of constant terms (right-hand parts)
static PT_vector_T PD_c;			// Gradient of Objective Function
//========================== Algorithm variables ===============================
static int PD_meq;			// Number of total constraints being equations
static int PD_meq_basis;	// Number of basic constraints being equations
static int PD_subspaceDim;	// Dimension of of support subspace (PD_n = PD_subspaceDim + PD_meq_basis)
static int PD_iterNo;		// Number of iterations
static double PD_objF_v;	// Objective function value in curerent point
static double PD_relativeError;
//========================== Algorithm structures ================================
static PT_vector_T PD_hi;			// Higher bound
static PT_vector_T PD_lo;			// Lower bound
static PT_column_T PD_norm_a;		// Column of norms of matrix rows

static PT_vector_T PD_v;			// Current vertex
static double PD_A0[PP_N][PP_N];	// A_basis_v * v = basis_b
static double PD_A0I[PP_N][PP_N];	// AI_basis_v is the inverse matrix to A_basis_v
static PT_column_T PD_u;			// Dual point

static int PD_neHyperplanes[PP_MM];	// Index of all boundary hyperplanes
static int PD_mne;					// Number of all boundary hyperplanes

static int PD_neHyperplanes_v[PP_MM];	// Total index of all boundary hyperplanes that include vertex v
static int PD_mneh_v;					// Total number of all boundary hyperplanes that include vertex v

static int PD_basis_v[PP_MM];	// Basis index of hyperplanes that include vertex v
static int PD_m_v;				// Total number of hyperplanes that include vertex v
static PT_matrix_T PD_D;		// Auxiliary matrix used in function Matrix_Rank(*)
//========================== Input/Output =====================================
static string PD_problemName;