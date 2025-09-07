/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Simplex (No MPI)
Module: Problem-bsfCode.cpp (Implementation of Problem Code)
Prefix:	PC_bsf	- BSF Predefined Problem Functions
		SF		- Shared Functionc
		PF		- Private Functions
Author: Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;
using namespace SF;
using namespace PF;

//---------------------------------- BSF Predefined Problem Functions -----------------------------

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	Vector_Copy(parameterIn.v, parameterOutP->v);
	Vector_i_Copy(parameterIn.basis_v, parameterOutP->basis_v);
}

void PC_bsf_Init(bool* success) {

	#ifdef PP_DEBUG
	if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
		cout << "PP_MM = " << PP_MM << "\tPP_N = " << PP_N << endl;
	#endif // PP_DEBUG

	PD_problemName = PP_PROBLEM_NAME;
	PD_m = 0;
	PD_n = 0;

	#ifdef PP_MPS_FORMAT
	* success = MPS___Load_Problem();
	#else
	* success = MTX__Load_Problem();
	#endif // PP_MPS_FORMAT

	if (*success == false)
		return;

	MakeColumnOfNorms(PD_A, PD_norm_a);

	for (int i = 0; i < PD_m; i++)
		if (PD_norm_a[i] < PP_EPS_ZERO) { //Degenerate equation!
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "PC_bsf_Init error: Equation/inequality " << i << " is degenerate with pecision of PP_EPS_ZERO = "
				<< PP_EPS_ZERO << "!!!" << endl;
			*success = false;
			return;
		}

	// Index of constraints being equations
	PD_meq_total = 0;
	for (int i = 0; i < PD_m; i++)
		if (PD_isEquation[i]) {
			PD_basis_v[PD_meq_total] = i;
			PD_meq_total++;
		}

	PD_meq_basic = PD_meq_total;

	if (PD_meq_total > 0) {
		int rankEq;
		Matrix__Rank(PD_basis_v, PD_meq_total, PP_EPS_ZERO, &rankEq);
		if (PD_meq_total != rankEq) {
			#ifdef PP_BASIC_VECTORS_ONLY	
			PD_meq_basic = PD_meq_total;
			List_i_Basis(PD_basis_v, &PD_meq_basic, PP_EPS_ZERO);
			#else  // !PP_BASIC_VECTORS_ONLY
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Number of equations = " << PD_meq_basic << " < rank = " << rankEq << ".\nYou should define the option PP_BASIC_VECTORS_ONLY." << endl;
			*success = false;
			return;
			#endif // PP_BASIC_VECTORS_ONLY
		}
	}

	PD_supportSubspaceDim = PD_n - PD_meq_basic;	// Dimension of the subspace of intersection of support hyperplanes

	// Index of all boundary hyperplanes
	PD_mneh = 0;
	for (int i = 0; i < PD_m; i++)
		if (!PD_isEquation[i]) {
			PD_neHyperplanes[PD_mneh] = i;
			PD_mneh++;
		}

	*success = MTX_LoadPoint(PD_v, PP_MTX_POSTFIX_V);
	if (*success == false)
		return;

	if (!PointBelongsToPolytope(PD_v, PP_EPS_ON_HYPERPLANE)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
			cout << "PC_bsf_Init error: Starting point does NOT belong to polytope with precision of PP_EPS_ON_HYPERPLANE = " <<
				PP_EPS_ON_HYPERPLANE << "!\t ObjF = " << ObjF(PD_v) << endl;
			double eps_on_polytope = PP_EPS_ON_HYPERPLANE;
			Tuning_Eps_PointBelongsToPolytope(PD_v, &eps_on_polytope);
			cout << "Starting point will belong to polytope with precision of " << eps_on_polytope << endl;
		}
		*success = false;
		return;
	}

	PD_iterNo = 0;
	PD_objF_v = ObjF(PD_v);

	MakeBasis_v(PD_v, PD_basis_v);
	if (PD_m_v != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "PC_bsf_ProcessResults error: v is NOT vertex with precision of PP_EPS_ON_HYPERPLANE = "
			<< PP_EPS_ON_HYPERPLANE << " and PP_EPS_ZERO = " << PP_EPS_ZERO << endl;
		*success = false;
		return;
	}
}

void PC_bsf_IterInit(PT_bsf_parameter_T parameter) {
	PreparationForIteration(parameter.basis_v);
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\t Basis in v: ";
	Vector_i_Print(parameter.basis_v);
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_JobDispatcher(PT_bsf_parameter_T* parameter, int* job, bool* exit, double t) {
	// Not used
}

void PC_bsf_MainArguments(int argc, char* argv[]) {
	// Not used
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* mapSuccess) {
	double lambda_min = PP_INFINITY;
	PT_vector_T lambda_DOT_y;
	static PT_vector_T y;			// Auxiliary vector
	int entryNumber = 0;

	if (OptimumIsFound(PP_EPS_ZERO)) {
		*mapSuccess = false;
		return;
	}

	reduceElem->i_star = -1;
	for (int i = 0; i < PD_m; i++) {
		if (PD_isEquation[i])
			continue;
		if (PD_u[i] < -PP_EPS_ZERO) {
			entryNumber++;
			if (entryNumber == mapElem->entryNumber) {
				reduceElem->i_star = i;
				break;
			}
		}
	}
		
	if (reduceElem->i_star == -1) {
		*mapSuccess = false;
		return;
	}

	if (mapElem->entryNumber == 1) {
		assert(reduceElem->i_star >= 0);
	}

	assert(!PD_isEquation[reduceElem->i_star]);

	bool ok = false;
	for (int j = 0; j < PD_n; j++)
		if (BSF_sv_parameter.basis_v[j] == reduceElem->i_star) {
			for (int i = 0; i < PD_n; i++)
				y[i] = -PD_AI_v[i][j];
			ok = true;
			break;
		}
	assert(ok);

	for (int j = 0; j < PD_m; j++) {
		double a_j_DOT_y = Vector_DotProduct(PD_A[j], y);
		if (a_j_DOT_y > PP_EPS_ZERO) {
			double lambda = (PD_b[j] - Vector_DotProduct(PD_A[j], BSF_sv_parameter.v)) / a_j_DOT_y;
			if (lambda < lambda_min) {
				lambda_min = lambda;
				reduceElem->j_star = j;
			}
		}
	}

	Vector_MultiplyByNumber(y, lambda_min, lambda_DOT_y);
	Vector_Addition(BSF_sv_parameter.v, lambda_DOT_y, reduceElem->v_nex);
	reduceElem->objF_nex = ObjF(reduceElem->v_nex);
} // end PC_bsf_MapF

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// Not used
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	if (BSF_sv_mpiMaster == 0)
		cout << "No MPI" << endl;
	else {
		#define CORES_PER_NODE 12
		cout << "Number of " << CORES_PER_NODE << "-core nodes: " << (BSF_sv_numOfWorkers + 1) / CORES_PER_NODE << endl;
		cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
	}

	cout << "=================================================== " << PP_METHOD_NAME << " =================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;

	#ifdef PP_MPS_FORMAT
	cout << "Input format: MPS" << endl;
	cout << "m = " << PD_m << " n = " << PD_n << endl;
	#else
	cout << "Input format: MTX (with elimination of free variables)" << endl;
	cout << "Before elimination: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After elimination:  m =\t" << PD_m << "\tn = " << PD_n << endl;
	#endif // PP_MPS_FORMAT

	#ifdef PP_BASIC_VECTORS_ONLY	
	cout << "Total number of equations: " << PD_meq_total << endl;
	cout << "Number of basic equations: " << PD_meq_basic << endl;
	#else
	cout << "Number of equations: " << PD_meq_total << endl;
	#endif // PP_BASIC_VECTORS_ONLY

	if (PD_meq_basic > 0)
		cout << "Subspace dimension: " << PD_supportSubspaceDim << endl;

	/**
	#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented" << endl;
	#else
	cout << "Map List is not Fragmented" << endl;
	#endif
	/**/

	cout << endl;

	cout << "PP_EPS_ZERO\t\t\t" << PP_EPS_ZERO << endl;
	cout << "PP_EPS_ON_HYPERPLANE\t\t" << PP_EPS_ON_HYPERPLANE << endl;

	#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	Print_Constraints();
	cout << "Obj Function:\t"; 	Vector_Print(PD_c); cout << endl;
	cout << "v =\t\t"; Vector_Print(PD_v); cout << endl;
	#endif // PP_MATRIX_OUTPUT

	cout << endl;
	cout << "ObjF = " << setprecision(PP_SETW) << ObjF(PD_v) << setprecision(PP_SETW / 2) << endl;
	cout << "Number of inequality hyperplanes including u0:\t" << Number_IncludingNeHyperplanes(PD_v, PP_EPS_ON_HYPERPLANE) << endl;

	/*DEBUG PC_bsf_ParametersOutput**
	#ifdef PP_DEBUG
	if (!PointBelongsToPolytope(PD_v, PP_EPS_ON_HYPERPLANE))
		cout << "v0 is outside feasible polytope!!!" << endl;
	else
		cout << "v0 belongs to feasible polytope." << endl;
	//cout << "Including hyperplanes:\t"; Print_HyperplanesIncludingPoint(PD_v, PP_EPS_ON_HYPERPLANE); cout << endl;
	#endif // PP_DEBUG /**/
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	cout << setprecision(PP_SETW / 2);

	cout << "=================================================" << endl;
	cout << "// Elapsed time: " << t << endl;
	cout << "// Number of iterations: " << PD_iterNo << endl;
	cout << "// Computed objective value: " << setprecision(24) << ObjF(PD_v) << endl;
	#ifdef PP_MAX_OBJ_VALUE
	cout << "// Maximal objective value:  " << PP_MAX_OBJ_VALUE << endl;
	cout << "// Relative error = " << setprecision(3) << RelativeError(PP_MAX_OBJ_VALUE, ObjF(PD_v)) << setprecision(PP_SETW / 2) << endl;
	#endif // PP_MAX_OBJ_VALUE
	cout << "// Distance to polytope: " << Distance_PointToPolytope(PD_v) << endl;
	cout << "=================================================" << endl;
	#ifdef PP_DEBUG
	//cout << "Solution point:\t"; Vector_Print(PD_v); cout << endl;
	#endif // PP_DEBUG

	#ifdef PP_SAVE_RESULT
	if (MTX_SavePoint(PD_v, PP_MTX_POSTFIX_SO))
		cout << "Calculated solution point is saved into file *.so" << endl;
	#endif // PP_SAVE_RESULT

} // end PC_bsf_ProblemOutput

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProcessResults(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	
	if (reduceCounter == 0) {
		if (PD_iterNo == 0)
			cout << "The initial vertex is optimal!" << endl;
		*exit = true;
		return;
	}

	PD_iterNo++;

	for (int j = 0; j < PD_n; j++) // Replacing basis
		if (PD_basis_v[j] == reduceResult->i_star) {
			PD_basis_v[j] = reduceResult->j_star;
			break;
		}

	Vector_Copy(reduceResult->v_nex, PD_v);
	Vector_Copy(reduceResult->v_nex, parameter->v);
	Vector_i_Copy(PD_basis_v, parameter->basis_v);

	cout << "_________________________________________________ " << PD_iterNo << " _____________________________________________________" << endl;
	cout << "ObjF = " << ObjF(PD_v) << endl;
	#ifdef PP_MATRIX_OUTPUT
	cout << "v =\t\t"; Vector_Print(PD_v); cout << endl;
	#endif // PP_MATRIX_OUTPUT
}

void PC_bsf_ProcessResults_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ProcessResults_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ProcessResults_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* exit) {
	// Not used
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x o y
	if (x->objF_nex >= y->objF_nex) {
		z->i_star = x->i_star;
		z->j_star = x->j_star;
		z->objF_nex = x->objF_nex;
		Vector_Copy((*x).v_nex, (*z).v_nex);
	}
	else {
		z->i_star = y->i_star;
		z->j_star = y->j_star;
		z->objF_nex = y->objF_nex;
		Vector_Copy((*y).v_nex, (*z).v_nex);
	}
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// Not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// Not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// Not used
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	Vector_Copy(PD_v, parameter->v);
	Vector_i_Copy(PD_basis_v, parameter->basis_v);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PD_m;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->entryNumber = i + 1;
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//---------------------------------- Shared Functions -------------------------
namespace SF {

	static inline unsigned long long BinomialCoefficient(int n, int k) { // https://habr.com/ru/articles/274689/
		unsigned long long res = 1;
		if (k > n / 2) k = n - k;
		if (k == 1)  return n;
		if (k == 0)  return 1;
		for (int i = 1; i <= k; ++i) {
			res *= n - k + i;
			res /= i;
		}
		return res;
	}

	static inline void Bitscale_Create(bool* bitscale, int m, int* hyperplanes, int mh) {
		for (int i = 0; i < m; i++)
			bitscale[i] = false;
		for (int ih = 0; ih < mh; ih++)
			bitscale[hyperplanes[ih]] = true;
	}

	static inline void Column_Zeroing(PT_vector_T u) {  // u = 0
		for (int i = 0; i < PD_m; i++) u[i] = 0;
	}

	static inline double Distance_PointToHalfspace_i(PT_vector_T x, int i) {
		double a_DoT_z_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];

		if (a_DoT_z_MinuS_b < 0) // Point belongs to halfspace
			return 0;

		return a_DoT_z_MinuS_b / PD_norm_a[i];
	}

	static inline double Distance_PointToHyperplane_i(PT_vector_T x, int i) {
		return fabs(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / PD_norm_a[i];
	}

	static inline double Distance_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_Norm(z);
	}

	static inline double Distance_PointToPolytope(PT_vector_T x) { // Measure of distance from point to polytope
		double maxDistance = 0;
		double distance;

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				distance = Distance_PointToHyperplane_i(x, i);
			else
				distance = Distance_PointToHalfspace_i(x, i);
			if (distance > 0)
				maxDistance = PF_MAX(maxDistance, distance);
		}
		return maxDistance;
	}

	static inline double DistanceSQR_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_NormSquare(z);
	}

	static inline void Flat_BipProjection(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps_projection, int maxProjectingIter, PT_vector_T w, int* exitCode) {
		PT_vector_T p;
		PT_vector_T r;
		int iterCount = 0;
		double length_r;

		Vector_Copy(v, w);
		*exitCode = 1;

		do {
			Vector_Zeroing(r);
			for (int i = 0; i < m_flat; i++) {
				OrthogonalProjectingVectorOntoHyperplane_i(w, flatHyperplanes[i], p);
				Vector_PlusEquals(r, p);
			}

			Vector_DivideEquals(r, m_flat);

			#ifdef PP_DEBUG
			PT_vector_T w_prev;
			Vector_Copy(w, w_prev);
			#endif // PP_DEBUG

			Vector_PlusEquals(w, r);

			#ifdef PP_DEBUG
			double dist = Distance_PointToPoint(w, w_prev);
			if (dist < PF_DBL_EPSILON * 10) { // Significand bit depth is exceeded!
				*exitCode = -1;
				return;
			}
			#endif // PP_DEBUG

			iterCount++;
			if (maxProjectingIter > 0) {
				if (iterCount > maxProjectingIter) {
					*exitCode = -2;
					break;
				}
			}

			length_r = Vector_Norm(r);

			/*DEBUG Flat_BipProjection**
			#ifdef PP_DEBUG
			if (iterCount % PP_PROJECTION_COUNT == 0)
				//if (BSF_sv_mpiRank == 0)
				cout << "Worker " << BSF_sv_mpiRank << ": \t Length of r = " << length_r << endl;
			#endif // PP_DEBUG /**/

		} while (length_r >= eps_projection);

		/*DEBUG Flat_BipProjection**
		#ifdef PP_DEBUG
		cout << "Flat_BipProjection: iterCount = " << iterCount << endl;
		#endif // PP_DEBUG /**/
	}

	static inline void Flat_MaxProjection(int* flatHyperplanes, int m_flat, PT_vector_T v, double eps_projection, int maxProjectingIter, PT_vector_T w, int* exitCode) {
		PT_vector_T p;
		PT_vector_T p_max;
		double max_length;
		int iterCount = 0;

		Vector_Copy(v, w);
		*exitCode = 1;

		do {
			max_length = 0;
			Vector_Zeroing(p_max);
			for (int i = 0; i < m_flat; i++) {
				OrthogonalProjectingVectorOntoHyperplane_i(w, flatHyperplanes[i], p);
				double norm_p = Vector_Norm(p);
				if (norm_p > max_length) {
					max_length = norm_p;
					Vector_Copy(p, p_max);
				}
			}

			Vector_PlusEquals(w, p_max);

			iterCount++;
			if (maxProjectingIter > 0) {
				if (iterCount > maxProjectingIter) {
					*exitCode = -2;
					break;
				}
			}

			/*DEBUG Flat_MaxProjection**
			#ifdef PP_DEBUG
			if (iterCount % PP_PROJECTION_COUNT == 0)
				//if (BSF_sv_mpiRank == 0)
				cout << "Worker " << BSF_sv_mpiRank << ": \t max_length = " << max_length << endl;
			#endif // PP_DEBUG /**/

		} while (max_length >= eps_projection);

		/*DEBUG Flat_MaxProjection**
		#ifdef PP_DEBUG
		cout << "Flat_MaxProjection: iterCount = " << iterCount << endl;
		#endif // PP_DEBUG /**/
	}

	static inline void JumpingOnPolytope(PT_vector_T startPoint, PT_vector_T jumpVector, PT_vector_T finishPoint, double eps_jump_vector_len, double eps_on_hyperplane, double eps_zero, bool* parallelHPlanes, int* exitCode) {
		PT_vector_T o;		// Oblique projection vector
		PT_vector_T o_min;	// Oblique projection vector with minimum length
		double length_o;
		double a_DoT_d;
		double norm_d;
		double norm_a_DoT_norm_d;
		int location_z;
		double a_DoT_z_MinuS_b;
		double* d = jumpVector;		// Direction vector
		double* z = startPoint;
		double minLength_o = PP_INFINITY;

		*exitCode = 1;

		/*DEBUG JumpingOnPolytope**
		#ifdef PP_DEBUG
		cout << "d =\t"; Vector_Print(d); cout << endl;
		#endif // PP_DEBUG /**/

		norm_d = Vector_Norm(d);
		if (norm_d < eps_jump_vector_len) {
			/*DEBUG JumpingOnPolytope**
			#ifdef PP_DEBUG
			cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: norm_d < eps_jump_vector_len => return" << endl;
			#endif // PP_DEBUG /**/
			Vector_Copy(startPoint, finishPoint);
			*exitCode = -1; // Tiny length of direction vector
			return;
		}

		Vector_Zeroing(o_min);

		for (int i = 0; i < PD_m; i++) {
			if (parallelHPlanes[i])
				continue;

			a_DoT_d = Vector_DotProduct(PD_A[i], d);
			norm_a_DoT_norm_d = a_DoT_d / (PD_norm_a[i] * norm_d);

			location_z = PointLocation_i(z, i, eps_on_hyperplane, &a_DoT_z_MinuS_b);
			switch (location_z) {
			case PP_ON_HYPERPLANE:

				if (fabs(norm_a_DoT_norm_d) < eps_zero) { // Vector d is parallel to hyperplane
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i
						<< ") \tStart and finish points belong to hyperplane. => continue" << endl;
					#endif // PP_DEBUG /**/
					continue;
				}

				if (norm_a_DoT_norm_d >= eps_zero) {
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i
						<< ") \tStart point belong to hyperplane, finish point is outside half-space. => return" << endl;
					#endif // PP_DEBUG /**/
					Vector_Copy(startPoint, finishPoint);
					*exitCode = -2; // Start point belong to hyperplane, finish point is outside half-space
					return;
				}

			case PP_INSIDE_HALFSPACE:

				if (fabs(norm_a_DoT_norm_d) < eps_zero) { // Vector d is parallel to hyperplane
					/*DEBUG JumpingOnPolytope**
				#ifdef PP_DEBUG
				cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i << ") \tVector d is parallel to hyperplane. => continue" << endl;
				#endif // PP_DEBUG /**/
					continue;
				}

				if (norm_a_DoT_norm_d < 0) { // Vector looks inside of half-space
					/*DEBUG JumpingOnPolytope**
					#ifdef PP_DEBUG
					cout << "Worker " << BSF_sv_mpiRank << ": JumpingOnPolytope: " << i << ") \tVector looks inside of half-space. => continue" << endl;
					#endif // PP_DEBUG /**/
					continue;
				}

				// norm_a_DoT_norm_d > 0	// Vector looks out from half-space
				// Oblique projection vector: o = -(<a,z> - b)d/<a, d>
				Vector_MultiplyByNumber(d, -a_DoT_z_MinuS_b / a_DoT_d, o);
				length_o = Vector_Norm(o);
				if (minLength_o > length_o) {
					minLength_o = length_o;
					Vector_Copy(o, o_min);
				}
				/*DEBUG JumpingOnPolytope**
				#ifdef PP_DEBUG
				PT_vector_T u_next;
				Vector_Addition(startPoint, o, u_next);
				cout << "Worker " << BSF_sv_mpiRank << ": " << i << ") \tJump: length_o = " << length_o << ":\tu_next = "; Vector_Print(u_next);
				cout << "\tObjF(u_next) = " << ObjF(u_next) << endl;
				#endif // PP_DEBUG /**/
				break;
			case PP_OUTSIDE_HALFSPACE:
				cout << "JumpingOnPolytope error: Point is outside halfspace!" << endl;
				assert(false);
			default:
				assert(false);
			}
		}
		Vector_Addition(startPoint, o_min, finishPoint);
	}

	static void List_i_Basis(int* list_i, int* mi, double eps_zero) {
		int rank;
		int meq_total = *mi;
		int firstForProbe = 0;
		int probeRank;

		Matrix__Rank(list_i, *mi, eps_zero, &rank);

		if (*mi > rank)
			for (int check_count = 0; check_count < meq_total; check_count++) { // always check the last
				Matrix__Rank(list_i, (*mi) - 1, PP_EPS_ZERO, &probeRank);
				if (probeRank == rank)
					(*mi)--;
				else {
					int buf = list_i[firstForProbe];
					list_i[firstForProbe] = list_i[(*mi) - 1];
					list_i[(*mi) - 1] = buf;
					firstForProbe++;
				}
			}
	}

	static inline void MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a) {
		for (int i = 0; i < PD_m; i++)
			norm_a[i] = Vector_Norm(A[i]);
	}

	static inline void MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps_on_hyperplane) {
		int mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				notIncludingHalfspacesList[mo] = i;
				mo++;
			}
		if (mo < PD_m)
			notIncludingHalfspacesList[mo] = -1;
	}

	static inline void MakeNeHyperplaneList_v(PT_vector_T u, int* neHyperplanes_all, int mneh_all, int* neHyperplanes_v, int* mneh_v, double eps_on_hyperplane) {
		// List of hyperplanes that are not equations and include point u.
		*mneh_v = 0;
		for (int i = 0; i < mneh_all; i++) {
			if (PointBelongsToHyperplane_i(u, neHyperplanes_all[i], eps_on_hyperplane)) {
				neHyperplanes_v[*mneh_v] = neHyperplanes_all[i];
				(*mneh_v)++;
			}
		}
	}

	static void Matrix__Rank(int* list_i, int mi, double eps_zero, int* rank) {
		double factor;
		int maxRank = PF_MIN(PD_n, mi);

		if (mi == 1) {
			*rank = 1;
			return;
		}

		Matrix_CopyConstraintsToD(list_i, mi);

		/* Debug Matrix__Rank **
		#define WIDTH 9
		cout << setprecision(2);
		cout << "PD_n = " << PD_n << "\t mi = " << mi << endl;
		cout << "------------------" << endl;
		for (int i = 0; i < mi; i++) {
			cout << i << ")\t";
			for (int j = 0; j < PD_n; j++)
				cout << setw(WIDTH) << PD_D[i][j];
			cout << endl;
		}
		cout << "------------------" << endl;
		/* Ebd Debug */

		*rank = 0;
		for (int j = 0; j < maxRank; j++) {

			if (fabs(PD_D[j][j]) <= eps_zero) {

				int i_ne_0 = -1; // Search for a non-zero below
				for (int i = j + 1; i < mi; i++) {

					if (i_ne_0 == -1) {
						if (fabs(PD_D[i][j]) > eps_zero)
							i_ne_0 = i;
					}
					else {
						if (fabs(PD_D[i][j]) > fabs(PD_D[i_ne_0][j]))
							i_ne_0 = i;
					}
				}

				if (i_ne_0 > -1) {
					for (int lj = 0; lj < PD_n; lj++) { // Permute the rows
						double buf = PD_D[j][lj];
						PD_D[j][lj] = PD_D[i_ne_0][lj];
						PD_D[i_ne_0][lj] = buf;
					}

					/* Debug Matrix__Rank **
					cout << "------------------" << endl;
					for (int i = 0; i < mi; i++) {
						cout << i << ")\t";
						for (int j = 0; j < PD_n; j++)
							cout << setw(WIDTH) << PD_D[i][j];
						cout << endl;
					}
					cout << "------------------" << endl;
					/* End Debug */
				}
				else {

					int j_ne_0 = -1; // Search for a non-zero on the right
					for (int lj = j + 1; lj < PD_n; lj++) {
						if (j_ne_0 == -1) {
							if (fabs(PD_D[j][lj]) > eps_zero)
								j_ne_0 = lj;
						}
						else {
							if (fabs(PD_D[j][lj]) > fabs(PD_D[j][j_ne_0]))
								j_ne_0 = lj;
						}
					}

					if (j_ne_0 == -1)
						continue;

					for (int i = 0; i < mi; i++) { // Permute the columns
						double buf = PD_D[i][j];
						PD_D[i][j] = PD_D[i][j_ne_0];
						PD_D[i][j_ne_0] = buf;
					}

					/* Debug Matrix__Rank **
					cout << "------------------" << endl;
					for (int i = 0; i < mi; i++) {
						cout << i << ")\t";
						for (int j = 0; j < PD_n; j++)
							cout << setw(WIDTH) << PD_D[i][j];
						cout << endl;
					}
					cout << "------------------" << endl;
					/* End Debug */
				}
			}


			if (PD_D[j][j] != 1) {
				factor = PD_D[j][j];
				for (int lj = j; lj < PD_n; lj++)
					if (PD_D[j][lj] != 0)
						PD_D[j][lj] = PD_D[j][lj] / factor;

				/* Debug Matrix__Rank **
				cout << "------------------" << endl;
				for (int i = 0; i < mi; i++) {
					cout << i << ")\t";
					for (int j = 0; j < PD_n; j++)
						cout << setw(WIDTH) << PD_D[i][j];
					cout << endl;
				}
				cout << "------------------" << endl;
				/* End Debug */
			}

			for (int i = j + 1; i < mi; i++) {
				factor = PD_D[i][j];
				for (int lj = j; lj < PD_n; lj++)
					if (PD_D[j][lj] != 0)
						PD_D[i][lj] = PD_D[i][lj] - PD_D[j][lj] * factor;
			}
			*rank = *rank + 1;

			/* Debug Matrix__Rank **
			cout << "Rank = " << *rank << endl;
			cout << "------------------" << endl;
			for (int i = 0; i < mi; i++) {
				cout << i << ")\t";
				for (int j = 0; j < PD_n; j++)
					cout << setw(WIDTH) << PD_D[i][j];
				cout << endl;
			}
			cout << "------------------" << endl;
			/* End Debug */
		}

		/* Debug Matrix__Rank **
		cout << "Rank = " << *rank << endl;
		cout << "------------------" << endl;
		for (int i = 0; i < mi; i++) {
			cout << i << ")\t";
			for (int j = 0; j < PD_n; j++)
				cout << setw(WIDTH) << PD_D[i][j];
			cout << endl;
		}
		cout << "------------------" << endl;
		/* End Debug */
	} // End Matrix__Rank

	static void Matrix_CopyConstraintsToD(int* list_i, int count) {
		assert(count <= PD_m);
		for (int i = 0; i < count; i++)
			for (int j = 0; j < PD_n; j++)
				PD_D[i][j] = PD_A[list_i[i]][j];
	}

	static bool MPS___Load_Problem() {
		const char* mpsFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MPS_file;
		string word;
		int n_row = 0;
		int n_col = 0;
		int n_vp = 0;
		int n_fx = 0;

		PT_MPS_row_T* rows = (PT_MPS_row_T*)calloc(PP_MAX_NUMBER_OF_ROWS, sizeof(PT_MPS_row_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'rows'." << endl;
			return false;
		}

		PT_MPS_column_T* columns = (PT_MPS_column_T*)calloc(PP_MAX_NUMBER_OF_COLS, sizeof(PT_MPS_column_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'columns'." << endl;
			return false;
		}

		double* loBounds = (double*)calloc(PP_N, sizeof(double));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'loBounds'." << endl;
			return false;
		}

		PT_MPS_vpBound_T* upBounds = (PT_MPS_vpBound_T*)calloc(PP_N, sizeof(PT_MPS_vpBound_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'upBounds'." << endl;
			return false;
		}

		PT_MPS_fxVariable_T* fxVariables = (PT_MPS_fxVariable_T*)calloc(PP_N, sizeof(PT_MPS_fxVariable_T));
		if (rows == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: Can't allocate memory for array 'fxVariables'." << endl;
			return false;
		}

		MPS_file = PP_PATH;
		MPS_file += PP_MPS_PREFIX;
		MPS_file += PD_problemName;
		MPS_file += PP_MPS_EXTENSION;
		mpsFile = MPS_file.c_str();
		stream = fopen(mpsFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mpsFile << "'." << endl;
			return false;
		}

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "NAME")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'NAME'." << endl;
			return false;
		}

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "ROWS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'ROWS'." << endl;
			return false;
		}

		if (!MPS__ReadRows(stream, rows, &n_row))
			return false;

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "COLUMNS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'COLUMNS'." << endl;
			return false;
		}

		if (!MPS__ReadColumns(stream, columns, &n_col))
			return false;

		MPS__SetColumnIndexes(columns, n_col, &PD_n);

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "RHS")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem: Syntax error '" << word << "', expected 'RHS'." << endl;
			return false;
		}

		if (!MPS__ReadRHS(stream, rows, n_row))
			return false;

		MPS__SkipComment(stream);

		if (MPS__ReadKeyWord(stream, &word, "BOUNDS"))
			if (!MPS__ReadBounds(stream, columns, n_col, loBounds, upBounds, &n_vp, fxVariables, &n_fx))
				return false;

		MPS__SkipComment(stream);

		if (!MPS__ReadKeyWord(stream, &word, "ENDATA")) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS___Load_Problem error: '" << word << "', expected 'ENDATA'." << endl;
			return false;
		}

		if (!MPS__MakeProblem(rows, n_row, columns, n_col, loBounds, upBounds, n_vp, fxVariables, n_fx))
			return false;

		free(rows);
		free(columns);
		free(loBounds);
		free(upBounds);
		return true;
	}

	static inline bool MPS__MakeProblem(PT_MPS_row_T* row, int n_row, PT_MPS_column_T* column, int n_col, double* loBound,
		PT_MPS_vpBound_T* upBound, int n_vp, PT_MPS_fxVariable_T* fxVariable, int n_fx)
	{
		PD_m = 0;

		for (int i = 0; i < PP_MM; i++) {
			for (int j = 0; j < PP_N; j++)
				PD_A[i][j] = 0;
			PD_b[i] = 0;
			PD_isEquation[i] = false;
		}

		for (int j = 0; j < PP_N; j++)
			PD_c[j] = 0;

		for (int i_row = 0; i_row < n_row; i_row++) {
			switch (row[i_row].type) {
			case 'E':
				if (!MPS_AddEquation(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'G':
				if (!MPS_AddInequality_G(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'L':
				if (!MPS_AddInequality_L(row[i_row].name, row[i_row].RHS_value, column, n_col))
					return false;
				break;
			case 'N':
				if (!MPS_AddObjectiveFunction(row[i_row].name, column, n_col))
					return false;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__MakeProblem error:Unpredictable row type " << row[i_row].type << endl;
				return false;
			}
		}

		if (PD_m != PP_M) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__MakeProblem error: Number of constraints in mps-file = " << PD_m << " not equal to PP_M = " << PP_M << "." << endl;
			return false;
		}

		if (PD_n != PP_N) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__MakeProblem error: Number of variables in mps-file = " << PD_n << " not equal to PP_N = " << PP_N << "." << endl;
			return false;
		}

		for (int j = 0; j < PD_n; j++) { // Adding lower bounds
			PD_A[PD_m + j][j] = -1;
			if (loBound[j] > 0)
				PD_b[PD_m + j] = -loBound[j];
		}
		PD_m += PD_n;

		for (int j_vp = 0; j_vp < n_vp; j_vp++) { // Adding upper bounds
			PD_A[PD_m][upBound[j_vp].varIndex] = 1;
			PD_b[PD_m] = upBound[j_vp].value;
			PD_m++;
			assert(PD_m <= PP_MM);
		}

		for (int j_fx = 0; j_fx < n_fx; j_fx++) { // Adding fixed variables
			PD_A[PD_m][fxVariable[j_fx].varIndex] = 1;
			PD_b[PD_m] = fxVariable[j_fx].value;
			PD_isEquation[PD_m] = true;
			PD_m++;
			assert(PD_m <= PP_MM);
		}

		return true;
	}

	static inline bool MPS__ReadBounds(FILE* stream, PT_MPS_column_T* column, int n_col, double* loBound,
		PT_MPS_vpBound_T* upBound, int* n_vp, PT_MPS_fxVariable_T* fxVariable, int* n_fx)
	{
		char ch, typeCh1, typeCh2;
		PT_MPS_name_T baundName;
		PT_MPS_name_T varName;
		int j_var;
		double boundValue;
		fpos_t pos;	// Position in the input stream

		assert(*n_vp == 0);
		assert(*n_fx == 0);

		for (int j_lo = 0; j_lo < PD_n; j_lo++)
			loBound[j_lo] = 0;

		ch = getc(stream);
		if (ch != ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadBounds error: Expected 'space'!" << endl;
			return false;
		}

		while (ch == ' ') {
			typeCh1 = getc(stream);
			typeCh2 = getc(stream);
			ch = getc(stream);

			MPS_ReadName(stream, baundName);
			MPS_SkipSpaces(stream);
			MPS_ReadName(stream, varName);

			j_var = -1;
			for (int j_col = 0; j_col < n_col; j_col++)
				if (MPS_SameNames(column[j_col].varName, varName)) {
					j_var = column[j_col].j;
					break;
				}
			if (j_var == -1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadBounds Error: Variable " << varName << " was not found." << endl;
				return false;
			}

			MPS_SkipSpaces(stream);

			MPS_ReadValue(stream, &boundValue);

			switch (typeCh1) {
			case 'L':
				loBound[j_var] = boundValue;
				break;
			case 'U':
				upBound[*n_vp].value = boundValue;
				upBound[*n_vp].varIndex = j_var;
				(*n_vp)++;
				break;
			case 'F':
				fxVariable[*n_fx].value = boundValue;
				fxVariable[*n_fx].varIndex = j_var;
				(*n_fx)++;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadBounds Error: unpredictable bound type " << typeCh1 << endl;
				return false;
			}

			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
		}

		fsetpos(stream, &pos);
		return true;
	}

	static inline bool MPS__ReadColumns(FILE* stream, PT_MPS_column_T* column, int* n_col) {
		char ch;
		fpos_t pos;	// Position in the input stream

		*n_col = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		fsetpos(stream, &pos);
		if (ch != ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadColumns error: Expected 'space'!" << endl;
			return false;
		}

		while (ch == ' ') {
			fgetpos(stream, &pos);
			ch = getc(stream);
			fsetpos(stream, &pos);
			if (ch != ' ')
				return true;

			if (!MPS_ReadColumnLine(stream, column, n_col))
				return false;
		}

		return true;
	}

	static inline bool MPS__ReadKeyWord(FILE* stream, string* word, string pattern) {
		char str[80] = { '\0' };
		fpos_t pos;	// Position in the input stream

		fgetpos(stream, &pos);

		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "MPS__ReadKeyWord error: Failure to read string!" << endl;
			return false;
		}

		*word = str;

		if (*word != pattern) {
			fsetpos(stream, &pos);
			return false;
		}

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS__ReadRHS(FILE* stream, PT_MPS_row_T* row, int n_row) {
		PT_MPS_name_T RHS_name = { '\0' };
		fpos_t pos;	// Position in the input stream
		char ch;

		for (int i_row = 0; i_row < n_row; i_row++)
			row[i_row].RHS_value = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		if (ch == 'B') {
			fsetpos(stream, &pos);
			return true;
		}
		while (ch == ' ') {
			fsetpos(stream, &pos);
			MPS_ReadRHS_line(stream, row, n_row, RHS_name);

			fgetpos(stream, &pos);
			ch = getc(stream);
		}
		fsetpos(stream, &pos);
		return true;
	}

	static inline bool MPS__ReadRows(FILE* stream, PT_MPS_row_T* row, int* n_row) {
		char ch;
		int j;
		fpos_t pos;	// Position in the input stream

		*n_row = 0;

		fgetpos(stream, &pos);
		ch = getc(stream);
		while (ch == ' ') {
			j = 1;

			ch = getc(stream);
			j++;

			if (ch == ' ') {
				ch = getc(stream);
				j++;
			}

			if (ch <= ' ') {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: Syntax error in row " << *n_row + 1 << endl;
				return false;
			}

			switch (ch) {
			case 'E':
			case 'G':
			case 'L':
			case 'N':
				row[*n_row].type = ch;
				break;
			default:
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: Syntax error - unpredictable row type " << ch << endl;
				return false;
			}

			while (j < 4) {
				ch = getc(stream);
				if (ch < ' ') {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS__ReadRows: Syntax error - not ASCII symbol." << endl;
					return false;
				}
				j++;
			}

			MPS_SkipSpaces(stream);

			MPS_ReadName(stream, row[*n_row].name);

			(*n_row)++;
			if (*n_row > PP_MAX_NUMBER_OF_ROWS) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows: The number of row is greater than PP_MAX_NUMBER_OF_ROWS = " << PP_MAX_NUMBER_OF_ROWS << endl;
				return false;
			}

			if (!MPS_UniqueRowName(row, *n_row, row[*n_row - 1].name)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS__ReadRows:Error: The name '" << row[*n_row - 1].name << "' is not unique.";
				return false;
			}

			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
		}

		fsetpos(stream, &pos);
		return true;
	}

	static inline void MPS__SkipComment(FILE* stream) {
		fpos_t pos;	// Position in the input stream
		int res;

		fgetpos(stream, &pos);

		while (getc(stream) == '*') {
			while (getc(stream) != 10);
			res = fscanf(stream, "\n");
			fgetpos(stream, &pos);
		}

		fsetpos(stream, &pos);
	}

	static inline bool MPS_AddEquation(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddEquation error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << "." << endl;
					return false;
				}
				PD_A[PD_m][column[i_col].j] = column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddEquation error: Zero row " << rowName << " has non-zero RHS value." << endl;
				return false;
			}
			return true;
		}
		PD_b[PD_m] = RHS_value;
		PD_isEquation[PD_m] = true;
		PD_m++;
		assert(PD_m < PP_MM);
		return true;
	}

	static inline bool MPS_AddInequality_G(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddInequality_G error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << "." << endl;
					return false;
				}
				PD_A[PD_m][column[i_col].j] = -column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddInequality_G error: Zero row " << rowName << " has non-zero RHS value." << endl;
				return false;
			}
			return true;
		}
		PD_b[PD_m] = -RHS_value;
		PD_m++;
		assert(PD_m < PP_MM);
		return true;
	}

	static inline bool MPS_AddInequality_L(PT_MPS_name_T rowName, double RHS_value, PT_MPS_column_T* column, int n_col) {
		bool empty = true;

		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_A[PD_m][column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddInequality_L error: Coefficient redefinition of the variable " << column[i_col].varName
						<< " in row " << column[i_col].rowName << "." << endl;
					return false;
				}
				PD_A[PD_m][column[i_col].j] = column[i_col].value;
				empty = false;
			}
		if (empty) {
			if (RHS_value != 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_AddInequality_L error: Zero row " << rowName << " has non-zero RHS value." << endl;
				return false;
			}
			return true;
		}
		PD_b[PD_m] = RHS_value;
		PD_m++;
		assert(PD_m <= PP_MM);
		return true;
	}

	static inline bool MPS_AddObjectiveFunction(PT_MPS_name_T rowName, PT_MPS_column_T* column, int n_col) {
		for (int i_col = 0; i_col < n_col; i_col++)
			if (MPS_SameNames(column[i_col].rowName, rowName)) {
				if (!PD_c[column[i_col].j] == 0) {
					if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
						cout << "MPS_AddObjectiveFunction error: Coefficient redefinition of the objective function for " << column[i_col].varName << "." << endl;
					return false;
				}
				#ifdef PP_OPT_MIN
				PD_c[column[i_col].j] = column[i_col].value;
				#else
				PD_c[column[i_col].j] = -column[i_col].value;
				#endif
			}
		return true;
	}

	static inline void MPS_CopyName(char* name_x, char* name_y) {
		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			name_y[p] = name_x[p];
	}

	static inline void MPS_NewLine(FILE* stream) {
		char ch;
		fpos_t pos;	// Position in the input stream

		fgetpos(stream, &pos);
		if (getc(stream) == -1) // EOF
			return;
		fsetpos(stream, &pos);

		do {
			fgetpos(stream, &pos);
		} while (getc(stream) != '\n');

		fsetpos(stream, &pos);

		do {
			fgetpos(stream, &pos);
			ch = getc(stream);
		} while (ch == '\n' || ch == '\r');

		fsetpos(stream, &pos);
	}

	static inline bool MPS_ReadColumnLine(FILE* stream, PT_MPS_column_T* column, int* n_col) {
		fpos_t pos;	// Position in the input stream
		double value;

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, column[*n_col].varName))
			return false;

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, column[*n_col].rowName))
			return false;

		// Checking for duplicates
		for (int j_col = 0; j_col < *n_col; j_col++)
			if (MPS_SameNames(column[j_col].varName, column[*n_col].varName) && MPS_SameNames(column[j_col].rowName, column[*n_col].rowName)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadColumnLine error: Redefinition of item varName = " << column[j_col].varName
					<< ", rowName = " << column[j_col].rowName << "." << endl;
				return false;
			}

		MPS_SkipSpaces(stream);

		if (!MPS_ReadValue(stream, &value))
			return false;

		column[*n_col].value = value;

		(*n_col)++;
		if (*n_col > PP_MAX_NUMBER_OF_COLS) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadColumnLine error: The number of column is greater than PP_MAX_NUMBER_OF_COLS = " << PP_MAX_NUMBER_OF_COLS << endl;
			return false;
		}

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			fsetpos(stream, &pos);
			MPS_NewLine(stream);
			return true;
		};
		fsetpos(stream, &pos);

		MPS_CopyName(column[*n_col - 1].varName, column[*n_col].varName);

		if (!MPS_ReadName(stream, column[*n_col].rowName))
			return false;

		// Checking for duplicates
		for (int j_col = 0; j_col < *n_col; j_col++)
			if (column[j_col].varName == column[*n_col].varName && column[j_col].rowName == column[*n_col].rowName) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadColumnLine error: Redefinition of item varName = " << column[j_col].varName
					<< ", rowName = " << column[j_col].rowName << "." << endl;
				return false;
			}

		MPS_SkipSpaces(stream);

		if (!MPS_ReadValue(stream, &value))
			return false;

		column[*n_col].value = value;

		(*n_col)++;
		if (*n_col > PP_MAX_NUMBER_OF_COLS) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS__ReadColumns: The number of column is greater than PP_MAX_NUMBER_OF_COLS = " << PP_MAX_NUMBER_OF_COLS << endl;
			return false;
		}

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS_ReadName(FILE* stream, char* name) {
		char ch;
		fpos_t pos;	// Position in the input stream

		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			name[p] = '\0';

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadName: Syntax error - not ASCII symbol." << endl;
			return false;
		}
		fsetpos(stream, &pos);

		MPS_SkipSpaces(stream);

		for (int p = 0; p < 8; p++) {
			fgetpos(stream, &pos);
			ch = getc(stream);
			if (ch <= ' ') {
				fsetpos(stream, &pos);
				break;
			}
			name[p] = ch;
		}

		return true;
	}

	static inline bool MPS_ReadRHS_line(FILE* stream, PT_MPS_row_T* row, int n_row, PT_MPS_name_T RHS_name)
	{ // RHS_name must be equal to {'\0'} when calling MPS_ReadRHS_line at the first time
		fpos_t pos;	// Position in the input stream
		char ch;
		PT_MPS_name_T next_RHS_name;
		PT_MPS_name_T rowName;
		double RHS_value;
		int rowIndex;

		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++)
			next_RHS_name[p] = '\0';

		for (int p = 0; p < 4; p++) {
			ch = getc(stream);
			if (ch != ' ') {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadRHS_line: Syntax error, expected ' '" << endl;
				return false;
			}
		}

		int p = 0;
		fgetpos(stream, &pos);
		while (ch == ' ') {
			fgetpos(stream, &pos);
			ch = getc(stream);
			p++;
		}
		fsetpos(stream, &pos);

		if (p > 8)
			next_RHS_name[0] = ' ';
		else
			if (!MPS_ReadName(stream, next_RHS_name))
				return false;

		if (RHS_name[0] != '\0') {
			if (!MPS_SameNames(RHS_name, next_RHS_name)) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadRHS_line: Syntax error, Multiple RHS names." << endl;
				return false;
			}
		}
		else
			MPS_CopyName(next_RHS_name, RHS_name);

		MPS_SkipSpaces(stream);

		if (!MPS_ReadName(stream, rowName))
			return false;

		rowIndex = MPS_SearchRowByName(row, n_row, rowName);

		if (rowIndex < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Syntax error, non-existent row name '" << rowName << "'." << endl;
			return false;
		}

		MPS_SkipSpaces(stream);

		if (fscanf(stream, "%lf", &RHS_value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Unexpected end of line!" << endl;
			return false;
		}
		row[rowIndex].RHS_value = RHS_value;

		if (row[rowIndex].type == 'N')
			if (RHS_value != 0)
			{
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "MPS_ReadRHS_line warning: RHS value for row of type 'N' is not equal to 0." << endl;
				return false;
			}

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (getc(stream) < ' ') {
			fsetpos(stream, &pos);
			MPS_NewLine(stream);
			fgetpos(stream, &pos);
			ch = getc(stream);
			fsetpos(stream, &pos);
			return true;
		};
		fsetpos(stream, &pos);

		if (!MPS_ReadName(stream, rowName))
			return false;

		rowIndex = MPS_SearchRowByName(row, n_row, rowName);

		if (rowIndex < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Syntax error, non-existent row name '" << rowName << "'." << endl;
			return false;
		}

		if (row[rowIndex].type == 'N') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Invalid row type 'N'." << endl;
			return false;
		}

		MPS_SkipSpaces(stream);

		if (fscanf(stream, "%lf", &RHS_value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Unexpected end of line!" << endl;
			return false;
		}
		row[rowIndex].RHS_value = RHS_value;

		MPS_SkipSpaces(stream);

		fgetpos(stream, &pos);
		if (char ch = getc(stream) > ' ') {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadRHS_line: Illegal ASCII character '" << ch << "' at end of line." << endl;
			return false;
		}
		fsetpos(stream, &pos);

		MPS_NewLine(stream);

		return true;
	}

	static inline bool MPS_ReadValue(FILE* stream, double* value) {
		if (fscanf(stream, "%lf", value) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MPS_ReadValue: Error: Non-ASCII character." << endl;
			return false;
		}
		return true;
	}

	static inline bool MPS_SameNames(PT_MPS_name_T name_x, PT_MPS_name_T name_y) {
		for (int p = 0; p < PP_MPS_NAME_LENGTH; p++) {
			if (name_x[p] == '\0' && name_y[p] == '\0')
				return true;
			if (name_x[p] != name_y[p])
				return false;
		}
		return true;
	}

	static inline void MPS__SetColumnIndexes(PT_MPS_column_T* column, int n_col, int* n) {

		*n = 0;

		for (int j_col = 0; j_col < n_col; j_col++)
			column[j_col].j = -1;

		for (int j_col = 0; j_col < n_col; j_col++) {
			if (column[j_col].j >= 0)
				continue;
			column[j_col].j = *n;
			(*n)++;
			assert(*n <= PP_N);
			for (int j = j_col + 1; j < n_col; j++) {
				if (column[j].j >= 0)
					continue;
				if (MPS_SameNames(column[j_col].varName, column[j].varName))
					column[j].j = column[j_col].j;
			}
		}
	}

	static inline int MPS_SearchColByName(PT_MPS_column_T* column, int j_col, PT_MPS_name_T name) {
		for (int j = 0; j < j_col; j++)
			if (MPS_SameNames(column[j].varName, name))
				return j;
		return -1;
	}

	static inline int MPS_SearchRowByName(PT_MPS_row_T* row, int n_row, PT_MPS_name_T name) {
		for (int i = 0; i < n_row; i++)
			if (MPS_SameNames(row[i].name, name))
				return i;
		return -1;
	}

	static inline void MPS_SkipSpaces(FILE* stream) {
		fpos_t pos;	// Position in the input stream

		do {
			fgetpos(stream, &pos);
		} while (getc(stream) == ' ');

		fsetpos(stream, &pos);
	}

	static inline bool MPS_UniqueRowName(PT_MPS_row_T* rows, int n_row, PT_MPS_name_T name) {
		int count = 0;

		for (int i = 0; i < n_row; i++)
			if (MPS_SameNames(rows[i].name, name)) {
				if (count > 0)
					return false;
				count = 1;
			}
		return true;
	}

	static void MTX_Conversion() { // Removing free variables

		for (int i = 0; i < PD_m; i++)
			PD_isEquation[i] = true;

		MTX_RemoveFreeVariables();

		for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
			for (int j = 0; j < PD_n; j++)
				PD_A[i + PD_m][j] = 0;
			PD_A[i + PD_m][i] = -1;
			if (PD_lo[i] == 0)
				PD_b[i + PD_m] = 0;
			else
				PD_b[i + PD_m] = -PD_lo[i];
		}
		PD_m += PD_n; assert(PD_m <= PP_MM);

		for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
			if (PD_hi[i] != PP_INFINITY) {
				for (int j = 0; j < PD_n; j++)
					PD_A[PD_m][j] = 0;
				PD_A[PD_m][i] = 1;
				PD_b[PD_m] = PD_hi[i];
				PD_m++; assert(PD_m <= PP_MM);
			}
		}

		/**
		cout << "-----------------------------------------------------" << endl;
		Print_Constraints();
		cout << "-----------------------------------------------------" << endl;
		cout << "PD_c: "; Vector_Print(PD_c); cout << endl;/**/
	}

	static bool MTX__Load_Problem() {

		//--------------- Reading A ------------------
		if (!MTX_Load_A())
			return false;

		//--------------- Reading b ------------------
		if (!MTX_Load_b())
			return false;

		//--------------- Reading lo ------------------
		if (!MTX_Load_lo())
			return false;

		//--------------- Reading c ------------------
		if (!MTX_Load_c())
			return false;

		//--------------- Reading hi ------------------
		if (!MTX_Load_hi())
			return false;

		//---------- Conversion to inequalities -----------
		MTX_Conversion();

		/**
		cout << "-----------------------------------------------------" << endl;
		Print_Constraints();
		cout << "-----------------------------------------------------" << endl;
		cout << "PD_c: "; Vector_Print(PD_c); cout << endl;/**/

		return true;
	}

	static inline bool MTX_Load_A() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		int non;	// Number of non-zero elements
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_A;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file " << mtxFile << endl;
			return false;
		}

		if (nor >= noc) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "" << endl;
			return false;
		}

		if (noc != PP_N) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MTX_Load_A error: PP_N must be = " << noc << "" << endl;
			return false;
		}

		if (nor != PP_M) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "MTX_Load_A error:  PP_M must be = " << nor << "" << endl;
			return false;
		}

		PD_m = nor;
		PD_n = noc;

		if (nor + noc > PP_MM) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Invalid input data: number of constraints m = " << nor + noc
				<< " must be < PP_MM + 1 =" << PP_MM + 1 << "" << endl;
			return false;
		}

		for (int k = 0; k < non; k++) {
			int i, j;

			if (fscanf(stream, "%d%d%s", &i, &j, str) < 3) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file'" << mtxFile << "'." << endl;
				return false;
			}

			i -= 1;
			j -= 1;
			if (i < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Negative row index in'" << mtxFile << "'." << endl;
				return false;
			}
			if (j < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Negative column index in'" << mtxFile << "'." << endl;
				return false;
			}
			PD_A[i][j] = strtod(str, &chr);
		}

		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_b() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_B;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (PD_m != nor) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'." << endl;
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'." << endl;
			return false;
		}

		for (int i = 0; i < PD_m; i++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_b[i] = strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_c() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_C;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'." << endl;
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'." << endl;
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file" << endl;
				return false;
			}
			PD_c[j] = -strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline bool MTX_Load_hi() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_HI;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'." << endl;
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'." << endl;
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout << "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_hi[j] = strtod(str, &chr);
		}
		fclose(stream);
		return true;
	}

	static inline bool MTX_Load_lo() {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += PP_MTX_POSTFIX_LO;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'." << endl;
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'." << endl;
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 1) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file '" << mtxFile << "'." << endl;
				return false;
			}
			PD_lo[j] = strtod(str, &chr);
		}

		fclose(stream);

		return true;
	}

	static inline bool MTX_LoadPoint(PT_vector_T x, string postfix) {
		int nor;	// Number of matrix rows
		int noc;	// Number of matrix columns
		const char* mtxFile;
		FILE* stream;// Input stream
		char str[80] = { '\0' };
		char* chr = str;
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += postfix;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "r+b");

		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		MTX_SkipComments(stream);
		if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}
		if (nor != PD_n) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of rows in'" << mtxFile << "'. Must be " << PD_n << "" << endl;
			return false;
		}
		if (noc != 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Incorrect number of columnws in'" << mtxFile << "'." << endl;
			return false;
		}

		for (int j = 0; j < PD_n; j++) {
			if (fscanf(stream, "%s", str) < 0) {
				if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
					cout
					<< "Unexpected end of file" << endl;
				return false;
			}
			x[j] = strtod(str, &chr);
		}
		fclose(stream);

		return true;
	}

	static inline void MTX_RemoveFreeVariables(void) {
		int freeVariable_j_col[PP_M];
		int m_freeVariable_j_col = 0;
		bool unique;

		// Find free variables
		for (int j_col = 0; j_col < PD_n; j_col++) { // Find zero element in PD_c
			if (PD_c[j_col] != 0)
				continue;
			if (PD_lo[j_col] != 0)
				continue;
			if (PD_hi[j_col] != PP_INFINITY)
				continue;

			for (int i_row = 0; i_row < PD_m; i_row++) { // Find PD_A i_row with nonzero element in column j_col
				if (PD_A[i_row][j_col] == 0)
					continue;

				// Check uniqueness in column j_col
				unique = true;
				for (int i = 0; i < PD_m; i++) {
					if (i == i_row)
						continue;
					if (PD_A[i][j_col] != 0) {
						unique = false;
						break;
					}
				}
				if (!unique)
					continue;

				freeVariable_j_col[m_freeVariable_j_col] = j_col;
				m_freeVariable_j_col++; assert(m_freeVariable_j_col <= PP_M);

				// Check sign of free variable
				if (PD_A[i_row][j_col] < 0) {
					// Change sign of inequality
					for (int j = 0; j < PD_n; j++)
						PD_A[i_row][j] = -PD_A[i_row][j];
					PD_b[i_row] = -PD_b[i_row];
				}

				PD_isEquation[i_row] = false;

				break;
			}
		}

		{// Eliminate columns with free variables
			static bool colToDeleteLable[PP_N];
			for (int j = 0; j < m_freeVariable_j_col; j++) {
				assert(freeVariable_j_col[j] < PP_N);
				colToDeleteLable[freeVariable_j_col[j]] = true;
			}

			for (int j = 0; j < PD_n; j++) {
				if (colToDeleteLable[j]) {
					for (int i = 0; i < PD_m; i++)
						PD_A[i][j] = PD_A[i][PD_n - 1];
					PD_c[j] = PD_c[PD_n - 1];
					colToDeleteLable[j] = colToDeleteLable[PD_n - 1];
					j--; assert(j >= 0);
					PD_n--; assert(PD_n >= 0);

					/**
					cout << "-----------------------------------------------------" << endl;
					Print_Constraints();
					cout << "-----------------------------------------------------" << endl;/**/
				}
			}
		}
	}

	static bool MTX_SavePoint(PT_vector_T x, string postfix) {
		const char* mtxFile;
		FILE* stream;// Input stream
		string MTX_file;

		MTX_file = PP_PATH;
		MTX_file += PP_MTX_PREFIX;
		MTX_file += PD_problemName;
		MTX_file += postfix;
		mtxFile = MTX_file.c_str();
		stream = fopen(mtxFile, "w");
		if (stream == NULL) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Failure of opening file '" << mtxFile << "'." << endl;
			return false;
		}

		fprintf(stream, "%d %d\n", PD_n, 1);

		for (int j = 0; j < PD_n; j++)
			fprintf(stream, "%.16f\n", x[j]);

		fclose(stream);
		return true;
	}

	static inline void MTX_SkipComments(FILE* stream) {
		fpos_t pos;	// Position in the input stream
		int res;
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
		while (getc(stream) == '%') {
			while (getc(stream) != 10);
			res = fscanf(stream, "\n");
			fgetpos(stream, &pos);
		}
		fsetpos(stream, &pos);
	}

	static inline int Number_IncludingNeHyperplanes(PT_vector_T x, double eps_on_hyperplane) {
		int number = 0;

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				number++;
		}
		return number;
	}

	static inline double ObjF(PT_vector_T x) {
		double s = 0;
		for (int j = 0; j < PD_n; j++)
			s += PD_c[j] * x[j];
		return s;
	}

	static inline void OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p) {
		double ns = Vector_NormSquare(PD_A[i]);
		Vector_MultiplyByNumber(PD_A[i], -(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / ns, p);
	}

	static inline bool PointBelongsToFlat(PT_vector_T x, int* hyperplaneList, int hyperplaneCount, double eps_on_hyperplane) { // If the point belongs to the flat
		for (int i = 0; i < hyperplaneCount; i++)
			if (!PointBelongsToHyperplane_i(x, hyperplaneList[i], eps_on_hyperplane))
				return false;
		return true;
	}

	static inline bool PointBelongsToHalfspace_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];

		/**/
		#ifdef PP_DEBUG
		if (distanceToHyperplane > eps_on_hyperplane && distanceToHyperplane < eps_on_hyperplane * 10) {
			cout << "PointBelongsToHalfspace_i: Distance from testing point is less than " << eps_on_hyperplane * 10 << ", but greater than " << eps_on_hyperplane << "!" << endl;
		}
		#endif // PP_DEBUG /**/

		if (distanceToHyperplane < eps_on_hyperplane)
			return true;
		if (PD_isEquation[i])
			return false;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline bool PointBelongsToHyperplane_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double dist = Distance_PointToHyperplane_i(x, i);

		/**/
		#ifdef PP_DEBUG
		if (dist > eps_on_hyperplane && dist < eps_on_hyperplane * 10) {
			cout << "PointBelongsToHyperplane_i: Distance from testing point is less than " << eps_on_hyperplane * 10 << ", but greater than " << eps_on_hyperplane << "!" << endl;
		}
		#endif // PP_DEBUG /**/

		if (dist < eps_on_hyperplane)
			return true;
		else
			return false;
	}

	static inline bool PointBelongsToPolytope(PT_vector_T x, double eps_on_hyperplane) { // If the point belongs to the polytope with prescigion of eps_on_hyperplane
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane))
				return false;
		return true;
	}

	static inline bool PointInsideHalfspace_i(PT_vector_T x, int i, double eps_on_hyperplane) {
		double a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_DoT_x_MinuS_b) / PD_norm_a[i];

		#ifdef PP_DEBUG
		if (distanceToHyperplane > eps_on_hyperplane && distanceToHyperplane < eps_on_hyperplane * 10) {
			cout << "PointInsideHalfspace_i: Distance from testing point is less than " << eps_on_hyperplane * 10 << ", but greater than " << eps_on_hyperplane << "!" << endl;
		}
		#endif // PP_DEBUG /**/

		if (distanceToHyperplane < eps_on_hyperplane)
			return false;
		if (a_DoT_x_MinuS_b < 0)
			return true;
		return false;
	}

	static inline bool PointIsBoundary(PT_vector_T x, double eps_on_hyperplane) {
		assert(PointBelongsToPolytope(x, eps_on_hyperplane));

		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane))
				return true;
		}
		return false;
	}

	static inline int PointLocation_i(PT_vector_T x, int i, double eps_on_hyperplane, double* a_DoT_x_MinuS_b) {
		*a_DoT_x_MinuS_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double dist = fabs(*a_DoT_x_MinuS_b) / PD_norm_a[i];

		/**/
		#ifdef PP_DEBUG
		if (dist > eps_on_hyperplane && dist < eps_on_hyperplane * 10) {
			cout << "PointLocation_i: Distance from testing point is less than " << eps_on_hyperplane * 10 << ", but greater than " << eps_on_hyperplane << "!" << endl;
		}
		#endif // PP_DEBUG /**/

		if (dist < eps_on_hyperplane)// <a,x> = b
			return PP_ON_HYPERPLANE;

		if (*a_DoT_x_MinuS_b < 0)								// <a,x> < b
			return PP_INSIDE_HALFSPACE;

		return PP_OUTSIDE_HALFSPACE;							// <a,x> > b
	}

	static inline void Print_Constraints() {
		for (int i = 0; i < PD_m; i++) {
			cout << i << ")";
			for (int j = 0; j < PD_n; j++)
				cout << setw(PP_SETW) << PD_A[i][j];
			cout << "\t" << (PD_isEquation[i] ? "==" : "<=") << setw(PP_SETW) << PD_b[i] << endl;
		}
	}

	static inline void Print_HalfspacesIncludingPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HalfspacesOutOfPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (!PointBelongsToHalfspace_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline void Print_HyperplanesIncludingPoint(PT_vector_T x, double eps_on_hyperplane) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsToHyperplane_i(x, i, eps_on_hyperplane)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	static inline double RelativeError(double trueValue, double calculatedValue) {
		if (trueValue == 0)
			return fabs(calculatedValue - trueValue);
		else
			return fabs(calculatedValue - trueValue) / fabs(trueValue);
	}

	static inline void Tuning_Eps_PointBelongsToFlat(PT_vector_T x, int* hyperplaneList, int hyperplaneCount, double* eps) {
		while (!PointBelongsToFlat(x, hyperplaneList, hyperplaneCount, *eps))
			(*eps) *= 2;
	}

	static inline void Tuning_Eps_PointBelongsToPolytope(PT_vector_T x, double* eps) {
		while (!PointBelongsToPolytope(x, *eps))
			(*eps) *= 2;
	}

	static inline int TWIDDLE__BinomialCoefficient(int n, int k, int* p) { // |p|=n+2
		int x, y, z;
		bool done = false;
		int B;

		assert(n <= PP_MM);

		TWIDDLE_Make_p(p, n, k);
		B = 0;
		while (!done) {
			TWIDDLE_Run(&x, &y, &z, p, &done);
			if (B == PF_INT_MAX) cout << "TWIDDLE__BinomialCoefficient warning: value of integer variable B has exceeded PF_INT_MAX = " << PF_INT_MAX << endl;
			B++;
		}
		return B;
	}

	// Combination of m out of n
	static inline void TWIDDLE__CodeToSubset(int code, int* a, int* c, int n, int m, int* p, bool* done) {
		static int x, y, z;

		TWIDDLE_Make_p(p, n, m);
		for (int k = 0; k < m; k++)
			c[k] = a[n - m + k];

		if (code == 0) return;

		for (int i = 0; i < code; i++) {
			TWIDDLE_Run(&x, &y, &z, p, done);
			if (*done)
				return;
			c[z - 1] = a[x - 1];
		}
	}

	static inline void TWIDDLE_Make_p(int* p, int n, int m) {
		// p - auxiliary integer array for generating all combinations of m out of n objects.
		assert(n >= m && m > 0);
		p[0] = n + 1;
		p[n + 1] = -2;
		for (int j = 1; j <= n - m; j++)
			p[j] = 0;
		for (int j = n - m + 1; j <= n; j++)
			p[j] = j - n + m;
	}

	static inline void TWIDDLE_Run // https://doi.org/10.1145/362384.362502
	(int* x, int* y, int* z, int* p, bool* done) {
		int i, j, k;
		j = 0;
		*done = false;

		do {
			j++;
		} while (p[j] <= 0);

		if (p[j - 1] == 0) {
			i = j - 1;
			while (i != 1) {
				p[i] = -1;
				i -= 1;
			}
			p[j] = 0;
			p[1] = *x = *z = 1;
			*y = j;
			return;
		}

		if (j > 1)
			p[j - 1] = 0;

		do {
			j++;
		} while (p[j] > 0);

		i = k = j - 1;

		i++;
		while (p[i] == 0) {
			p[i] = -1;
			i++;
		}

		if (p[i] == -1) {
			p[i] = *z = p[k];
			*x = i;
			*y = k;
			p[k] = -1;
			return;
		}

		if (i == p[0]) {
			*done = true;
			return;
		}

		*z = p[j] = p[i];
		p[i] = 0;
		*x = j;
		*y = i;
	}

	static inline void Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint) {
		for (int j = 0; j < PD_n; j++)
			shiftedPoint[j] = point[j] + shiftVector[j] * factor;
	}

	static inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] + y[j];
	}

	static inline void Vector_Copy(PT_vector_T fromVector, PT_vector_T toVector) { // toPoint = fromPoint
		for (int j = 0; j < PD_n; j++)
			toVector[j] = fromVector[j];
	}

	static inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] / r;
	}

	static inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
		for (int j = 0; j < PD_n; j++)
			x[j] /= r;
	}

	static inline double Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
		double sum = 0;
		for (int j = 0; j < PD_n; j++)
			sum += x[j] * y[j];
		return sum;
	}

	static inline bool Vector_Equal(PT_vector_T x, PT_vector_T y) { // x = y
		for (int j = 0; j < PD_n; j++)
			if (x[j] != y[j])
				return false;
		return true;
	}

	static inline void Vector_i_Copy(PT_vector_i_T fromVector, PT_vector_i_T toVector) { // toVector = fromVector
		for (int j = 0; j < PD_n; j++)
			toVector[j] = fromVector[j];
	}

	static inline void Vector_i_Print(PT_vector_i_T x) {
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << x[j];
		cout << endl;
	}

	static inline void Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector) {
		double norm_x = Vector_Norm(x);
		if (norm_x == 0)
			Vector_Zeroing(likeVector);
		else
			Vector_MultiplyByNumber(x, lengthOfLikeVector / norm_x, likeVector);
	}

	static inline void Vector_MakeMinus_e(PT_vector_T minus_e) {
		for (int j = 0; j < PD_n; j++)
			minus_e[j] = -1;
	}

	static inline void Vector_Median(PT_vector_T x, double length) {
		for (int i = 0; i < PD_n; i++)
			x[i] = length;
	}

	static inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
		for (int j = 0; j < PD_n; j++)
			equalPoint[j] -= minusVector[j];
	}

	static inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] * r;
	}

	static inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
		for (int j = 0; j < PD_n; j++)
			x[j] *= r;
	}

	static inline double Vector_Norm(PT_vector_T x) {
		double norm_x = sqrt(Vector_NormSquare(x));
		return norm_x;
	}

	static inline double Vector_NormSquare(PT_vector_T x) {
		double sum = 0;

		for (int j = 0; j < PD_n; j++) {
			sum += x[j] * x[j];
		}
		return sum;
	}

	static inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
		for (int j = 0; j < PD_n; j++)
			equalVector[j] += plusVector[j];
	}

	static inline void Vector_Print(PT_vector_T x) {
		for (int j = 0; j < PD_n; j++) 
			cout << setw(PP_SETW) << x[j];
		cout << endl;
	}

	static inline void Vector_Random(PT_vector_T x, int seed) {
		srand(seed);
		for (int i = 0; i < PD_n; i++) {
			x[i] = 2 * rand() - RAND_MAX;
		}
	}

	static inline void Vector_Round(PT_vector_T x, double eps_round) {
		double floorValue;
		double fractionalPart;
		double sign;
		double absValue;

		if (eps_round == 0)
			return;

		for (int j = 0; j < PD_n; j++) {
			if (fabs(x[j]) < eps_round) {
				x[j] = 0;
				continue;
			}
			absValue = fabs(x[j]);
			sign = x[j] > 0 ? 1 : -1;
			floorValue = floor(absValue);
			fractionalPart = absValue - floorValue;
			if (1 - fractionalPart < eps_round) {
				x[j] = sign * (floorValue + 1);
				continue;
			}
			if (fractionalPart < eps_round)
				x[j] = sign * floorValue;
		}
	}

	static inline void Vector_SetValue(PT_vector_T x, double v) {  // x = (v,...,v)
		for (int j = 0; j < PD_n; j++) x[j] = v;
	}

	static inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] - y[j];
	}

	static inline void Vector_Zeroing(PT_vector_T x) {  // x = 0
		for (int j = 0; j < PD_n; j++) x[j] = 0;
	}

}

//---------------------------------- Private Functions -------------------------
namespace PF {
	using namespace SF;

	static inline void _Make_A_v(int* basis_v) {
		for (int i = 0; i < PD_n; i++)
			for (int j = 0; j < PD_n; j++)
				PD_A_v[i][j] = PD_A[basis_v[i]][j];
	}

	static inline void _Make_AI_v(double eps_zero, bool* success) {
		double factor;
		int n = PD_n;

		*success = true;

		for (int i = 0; i < n; i++)   // Make a copy of original matrix
			for (int j = 0; j < n; j++) 
					PD_D[i][j] = PD_A_v[i][j];

		for (int i = 0; i < n; i++) { // Make identity matrix 
			for (int j = 0; j < n; j++) 
				PD_AI_v[i][j] = 0;
			PD_AI_v[i][i] = 1;
		}

		for (int j = 0; j < n; j++) {

			if (fabs(PD_D[j][j]) <= eps_zero) { // Search for a non-zero below
				int i_ne_0 = -1;
				for (int i = j + 1; i < n; i++) {
					if (i_ne_0 == -1) {
						if (fabs(PD_D[i][j]) > eps_zero)
							i_ne_0 = i;
					}
					else {
						if (fabs(PD_D[i][j]) > fabs(PD_D[i_ne_0][j]))
							i_ne_0 = i;
					}
				}

				if (i_ne_0 > -1) { // Permute the rows
					for (int l = 0; l < n; l++) {
						double buf;
						buf = PD_D[j][l];
						PD_D[j][l] = PD_D[i_ne_0][l];
						PD_D[i_ne_0][l] = buf;
						buf = PD_AI_v[j][l];
						PD_AI_v[j][l] = PD_AI_v[i_ne_0][l];
						PD_AI_v[i_ne_0][l] = buf;
					}
				}
				else {
					int j_ne_0 = -1; // Search for a non-zero on the right
					for (int lj = j + 1; lj < n; lj++) {
						if (j_ne_0 == -1) {
							if (fabs(PD_D[j][lj]) > eps_zero)
								j_ne_0 = lj;
						}
						else {
							if (fabs(PD_D[j][lj]) > fabs(PD_D[j][j_ne_0]))
								j_ne_0 = lj;
						}
					}

					if (j_ne_0 == -1) {
						/* Debug _Make_AI_v**
						cout << "Row " << j << endl;
						for (int jj = j + 1; jj < n; jj++)
							cout << PD_D[j][jj] << "\t";
						cout <<"\nColumn " << j << endl;
						for (int ii = j + 1; ii < n; ii++)
							cout << PD_D[ii][j] << "\t";
						cout << endl;
						/* End Debug _Make_AI_v*/
						*success = false;
						return;
					}

					for (int i = 0; i < n; i++) { // Permute the columns
						double buf;
						buf = PD_D[i][j];
						PD_D[i][j] = PD_D[i][j_ne_0];
						PD_D[i][j_ne_0] = buf;
						buf = PD_AI_v[i][j];
						PD_AI_v[i][j] = PD_AI_v[i][j_ne_0];
						PD_AI_v[i][j_ne_0] = buf;
					}
				}
			}

			factor = PD_D[j][j];
			for (int l = 0; l < n; l++) { // make 1
				PD_D[j][l] = PD_D[j][l] / factor;
				PD_AI_v[j][l] = PD_AI_v[j][l] / factor;
			}

			for (int i = j + 1; i < n; i++) {
				factor = PD_D[i][j];
				for (int l = 0; l < n; l++) {
					PD_D[i][l] = PD_D[i][l] - PD_D[j][l] * factor;
					PD_AI_v[i][l] = PD_AI_v[i][l] - PD_AI_v[j][l] * factor;
				}
			}
		}

		for (int i = n - 1; i > 0; i--)
			for (int k = i - 1; k >= 0; k--) {
				factor = PD_D[k][i];
				for (int j = 0; j < n; j++) {
					PD_D[k][j] = PD_D[k][j] - PD_D[i][j] * factor;
					PD_AI_v[k][j] = PD_AI_v[k][j] - PD_AI_v[i][j] * factor;
				}
			}
	}

	static inline void _Make_cAI_v(PT_vector_T cAI_v){
		for (int j = 0; j < PD_n; j++) {
			cAI_v[j] = 0;
			for (int i = 0; i < PD_n; i++)
				cAI_v[j] = cAI_v[j] + PD_c[i] * PD_AI_v[i][j];
		}
	}

	static inline void _Make_u(int* basis_v, PT_vector_T cAI_v, double eps_zero) {
		Column_Zeroing(PD_u);
		for (int i = 0; i < PD_n; i++)
			if (fabs(cAI_v[i]) > eps_zero)
				PD_u[basis_v[i]] = cAI_v[i];
			else
				PD_u[basis_v[i]] = 0;
	}

	static inline void MakeBasis_v(PT_vector_T v, PT_vector_i_T basis_v) {
		MakeNeHyperplaneList_v(v, PD_neHyperplanes, PD_mneh, PD_neHyperplanes_v, &PD_mneh_v, PP_EPS_ON_HYPERPLANE);

		PD_m_v = PD_meq_basic + PD_mneh_v;
		assert(PD_m_v <= PP_MM);

		for (int i = PD_meq_basic; i < PD_m_v; i++)
			basis_v[i] = PD_neHyperplanes_v[i - PD_meq_basic];

		List_i_Basis(basis_v, &PD_m_v, PP_EPS_ZERO);
	}

	static inline bool OptimumIsFound(double eps_zero) {
		for (int i = 0; i < PD_m; i++) {
			if (PD_isEquation[i])
				continue;
			if (PD_u[i] < -eps_zero)
				return false;
		}
		return true;
	}

	static inline bool PointIsVertex(PT_vector_T v, double eps_on_hyperplane) {
		int rank;
		int basis_v[PP_MM];

		MakeNeHyperplaneList_v(v, PD_neHyperplanes, PD_mneh, basis_v, &PD_m_v, PP_EPS_ON_HYPERPLANE);
		assert(PD_m_v <= PP_MM);
		Matrix__Rank(basis_v, PD_m_v, PP_EPS_ZERO, &rank);
		assert(PD_supportSubspaceDim > 0 && PD_supportSubspaceDim <= PD_n);
		if (rank >= PD_supportSubspaceDim)
			return true;
		else
			return false;
	}

	static inline void PreparationForIteration(PT_vector_i_T basis_v) {
		PT_vector_T cAI_v;

		/* PreparationForIteration */
		#ifdef PP_DEBUG
		int rank;
		Matrix__Rank(basis_v, PD_n, PP_EPS_ZERO, &rank);
		assert(rank == PD_n);
		#endif PP_DEBUG /**/

		_Make_A_v(basis_v);

		bool success;
		_Make_AI_v(PP_EPS_ZERO, &success);
		assert(success);
		/* PreparationForIteration **
		#ifdef PP_DEBUG
		for (int i = 0; i < PD_n; i++)
			for (int j = 0; j < PD_n; j++) {
				double s = 0;
				for (int k = 0; k < PD_n; k++)
					s = s + PD_A_v[i][k] * PD_AI_v[k][j];
				if (i == j)
					assert(fabs(s - 1) <= PP_EPS_ZERO);
				else
					assert(fabs(s) <= PP_EPS_ZERO);
			}
		#endif PP_DEBUG /**/

		_Make_cAI_v(cAI_v);
		_Make_u(basis_v, cAI_v, PP_EPS_ZERO);
	}
}