#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <math.h> 

/*
	Constant Parameters
*/
const int NS = 9;
const int mu_max = 5;
const double PI_VAL = 3.14159265358979323846;

typedef std::vector<int> vi;
typedef std::vector<double> vd;
typedef std::vector<vi> vvi;
typedef std::vector<vd> vvd;
typedef std::vector<vvd> vvvd;

/*
	comp_vector: Compare two vectors in lexicographic order + level of precesion
		Example:
		We compare the level of precesion first by adding up all the elements.
		
			(1) a = (1, 2, 2, 2, 3) and b = (3, 1, 1, 1, 1)
				then sum(a) = 10 and sum(b) = 7, so b precedes a.
			(2) a = (2, 2, 1, 1, 1, 2) and b = (2, 2, 2, 1, 1, 1)
				sum(a) = 9 = sum(b), then we compare the lexicographic order.
				a precedes b because 1st and 2nd elements are same,
				but third element of a is less than that of b.
*/
template <typename T> bool comp_vectors(const std::vector<T> &invec1, 
	const std::vector<T> &invec2) {
	
	T sum1 = std::accumulate(invec1.begin(), invec1.end(), 0);
	T sum2 = std::accumulate(invec2.begin(), invec2.end(), 0);

	// If level of precesion is the same, compare two vectors in lexicographic order
	if (sum1 == sum2) {
		int ii = 0;
		bool tempbool = true;
		bool tempbool2 = true;

		while (tempbool2 && ii < NS) {
			tempbool2 = invec1[ii] == invec2[ii];
			tempbool = invec1[ii] < invec2[ii];
			ii += 1;
		}

		return tempbool;
	}
	// If level of precesion is different, then the one with less precesion comes first
	else {
		return sum1 < sum2;
	}
}

template <typename T> bool equiv_vectors(const std::vector<T> &in_vec1, const std::vector<T> &in_vec2) {
	bool tempbool = true;
	for (int ii = 0; ii < NS; ii++) {
		tempbool = tempbool & (in_vec1[ii] == in_vec2[ii]);
	}

	return tempbool;
}

int count_points(vi in_vec) {
	int out = 1;

	for (int ii = 0; ii < NS; ii++) {
		switch (in_vec[ii]) {
		case 1: out *= 1; break;
		case 2: out *= 2; break;
		default: 
			for (int ii = 0; ii < in_vec[ii] - 2; ii++) {
				out *= 2;
			}
			break;
		}
	}

	return out;
}

/*
	cartesian_1t1: Outer product between vector & matrix

	Warning: function overloaded with the one below (vector & vector)
*/
vvi cartesian_1t1(const vi& in_vec, const vvi& in_mat) {
	int N1 = static_cast<int>(in_vec.size());
	int N2 = static_cast<int>(in_mat.size());
	int C2 = static_cast<int>(in_mat[0].size());

	vvi RetMat = vvi(N1 * N2, vi(C2 + 1));

	for (int ii = 0; ii < N1; ii++) {
		for (int jj = 0; jj < N2; jj++) {
			RetMat[ii * N2 + jj][0] = in_vec[ii];

			for (int cc = 0; cc < C2; cc++) {
				RetMat[ii * N2 + jj][cc + 1] = in_mat[jj][cc];
			}
		}
	}

	return RetMat;
}

/*
	cartesian_1t1: Cartesian product between vector & vector (vector of vectors)

	Warning: function overloaded with the one above (vector & matrix)
*/
vvi cartesian_1t1(const vi &in_vec, const vi &in_mat) {
	int N1 = static_cast<int>(in_vec.size());
	int N2 = static_cast<int>(in_mat.size());

	vvi RetMat = vvi(N1 * N2, vi(2));

	for (int ii = 0; ii < N1; ii++) {
		for (int jj = 0; jj < N2; jj++) {
			RetMat[ii * N2 + jj][0] = in_vec[ii];
			RetMat[ii * N2 + jj][1] = in_mat[jj];
		}
	}

	return RetMat;
}

/*
	cartesian_product: Cartesian product between sequence of vectors
*/
vvi cartesian_product(const vvi& invec) {
	int NIter = static_cast<int>(invec.size());

	vvi RetMat = cartesian_1t1(invec[NIter - 2], invec[NIter - 1]);

	if (NIter > 2) {
		for (int ii = 3; ii <= NIter; ii++) {
			RetMat = cartesian_1t1(invec[NIter - ii], RetMat);
		}
	}

	return RetMat;
}

/*
	MakeIndexes: make indexes when given precesion levels
*/
vvi MakeIndexes(const vi& in_levels) {
	vvi RetMat(NS);

	for (int ii = 0; ii < NS; ii++) {
		switch (in_levels[ii]) {
		case 1: RetMat[ii] = { 1 }; break;
		case 2: RetMat[ii] = { 2, 3 }; break;
		default:
			int tempval = 1;

			for (int jj = 0; jj < in_levels[ii] - 2; jj++) {
				tempval *= 2;
			}

			int sta_val = tempval + 2;	// sta_val = 2^(in_levels[ii] - 2) + 2;
			int end_val = tempval*2 + 1;	// sta_val = 2^(in_levels[ii] - 1) + 1;
			int NN = end_val - sta_val + 1;

			RetMat[ii] = vi(NN, 0);

			for (int jj = 0; jj < NN; jj++) {
				RetMat[ii][jj] = sta_val + jj;
			}
			break;
		}
	}

	return RetMat;
}

vvi Smolyak_Elem_Isotrop() {
	vvi Smol_rule = vvi(1, vi(NS,1));
	vvi prev_incr = Smol_rule;

	int mm_end = 1;
	int mm_sta = 0;
	size_t mm_new = 0;

	vvi::iterator it;

	for (int ii = 1; ii < mu_max+1; ii++) {
		vvi temp_incr((mm_end - mm_sta)*NS, vi(NS,0));

		// Increment from previous iteration for each dimension
		for (int jj = 0; jj < mm_end - mm_sta; jj++) {
			for (int kk = 0; kk < NS; kk++) {
				temp_incr[jj*NS + kk] = prev_incr[jj];
				temp_incr[jj*NS + kk][NS - kk - 1] += 1;
			}
		}
		
		// Remove duplicates
		sort(temp_incr.begin(), temp_incr.end(), comp_vectors<int>);
		temp_incr.erase(std::unique(temp_incr.begin(), temp_incr.end()), temp_incr.end());
		prev_incr = temp_incr;

		mm_sta = mm_end;
		mm_end += static_cast<int>(temp_incr.size());

		// Copy to Smol_rule
		Smol_rule.resize(mm_end);

		for (int jj = mm_sta; jj < mm_end; jj++) {
			Smol_rule[jj] = temp_incr[jj-mm_sta];
		}
	}

	int num_points = static_cast<int>(Smol_rule.size());
	mm_sta = 0;
	mm_end = 0;
	vvi Smolyak_elem_iso;
	vvi Indexes(NS);
	vvi zz;

	// Make cartesian product here
	for (int ii = 0; ii < num_points; ii++) {
		Indexes = MakeIndexes(Smol_rule[ii]);
		zz = cartesian_product(Indexes);

		mm_sta = mm_end;
		mm_end += static_cast<int>(zz.size());

		Smolyak_elem_iso.resize(mm_end);

		for (int jj = mm_sta; jj < mm_end; jj++) {
			Smolyak_elem_iso[jj] = zz[jj - mm_sta];
		}
	}

	return Smolyak_elem_iso;
}

// remove_duplicates: remove duplicates from in_vec and return vector with unique elemnts
//	- order is also retained
vd remove_duplicates(const vd& in_vec) {
	int vec_size = static_cast<int>(in_vec.size());
	int countval = 1, comp_idx;

	vd out_vec(vec_size);
	double smallnum = 1e-12;
	
	bool temp_bool;

	out_vec[0] = in_vec[0];

	for (int ii = 1; ii < vec_size; ii++) {
		comp_idx = 0;
		temp_bool = true;

		while (temp_bool && comp_idx < ii) {
			temp_bool = (std::abs(in_vec[ii] - in_vec[comp_idx]) > smallnum);
			comp_idx += 1;
		}

		if (temp_bool) {
			out_vec[countval] = in_vec[ii];
			countval += 1;
		}
	}

	out_vec.resize(countval);

	return out_vec;
}

vvd Smolyak_Grid(const vvi &Smol_elem) {
	int i_max = mu_max + 1, m_i = 1, mm_sta, mm_end = 1;
	double temp1, temp2;

	vd points_ld{ 0 };

	for (int ii = 2; ii < i_max+1; ii++) {
		int tempval = 1;

		for (int jj = 0; jj < ii - 1; jj++) {
			tempval *= 2;
		}

		m_i = tempval + 1;

		vd extrem_Cheb_1d(m_i);
		
		for (int jj = 0; jj < m_i; jj++) {
			temp1 = static_cast<double>(jj);
			temp2 = static_cast<double>(m_i - 1.0);

			extrem_Cheb_1d[jj] = -cos(PI_VAL * temp1 / temp2);
			
			if (abs(extrem_Cheb_1d[jj]) < 1e-12) {
				extrem_Cheb_1d[jj] = 0.0;
			}
			else if (1 - extrem_Cheb_1d[jj] < 1e-12) {
				extrem_Cheb_1d[jj] = 1.0;
			}
			else if (1 + extrem_Cheb_1d[jj] < 1e-12) {
				extrem_Cheb_1d[jj] = -1.0;
			}	
		}

		mm_sta = mm_end;
		mm_end += m_i;

		points_ld.resize(mm_end);

		for (int jj = mm_sta; jj < mm_end; jj++) {
			points_ld[jj] = extrem_Cheb_1d[jj - mm_sta];
		}

		int qweqweqwe = 0;
	}

	points_ld = remove_duplicates(points_ld);

	int mat_len = static_cast<int>(Smol_elem.size());
	vvd Smol_grid(mat_len, vd(NS, 0));

	for (int jp = 0; jp < mat_len; jp++) {
		for (int jd = 0; jd < NS; jd++) {
			// Matlab uses 1-based indexing so I take minus one.
			Smol_grid[jp][jd] = points_ld[Smol_elem[jp][jd] - 1];
		}
	}

	return Smol_grid;
}

vd Smolyak_Polynomial(const vvd &points, const vvi &Smol_elem) {
	int i_max = mu_max + 1;

	int m_i_max = 1;
	if (i_max > 1) {
		for (int ii = 0; ii < i_max - 1; ii++) {
			m_i_max *= 2;
		}

		m_i_max += 1;	// m_i_max = 2^(i_max - 1) - 1
	}

	long int numb_pts = static_cast<long int>(points.size());
	int idx_phi1;
	int idx_phikk;

	long int mat_size = numb_pts * NS * m_i_max;

	vd phi(mat_size, 1.0);

	for (int ii = 0; ii < numb_pts; ii++) {
		for (int jj = 0; jj < NS; jj++) {
			idx_phi1 = ii * NS * m_i_max + jj * m_i_max + 1;

			phi[idx_phi1] = points[ii][jj];

			for (int kk = 2; kk < m_i_max; kk++) {
				idx_phikk = ii * NS * m_i_max + jj * m_i_max + kk;

				phi[idx_phikk] = 2.0 * phi[idx_phi1] * phi[idx_phikk - 1] - phi[idx_phikk - 2];
			}
		}
	}

	int numb_terms = static_cast<int>(Smol_elem.size());
	mat_size = numb_pts * numb_terms;

	std::cout << mat_size << std::endl;

	vd Smol_bases(mat_size, 1.0);

	/*
	double* Smol_bases = new double[(__int64) numb_pts * numb_terms];

	for (int ii = 0; ii < numb_pts; ii++) {
		for (int jj = 0; jj < numb_terms; jj++) {
			Smol_bases[ii * numb_terms + jj] = 1.0;
		}
	}
	*/

	for (int ii = 0; ii < numb_pts; ii++) {
		for (int jt = 0; jt < numb_terms; jt++) {
			for (int jd = 0; jd < NS; jd++) {
				idx_phi1 = ii * NS * m_i_max + jd * m_i_max + Smol_elem[jt][jd] - 1;

				// Matlab uses 1-based indexing so I take minus one.
				Smol_bases[ii*numb_terms+jt] *= phi[idx_phi1];
			}
		}
	}

	return Smol_bases;
}