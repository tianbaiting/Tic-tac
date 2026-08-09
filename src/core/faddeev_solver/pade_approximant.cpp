#include "pade_approximant.h"

#include "utils/matrix_routines.h"

#include <cmath>
#include <vector>

cdouble pade_approximant(const cdouble* a_coeff_array,
                         std::size_t N,
                         std::size_t M,
                         cdouble z){
	const std::size_t dim = M + 1;
	std::vector<cdouble> P_array(dim * dim);
	std::vector<cdouble> Q_array(dim * dim);

	for (std::size_t row_idx = 0; row_idx < M; row_idx++){
		for (std::size_t col_idx = 0; col_idx < dim; col_idx++){
			const cdouble value = a_coeff_array[N - M + 1 + row_idx + col_idx];
			P_array[row_idx * dim + col_idx] = value;
			Q_array[row_idx * dim + col_idx] = value;
		}
	}

	for (std::size_t col_idx = 0; col_idx < dim; col_idx++){
		Q_array[M * dim + col_idx] = std::pow(z, M - col_idx);
		P_array[M * dim + col_idx] = 0.0;
		for (std::size_t j = M - col_idx; j < N + 1; j++){
			P_array[M * dim + col_idx] +=
				a_coeff_array[j - (M - col_idx)] * std::pow(z, j);
		}
	}

	const cdouble P_det = determinant(P_array.data(), (int)dim);
	const cdouble Q_det = determinant(Q_array.data(), (int)dim);
	return P_det / Q_det;
}
