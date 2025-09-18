#ifndef RINTINTIN_EIGEN_H
#define RINTINTIN_EIGEN_H
#include "rintintin_scratch.h"

rintintin_eigen rintintin_compute_eigen_m3(rintintin_mat3x3 *I);

typedef enum rintintin_matrix_classification
{
	RINTINTIN_POSITIVE		= 1 << 0,
	RINTINTIN_NEGATIVE		= 1 << 1,
	RINTINTIN_DEFINITE		= 1 << 2,
	RINTINTIN_SEMI_DEFINITE	= 1 << 3,
} rintintin_matrix_classification;

rintintin_matrix_classification rintintin_sym_mat3_classify(rintintin_symmetric_mat3 const* it);

#endif // RINTINTIN_EIGEN_H
