#include "rintintin_mesh.h"
#define USING_NAMESPACE_RINTINTIN 1
#define RINTINTIN_USE_MATH_H 0
#include "rintintin_eigen.h"
#include "rintintin_math.h"

void rintintin_compute_rotation_3x3(rintintin_mat3x3 * R, rintintin_eigen const* eig)
{
    if (!R || !eig) return;
    // Treat eigenvectors as columns, directly
    dvec3 c0 = rintintin_vec3_normalize(&eig->vectors[0]);
    dvec3 c1 = rintintin_vec3_normalize(&eig->vectors[1]);
    // Gram-Schmidt orthogonalize
    double dot01 = rintintin_vec3_dot(&c0, &c1);
    dvec3 tmp1 = {c1.x - dot01 * c0.x,
                 c1.y - dot01 * c0.y,
                 c1.z - dot01 * c0.z};
    c1 = rintintin_vec3_normalize(&tmp1);
    // Third column from cross → right-handed basis
    dvec3 c2 = rintintin_vec3_cross(&c0, &c1);
    R->vec[0] = c0;
    R->vec[1] = c1;
    R->vec[2] = c2;

    // Determinant check AFTER orthonormalization
    double det =
          c0.x * (c1.y * c2.z - c1.z * c2.y)
        - c0.y * (c1.x * c2.z - c1.z * c2.x)
        + c0.z * (c1.x * c2.y - c1.y * c2.x);
    if (det < 0.0) {
        // Flip third column
        R->vec[2].x = -R->vec[2].x;
        R->vec[2].y = -R->vec[2].y;
        R->vec[2].z = -R->vec[2].z;
    }

}

rintintin_vec4 rintintin_quat_from_3x3(rintintin_mat3x3 const * in)
{
	const double (*m)[3] = in->dbl;
	
	double fourXSquaredMinus1 = m[0][0] - m[1][1] - m[2][2];
	double fourYSquaredMinus1 = m[1][1] - m[0][0] - m[2][2];
	double fourZSquaredMinus1 = m[2][2] - m[0][0] - m[1][1];
	double fourWSquaredMinus1 = m[0][0] + m[1][1] + m[2][2];

	int biggestIndex = 0;
	double fourBiggestSquaredMinus1 = fourWSquaredMinus1;
	if(fourXSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourXSquaredMinus1;
		biggestIndex = 1;
	}
	if(fourYSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourYSquaredMinus1;
		biggestIndex = 2;
	}
	if(fourZSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourZSquaredMinus1;
		biggestIndex = 3;
	}

	double biggestVal = sqrt(fourBiggestSquaredMinus1 + 1.0) * 0.5;
	double mult = 0.25 / biggestVal;

	switch(biggestIndex)
	{
	case 0:
		return (rintintin_vec4){.w=biggestVal, .x=(m[1][2] - m[2][1]) * mult, .y=(m[2][0] - m[0][2]) * mult, .z=(m[0][1] - m[1][0]) * mult};
	case 1:
		return (rintintin_vec4){.w=(m[1][2] - m[2][1]) * mult, .x=biggestVal, .y=(m[0][1] + m[1][0]) * mult, .z=(m[2][0] + m[0][2]) * mult};
	case 2:
		return (rintintin_vec4){.w=(m[2][0] - m[0][2]) * mult, .x=(m[0][1] + m[1][0]) * mult, .y=biggestVal, .z=(m[1][2] + m[2][1]) * mult};
	case 3:
		return (rintintin_vec4){.w=(m[0][1] - m[1][0]) * mult, .x=(m[2][0] + m[0][2]) * mult, .y=(m[1][2] + m[2][1]) * mult, .z=biggestVal};
	default:
		return (rintintin_vec4){0, 0, 0, 1};
	}
	
	
}

rintintin_vec4 rintintin_compute_rotation_quat(rintintin_eigen const* eigen)
{
	rintintin_mat3x3 rotMatrix;
	rintintin_compute_rotation_3x3(&rotMatrix, eigen);
	return rintintin_quat_from_3x3(&rotMatrix);
}

rintintin_transform rintintin_transform_from_4x4(double (*matrix)[4][4])
{
   rintintin_transform result;
   
   // Extract translation from the last column
   result.translation.x = (*matrix)[0][3];
   result.translation.y = (*matrix)[1][3];
   result.translation.z = (*matrix)[2][3];
   
   // Extract the 3x3 upper-left submatrix (rotation + scale)
   double upper3x3[3][3];
   for (int i = 0; i < 3; i++) {
       for (int j = 0; j < 3; j++) {
           upper3x3[i][j] = (*matrix)[i][j];
       }
   }
   
   // Extract scale by computing column lengths
   result.scaling.x = sqrt(upper3x3[0][0] * upper3x3[0][0] + 
                          upper3x3[1][0] * upper3x3[1][0] + 
                          upper3x3[2][0] * upper3x3[2][0]);
   
   result.scaling.y = sqrt(upper3x3[0][1] * upper3x3[0][1] + 
                          upper3x3[1][1] * upper3x3[1][1] + 
                          upper3x3[2][1] * upper3x3[2][1]);
   
   result.scaling.z = sqrt(upper3x3[0][2] * upper3x3[0][2] + 
                          upper3x3[1][2] * upper3x3[1][2] + 
                          upper3x3[2][2] * upper3x3[2][2]);
   
   // Remove scale to get pure rotation matrix
   rintintin_mat3x3 mat;
   double (*rotation3x3)[3] = mat.dbl;
   for (int i = 0; i < 3; i++) {
       rotation3x3[i][0] = upper3x3[i][0] / result.scaling.x;
       rotation3x3[i][1] = upper3x3[i][1] / result.scaling.y;
       rotation3x3[i][2] = upper3x3[i][2] / result.scaling.z;
   }
   
   // Handle potential negative determinant (reflection)
   // If determinant is negative, we have a reflection - flip one axis
   double det = rotation3x3[0][0] * (rotation3x3[1][1] * rotation3x3[2][2] - rotation3x3[1][2] * rotation3x3[2][1]) -
                rotation3x3[0][1] * (rotation3x3[1][0] * rotation3x3[2][2] - rotation3x3[1][2] * rotation3x3[2][0]) +
                rotation3x3[0][2] * (rotation3x3[1][0] * rotation3x3[2][1] - rotation3x3[1][1] * rotation3x3[2][0]);
   
   if (det < 0) {
       // Flip the X axis and negate its scale
       result.scaling.x = -result.scaling.x;
       for (int i = 0; i < 3; i++) {
           rotation3x3[i][0] = -rotation3x3[i][0];
       }
   }
   
   // Convert rotation matrix to quaternion
   result.rotation = rintintin_quat_from_3x3(&mat);
   
   return result;
}


// Sort descending: values[0] is the largest eigenvalue, vectors[0] its eigenvector.
// For a covariance matrix this means vectors[0] points along the principal (longest)
// axis and vectors[2] along the minor (shortest) axis. Both call sites rely on this
// ordering so that the rotation built from vectors[0..2] as columns puts the major
// axis on local X, matching scale order (large, mid, small).
static rintintin_eigen rintintin_sort_eigen(rintintin_eigen * in)
{
	int _max = in->values[0] >= in->values[1] ? 0 : 1;
	int _min = 1 - _max;

	if (in->values[2] > in->values[_max]) _max = 2;
	else if (in->values[2] < in->values[_min]) _min = 2;
	//the three indices must sum to 0+1+2=3.
	int _mid = 3 - _min - _max;

	rintintin_eigen r = {
		.values={in->values[_max], in->values[_mid], in->values[_min]},
		.vectors={in->vectors[_max], in->vectors[_mid], in->vectors[_min]},
	};

	return r;
}

static void makeConsistent(rintintin_eigen* result) {
      // Choose the eigenvector component with largest absolute value
      // and make it positive (this matches Eigen's typical behavior)
      for (int i = 0; i < 3; i++) {
          dvec3 * v = &result->vectors[i];

          // Find component with largest absolute value
          int max_idx = 0;
          if (fabs(v->y) > fabs(v->x)) max_idx = 1;
          if (fabs(v->z) > fabs((&v->x)[max_idx])) max_idx = 2;

          // Make that component positive
          if ((&v->x)[max_idx] < 0) {
              v->x = -v->x;
              v->y = -v->y;
              v->z = -v->z;
          }
      }
  }
  
rintintin_eigen rintintin_compute_eigen(rintintin_symmetric_mat3 const*I) {
    dmat3 A = {.m={
		.x={I->xx, I->xy, I->xz},
		.y={I->xy, I->yy, I->yz},
		.z={I->xz, I->yz, I->zz}    
    }};
    
    return rintintin_compute_eigen_m3(&A);
}


rintintin_eigen rintintin_compute_eigen_m3(rintintin_mat3x3 const* I) {
    // Pure function: work on a local copy so the caller's matrix is untouched.
    double A[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            A[i][j] = I->dbl[i][j];

    double V[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    const int maxIterations = 50;

    // Frobenius-ish norm of the matrix, used as the only scale reference.
    double matrixScale = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            matrixScale += fabs(A[i][j]);
        }
    }

    if (matrixScale == 0.0) {
        return (rintintin_eigen) {
            .vectors = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
            .values  = {0, 0, 0}
        };
    }

    matrixScale /= 9.0; // average absolute value of matrix elements

    // Normalise the matrix to ~unit magnitude before iterating. This keeps the
    // off-diagonal convergence test scale-free: with elements ~1 the threshold
    // for "this off-diagonal is effectively zero" is just machine epsilon.
    // Eigenvalues are de-normalised at the end.
    double scale_factor = 1.0 / matrixScale;
    for (int i = 0; i < 9; ++i)
        A[i/3][i%3] *= scale_factor;

    // After normalisation A has |elements| ~1, so a relative tolerance
    // of ~1e-12 is "off-diagonal indistinguishable from zero in double
    // precision". The previous adaptive formula collapsed to ~1e-6, which
    // let Jacobi early-out in the middle of resolving a near-degenerate
    // 2D eigenspace (the rod/disk case) and produced a 90deg-off rotation.
    const double tolerance = 1e-12;

    for (int iter = 0; iter < maxIterations; iter++) {
        int p = 0, q = 1;
        double maxOffDiag = 0.0;

        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                double val = fabs(A[i][j]);
                if (val > maxOffDiag) {
                    maxOffDiag = val;
                    p = i;
                    q = j;
                }
            }
        }

        if (maxOffDiag < tolerance) {
            break;
        }

        double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = (tau >= 0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;

        double App = A[p][p];
        double Aqq = A[q][q];
        double Apq = A[p][q];

        A[p][p] = c * c * App - 2.0 * s * c * Apq + s * s * Aqq;
        A[q][q] = s * s * App + 2.0 * s * c * Apq + c * c * Aqq;
        A[p][q] = A[q][p] = 0.0;

        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                double Arp = A[r][p];
                double Arq = A[r][q];
                A[r][p] = A[p][r] = c * Arp - s * Arq;
                A[r][q] = A[q][r] = s * Arp + c * Arq;
            }
        }

        for (int i = 0; i < 3; i++) {
            double Vip = V[i][p];
            double Viq = V[i][q];
            V[i][p] = c * Vip - s * Viq;
            V[i][q] = s * Vip + c * Viq;
        }
    }

    rintintin_eigen result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = A[i][i] / scale_factor;
        result.vectors[i].x = V[0][i];
        result.vectors[i].y = V[1][i];
        result.vectors[i].z = V[2][i];
    }

    makeConsistent(&result);
    return rintintin_sort_eigen(&result);
}

rintintin_error_code rintintin_estimate_shapes(rintintin_inertia_estimation * dst, rintintin_metrics const* src, uint64_t no_items)
{
	if(dst == 0 || src == 0)
		return RINTINTIN_ERROR_NULL_POINTER;

    // Degenerate-volume fallback: small but nonzero so downstream math (e.g.
    // matrix inverse for OBB orientation drawing) does not divide by zero.
    const rintintin_vec3 degenerate_scale = {0.0001, 0.0001, 0.0001};

    for (uint64_t i = 0; i < no_items; i++) {
        // Eigendecompose the COVARIANCE (second_moment about centroid) directly.
        // Inertia and covariance share eigenvectors but their eigenvalue ordering
        // is inverted (I = tr(C)*Id - C), so using covariance lets the sort produce
        // values[0] = largest physical extent, matching rotation column order.
        rintintin_eigen eigen = rintintin_compute_eigen(&src[i].second_moment);

        dst[i].translation = src[i].centroid;
        dst[i].rotation = rintintin_compute_rotation_quat(&eigen);

        double volume = src[i].volume;

        // Degenerate input: no volume, or covariance lost to floating point.
        // Cannot recover semi-axes; emit a tiny non-zero scale.
        if (!(volume > 1e-30) || !(fabs(eigen.values[0]) > 1e-30)) {
            dst[i].scale = degenerate_scale;
            continue;
        }

        // ELLIPSOID fit: for a solid ellipsoid of uniform density with semi-axes
        // (a,b,c) and volume V, the covariance eigenvalues are lambda_i = V*a_i^2 / 5.
        // Inverted: a_i = sqrt(5 * lambda_i / V). Descending eigen sort means
        // values[0] is the major semi-axis, which aligns with local X (rotation
        // column 0), so scale order naturally matches rotation column order.
        double k = 5.0 / volume;
        double sx = sqrt(fabs(eigen.values[0]) * k);
        double sy = sqrt(fabs(eigen.values[1]) * k);
        double sz = sqrt(fabs(eigen.values[2]) * k);

        dst[i].scale = (rintintin_vec3){sx, sy, sz};
    }

    return RINTINTIN_SUCCESS;
}

/**
 * @brief Compute 3x3 inverse matrix
 */
int rintintin_compute_inverse_3x3(rintintin_mat3x3 * dst, rintintin_mat3x3 const* src)
{
	if(src == 0L || dst == 0L) return -2;

    const double (*m)[3] = src->dbl;
    double (*inv)[3] = dst->dbl;
    
    // Compute determinant
    double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    
    // Check for singular matrix
    if (fabs(det) < 1e-100) {
        // Handle singular matrix - could set to identity, NaN, or error
        // For now, set to zero matrix to indicate failure
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                inv[i][j] = 0.0;
            }
        }
        return -1;
    }
    
    double inv_det = 1.0 / det;
    
    // Compute inverse using cofactor matrix
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det;
    
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;
    
    return 0;
}



rintintin_matrix_classification rintintin_sym_mat3_classify(rintintin_symmetric_mat3 const* it) {
    if (!it) {
        return 0; // Invalid input
    }
    
    // Use Sylvester's criterion - check signs of leading principal minors
    // For 3x3 symmetric matrix:
    // M1 = xx
    // M2 = det([[xx, xy], [xy, yy]]) = xx*yy - xy*xy
    // M3 = det(full matrix)
    
    const double tolerance = 1e-12;
    
    // First leading principal minor
    double M1 = it->xx;
    
    // Second leading principal minor  
    double M2 = it->xx * it->yy - it->xy * it->xy;
    
    // Third leading principal minor (full determinant)
    double M3 = it->xx * (it->yy * it->zz - it->yz * it->yz)
              - it->xy * (it->xy * it->zz - it->yz * it->xz)
              + it->xz * (it->xy * it->yz - it->yy * it->xz);
    
    // Check for positive definite: all leading principal minors > 0
    if (M1 > tolerance && M2 > tolerance && M3 > tolerance) {
        return RINTINTIN_POSITIVE | RINTINTIN_DEFINITE;
    }
    
    // Check for negative definite: alternating signs starting with M1 < 0
    if (M1 < -tolerance && M2 > tolerance && M3 < -tolerance) {
        return RINTINTIN_NEGATIVE | RINTINTIN_DEFINITE;
    }
    
    // For semi-definite cases, we need to check if all eigenvalues have the same sign
    // This requires checking all principal minors of all sizes
    
    // All 1x1 principal minors
    double minor_xx = it->xx;
    double minor_yy = it->yy;
    double minor_zz = it->zz;
    
    // All 2x2 principal minors
    double minor_xy = it->xx * it->yy - it->xy * it->xy;  // top-left 2x2
    double minor_xz = it->xx * it->zz - it->xz * it->xz;  // (0,0),(2,2) 2x2
    double minor_yz = it->yy * it->zz - it->yz * it->yz;  // bottom-right 2x2
    
    // Check for positive semi-definite
    // All principal minors >= 0, at least one = 0
    if (minor_xx >= -tolerance && minor_yy >= -tolerance && minor_zz >= -tolerance &&
        minor_xy >= -tolerance && minor_xz >= -tolerance && minor_yz >= -tolerance &&
        M3 >= -tolerance) {
        
        // Check if at least one minor is zero (making it semi-definite, not definite)
        if (fabs(minor_xx) <= tolerance || fabs(minor_yy) <= tolerance || fabs(minor_zz) <= tolerance ||
            fabs(minor_xy) <= tolerance || fabs(minor_xz) <= tolerance || fabs(minor_yz) <= tolerance ||
            fabs(M3) <= tolerance) {
            return RINTINTIN_POSITIVE | RINTINTIN_SEMI_DEFINITE;
        }
    }
    
    // Check for negative semi-definite
    // This is trickier - we need (-1)^k * M_k >= 0 for all principal minors M_k of size k
    // For negative semi-definite: (-A) should be positive semi-definite
    if (-minor_xx >= -tolerance && -minor_yy >= -tolerance && -minor_zz >= -tolerance &&
        minor_xy >= -tolerance && minor_xz >= -tolerance && minor_yz >= -tolerance &&
        -M3 >= -tolerance) {
        
        // Check if at least one condition makes it semi- rather than definite
        if (fabs(minor_xx) <= tolerance || fabs(minor_yy) <= tolerance || fabs(minor_zz) <= tolerance ||
            fabs(minor_xy) <= tolerance || fabs(minor_xz) <= tolerance || fabs(minor_yz) <= tolerance ||
            fabs(M3) <= tolerance) {
            return RINTINTIN_NEGATIVE | RINTINTIN_SEMI_DEFINITE;
        }
    }
    
    // If none of the above, matrix is indefinite
    return 0;
}


// Argmax skin influence: returns the joint index this vertex is most strongly
// bound to, or -1 if no positive-weight joint exists. Used by both the
// covariance accumulation pass and the extents (min/max) pass so that the OBB
// rotation and extents are computed on the SAME point cloud (eliminating the
// rotation error caused by skinning weight bleed between chained bones).
static int rintintin_obb_argmax_joint(struct rintintin_vertex const* vert)
{
	int dom = -1;
	double dom_w = 0.0;
	for (int i = 0; i < 4; ++i) {
		double w = (&vert->weight.x)[i];
		if (vert->joint[i] >= 0 && w > dom_w) {
			dom_w = w;
			dom = i;
		}
	}
	return dom < 0 ? -1 : vert->joint[dom];
}

// Per-joint accumulator for the argmax-cluster PCA pass.
//
// Two phases, two views over the same 80-byte slot:
//   Phase A (pass 1, accumulating raw moments):
//     count + sum_xyz + raw outer-product sums (sxx..szz).
//   Phase B (post-prep, used by extents pass + refinement):
//     centroid_xyz (aliases sum_xyz, written in the prep stage) plus the
//     spent slots (was count + sxx..szz) repurposed as best_min_xyz / best_max_xyz
//     scratch for the refinement loop's save-and-revert pattern. count's slot
//     is unused in phase B (8 spare bytes).
typedef struct rintintin_obb_accum {
	union {
		struct {
			uint64_t count;
			double sum_x, sum_y, sum_z;
			double sxx, sxy, sxz, syy, syz, szz;
		};
		struct {
			uint64_t _pad;
			double cx, cy, cz;
			double bmin_x, bmin_y, bmin_z;
			double bmax_x, bmax_y, bmax_z;
		};
	};
} rintintin_obb_accum;

uint64_t rintintin_oriented_bounding_boxes_scratch_size(uint32_t no_joints)
{
	return (uint64_t)sizeof(rintintin_obb_accum) * no_joints;
}

static rintintin_error_code rintintin_obb_accum_cb(struct rintintin_vertex const* vert, void * user_data)
{
	struct rintintin_bounding_box_command * cmd = (struct rintintin_bounding_box_command*)user_data;
	rintintin_obb_accum * scratch = (rintintin_obb_accum*)cmd->scratch_space;

	int j = rintintin_obb_argmax_joint(vert);
	if (j < 0) return RINTINTIN_SUCCESS;

	rintintin_obb_accum * a = &scratch[j];
	double x = vert->position.x;
	double y = vert->position.y;
	double z = vert->position.z;

	a->count++;
	a->sum_x += x; a->sum_y += y; a->sum_z += z;
	a->sxx += x*x; a->sxy += x*y; a->sxz += x*z;
	a->syy += y*y; a->syz += y*z; a->szz += z*z;
	return RINTINTIN_SUCCESS;
}

static rintintin_error_code rintintin_oriented_bounding_box_cb(struct rintintin_vertex const* vert, void * user_data)
{
	struct rintintin_bounding_box_command * cmd = (struct rintintin_bounding_box_command*)user_data;

	int j = rintintin_obb_argmax_joint(vert);
	if (j < 0) return RINTINTIN_SUCCESS;

	rintintin_inertia_estimation* obb = &cmd->result[j];

	// Centroid source: argmax-cluster centroid (stored in scratch by pass 1) when
	// scratch is provided; otherwise fall back to the weighted mass centroid from
	// metrics. The two differ for bones with skinning bleed; using scratch keeps
	// the rotation and extents anchored on the same point cloud.
	rintintin_vec3 centroid;
	if (cmd->scratch_space != 0L) {
		rintintin_obb_accum const* a = &((rintintin_obb_accum const*)cmd->scratch_space)[j];
		centroid.x = a->sum_x;
		centroid.y = a->sum_y;
		centroid.z = a->sum_z;
	} else {
		centroid = cmd->metrics[j].centroid;
	}

	rintintin_vec3 centered = vec3_sub((dvec3*)&vert->position, &centroid);
	rintintin_vec3 local = rintintin_quat_conjugate_mul_vec3(&obb->rotation, &centered);

	// Update min/max in oriented space (using scale/translation as temporary storage)
	obb->scale.x = (local.x < obb->scale.x) ? local.x : obb->scale.x;
	obb->scale.y = (local.y < obb->scale.y) ? local.y : obb->scale.y;
	obb->scale.z = (local.z < obb->scale.z) ? local.z : obb->scale.z;

	obb->translation.x = (local.x > obb->translation.x) ? local.x : obb->translation.x;
	obb->translation.y = (local.y > obb->translation.y) ? local.y : obb->translation.y;
	obb->translation.z = (local.z > obb->translation.z) ? local.z : obb->translation.z;

	return RINTINTIN_SUCCESS;
}

// Hamilton product (q1 * q2): the resulting rotation applies q2 first, then q1
// when used in q v q* form. We use local-frame perturbation -- best * delta --
// so the perturbation axes track the current box orientation.
static rintintin_vec4 rintintin_obb_quat_mul(rintintin_vec4 a, rintintin_vec4 b)
{
	rintintin_vec4 r;
	r.w = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z;
	r.x = a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y;
	r.y = a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x;
	r.z = a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w;
	return r;
}

// Unit quaternion for a rotation of `theta` radians about the given local axis
// (0=x, 1=y, 2=z). delta_quat(axis, -theta) is the inverse rotation, used to
// revert a rejected trial without storing a backup of the rotation.
static rintintin_vec4 rintintin_obb_delta_quat(int axis, double theta)
{
	double h = theta * 0.5;
	double s = sin(h);
	double c = cos(h);
	rintintin_vec4 q = {.x = 0, .y = 0, .z = 0, .w = c};
	if (axis == 0) q.x = s;
	else if (axis == 1) q.y = s;
	else q.z = s;
	return q;
}

/// Compute oriented bounding boxes from the joint covariance + the joint's
/// argmax-skinned vertex cluster.
///
/// If cmd->scratch_space is non-NULL, the OBB rotation comes from the eigen
/// decomposition of the covariance recomputed from each joint's argmax cluster
/// (rotation and extents share the same point cloud -> self-consistent).
/// If NULL, falls back to using metrics->second_moment for the rotation; the
/// extents pass still uses argmax inclusion either way.
rintintin_error_code rintintin_oriented_bounding_boxes(struct rintintin_bounding_box_command * cmd)
{
	if(cmd == 0L
	|| cmd->meshes == 0L
	|| cmd->metrics == 0L
	|| cmd->result == 0L)
		return RINTINTIN_ERROR_NULL_POINTER;

	if(cmd->result_byte_length != sizeof(*cmd->result) * cmd->no_joints)
		return RINTINTIN_ERROR_INVALID_ARGUMENT;

	rintintin_obb_accum * scratch = (rintintin_obb_accum*)cmd->scratch_space;
	if (scratch != 0L) {
		uint64_t need = rintintin_oriented_bounding_boxes_scratch_size(cmd->no_joints);
		if (cmd->scratch_space_byte_length != need)
			return RINTINTIN_ERROR_INVALID_ARGUMENT;

		// Zero the accumulators (no allocation -- caller-provided buffer).
		for (uint32_t i = 0u; i < cmd->no_joints; ++i)
			scratch[i] = (rintintin_obb_accum){0};

		// Pass 1: accumulate raw moments per argmax joint.
		for (uint64_t mesh_idx = 0; mesh_idx < cmd->no_meshes; ++mesh_idx) {
			int ec = rintintin_visit_each_index(&cmd->meshes[mesh_idx], 0, 1,
												rintintin_obb_accum_cb, cmd);
			if (ec != RINTINTIN_SUCCESS) return ec;
		}
	}

	// Per joint: compute (or look up) the covariance, eigendecompose, set the
	// rotation, and initialise the min/max accumulators for the extents pass.
	for (uint32_t i = 0u; i < cmd->no_joints; ++i)
	{
		rintintin_inertia_estimation* obb = &cmd->result[i];
		rintintin_symmetric_mat3 cov;

		if (scratch != 0L) {
			rintintin_obb_accum * a = &scratch[i];
			if (a->count > 0) {
				double inv = 1.0 / (double)a->count;
				double cx = a->sum_x * inv;
				double cy = a->sum_y * inv;
				double cz = a->sum_z * inv;
				cov.xx = a->sxx * inv - cx*cx;
				cov.yy = a->syy * inv - cy*cy;
				cov.zz = a->szz * inv - cz*cz;
				cov.xy = a->sxy * inv - cx*cy;
				cov.xz = a->sxz * inv - cx*cz;
				cov.yz = a->syz * inv - cy*cz;
				// Replace running sums with the cluster centroid for pass 2.
				a->sum_x = cx;
				a->sum_y = cy;
				a->sum_z = cz;
			} else {
				// No argmax vertices for this joint: covariance is degenerate.
				// Use identity rotation; final extents will be zero.
				cov = (rintintin_symmetric_mat3){0};
				a->sum_x = cmd->metrics[i].centroid.x;
				a->sum_y = cmd->metrics[i].centroid.y;
				a->sum_z = cmd->metrics[i].centroid.z;
			}
		} else {
			cov = cmd->metrics[i].second_moment;
		}

		rintintin_eigen eigen = rintintin_compute_eigen(&cov);
		rintintin_mat3x3 rot_matrix;
		rintintin_compute_rotation_3x3(&rot_matrix, &eigen);

		obb->scale       = (dvec3){ 1e200,  1e200,  1e200};   // min
		obb->translation = (dvec3){-1e200, -1e200, -1e200};   // max
		obb->rotation    = rintintin_quat_from_3x3(&rot_matrix);
	}

	// Pass 2: argmax extents in the PCA-rotated frame.
	for (uint64_t mesh_idx = 0; mesh_idx < cmd->no_meshes; ++mesh_idx) {
		int ec = rintintin_visit_each_index(&cmd->meshes[mesh_idx], 0, 1,
											rintintin_oriented_bounding_box_cb, cmd);
		if (ec != RINTINTIN_SUCCESS) return ec;
	}

	// Optional refinement: PCA gives the best linear-Gaussian fit, which is
	// suboptimal for bent / asymmetric vertex clusters. Coordinate-descent on
	// box volume (one mesh pass per perturbation) finds a local minimum that's
	// usually visibly tighter. Requires scratch (we save min/max backups there
	// to allow rejected trials to revert without an extra recovery pass).
	if (scratch != 0L)
	{
		const double REFINE_STEP_INIT  = 0.20;    // ~11.5 deg starting perturbation
		const double REFINE_STEP_FLOOR = 0.005;   // ~0.29 deg, stop when step shrinks below
		const int    REFINE_MAX_PASSES = 64;      // hard cap on mesh passes

		double step = REFINE_STEP_INIT;
		int passes = 0;
		while (step >= REFINE_STEP_FLOOR && passes < REFINE_MAX_PASSES)
		{
			int improved_this_round = 0;
			for (int axis = 0; axis < 3 && passes < REFINE_MAX_PASSES; ++axis)
			{
				for (int sign = -1; sign <= 1 && passes < REFINE_MAX_PASSES; sign += 2)
				{
					double theta = step * (double)sign;
					rintintin_vec4 d     = rintintin_obb_delta_quat(axis,  theta);
					rintintin_vec4 d_inv = rintintin_obb_delta_quat(axis, -theta);

					// Save current best (min, max) into scratch and apply the
					// perturbation in the local frame: trial_q = best_q * delta.
					// Reset result's min/max accumulators for the trial pass.
					for (uint32_t j = 0; j < cmd->no_joints; ++j) {
						rintintin_inertia_estimation* obb = &cmd->result[j];
						rintintin_obb_accum* a = &scratch[j];
						a->bmin_x = obb->scale.x;
						a->bmin_y = obb->scale.y;
						a->bmin_z = obb->scale.z;
						a->bmax_x = obb->translation.x;
						a->bmax_y = obb->translation.y;
						a->bmax_z = obb->translation.z;
						obb->rotation = rintintin_obb_quat_mul(obb->rotation, d);
						obb->scale       = (dvec3){ 1e200,  1e200,  1e200};
						obb->translation = (dvec3){-1e200, -1e200, -1e200};
					}

					// Trial mesh pass: existing callback fills result[j].scale/translation
					// with trial min/max in the perturbed local frame.
					for (uint64_t mesh_idx = 0; mesh_idx < cmd->no_meshes; ++mesh_idx) {
						int ec = rintintin_visit_each_index(&cmd->meshes[mesh_idx], 0, 1,
															rintintin_oriented_bounding_box_cb, cmd);
						if (ec != RINTINTIN_SUCCESS) return ec;
					}
					++passes;

					// Per joint: compare volumes, accept or revert.
					for (uint32_t j = 0; j < cmd->no_joints; ++j) {
						rintintin_inertia_estimation* obb = &cmd->result[j];
						rintintin_obb_accum* a = &scratch[j];

						double tex = obb->translation.x - obb->scale.x;
						double tey = obb->translation.y - obb->scale.y;
						double tez = obb->translation.z - obb->scale.z;
						double bex = a->bmax_x - a->bmin_x;
						double bey = a->bmax_y - a->bmin_y;
						double bez = a->bmax_z - a->bmin_z;

						// Reject any degenerate/invalid trial (joint with no
						// argmax vertices ends up with min/max still +-inf).
						int trial_valid = (tex > 0 && tey > 0 && tez > 0);
						int best_valid  = (bex > 0 && bey > 0 && bez > 0);

						double trial_vol = trial_valid ? tex * tey * tez : 1e300;
						double best_vol  = best_valid  ? bex * bey * bez : 1e300;

						if (trial_valid && trial_vol < best_vol) {
							improved_this_round = 1;
							// Accept: result is already the trial state.
						} else {
							// Revert rotation and restore (min, max).
							obb->rotation = rintintin_obb_quat_mul(obb->rotation, d_inv);
							obb->scale.x = a->bmin_x; obb->scale.y = a->bmin_y; obb->scale.z = a->bmin_z;
							obb->translation.x = a->bmax_x; obb->translation.y = a->bmax_y; obb->translation.z = a->bmax_z;
						}
					}
				}
			}
			if (!improved_this_round)
				step *= 0.5;
		}
	}

	// Finalize: convert min/max into (half-extents, world-space center).
	for (uint32_t i = 0u; i < cmd->no_joints; ++i)
	{
		rintintin_inertia_estimation* obb = &cmd->result[i];

		rintintin_vec3 centroid;
		if (scratch != 0L) {
			rintintin_obb_accum const* a = &scratch[i];
			centroid.x = a->sum_x;
			centroid.y = a->sum_y;
			centroid.z = a->sum_z;
		} else {
			centroid = cmd->metrics[i].centroid;
		}

		dvec3 local_min = obb->scale;
		dvec3 local_max = obb->translation;

		if (local_min.x < local_max.x)
		{
			rintintin_vec3 center_local = {
				(local_max.x + local_min.x) * 0.5,
				(local_max.y + local_min.y) * 0.5,
				(local_max.z + local_min.z) * 0.5,
			};

			obb->scale.x = (local_max.x - local_min.x) * 0.5;
			obb->scale.y = (local_max.y - local_min.y) * 0.5;
			obb->scale.z = (local_max.z - local_min.z) * 0.5;

			rintintin_vec3 center_world = rintintin_quat_mul_vec3(&obb->rotation, &center_local);
			obb->translation = vec3_add(&centroid, &center_world);
		}
		else
		{
			// No vertices touched this joint.
			obb->scale = (dvec3){0, 0, 0};
			obb->translation = centroid;
		}
	}

	return RINTINTIN_SUCCESS;
}
