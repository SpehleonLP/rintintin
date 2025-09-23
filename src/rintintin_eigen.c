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


static rintintin_eigen rintintin_sort_eigen(rintintin_eigen * in)
{
	int _min = in->values[0] < in->values[1]? 0 : 1;
	int _max = 1 - _min;
	
	_min = in->values[_min] <= in->values[2]? _min : 2;
	_max = in->values[_max] > in->values[2]? _max : 2;
	//the three indices must sum to 0+1+2=3.
	int _mid = 3 - _min - _max;
	
	if(_mid == -1)
	{
		int break_point = 0;
		--break_point;
	
	}
	
	rintintin_eigen r = {
		.values={in->values[_min], in->values[_mid], in->values[_max]},
		.vectors={in->vectors[_min], in->vectors[_mid], in->vectors[_max]},
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
  
#define getElement(a, b, c) (a[b][c])

rintintin_eigen rintintin_compute_eigen(rintintin_symmetric_mat3 const*I) {
    dmat3 A = {.m={
		.x={I->xx, I->xy, I->xz},
		.y={I->xy, I->yy, I->yz},
		.z={I->xz, I->yz, I->zz}    
    }};
    
    return rintintin_compute_eigen_m3(&A);
}


rintintin_eigen rintintin_compute_eigen_m3(rintintin_mat3x3 *I) {
    double (*A)[3] = I->dbl;
    
    double V[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    
    const int maxIterations = 50;
    
    // Compute adaptive tolerance based on matrix scale
    double matrixScale = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            matrixScale += fabs(getElement(A, i, j));
        }
    }
    
    if(matrixScale == 0.0)
    {
		return (rintintin_eigen)
		{
			.vectors = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
			.values = {0, 0, 0}
		};
    }
    
    matrixScale /= 9.0; // Average absolute value of matrix elements
    
	// Scale matrix up to reasonable range
	double scale_factor = 1.0;
	if (matrixScale < 1e-6) {  // or whatever threshold makes sense
		scale_factor = 1e6 / matrixScale;  // bring into ~1.0 range
	}

	for(int i = 0; i < 9; ++i)		
		A[i/3][i%3] *= scale_factor;

    // Adaptive tolerance: scale with matrix magnitude, but keep reasonable bounds
    double tolerance = fmax(1e-12, fmin(1e-6, matrixScale * scale_factor * 1e-6));
  
	int iter = 0;
    for (iter = 0; iter < maxIterations; iter++) {
        int p = 0, q = 1;
        double maxOffDiag = 0.0;
        
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 3; j++) {
                double val = fabs(getElement(A, i, j));
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
        
        if (fabs(getElement(A, p, q)) < 1e-15) continue;
        
        double tau = (getElement(A, q, q) - getElement(A, p, p)) / (2.0 * getElement(A, p, q));
        double t = (tau >= 0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
        double c = 1.0 / sqrt(1.0 + t * t);
        double s = t * c;
        
        double App = getElement(A, p, p);
        double Aqq = getElement(A, q, q);
        double Apq = getElement(A, p, q);
        
        getElement(A, p, p) = c * c * App - 2.0 * s * c * Apq + s * s * Aqq;
        getElement(A, q, q) = s * s * App + 2.0 * s * c * Apq + c * c * Aqq;
        getElement(A, p, q) = getElement(A, q, p) = 0.0;
        
        for (int r = 0; r < 3; r++) {
            if (r != p && r != q) {
                double Arp = getElement(A, r, p);
                double Arq = getElement(A, r, q);
                getElement(A, r, p) = getElement(A, p, r) = c * Arp - s * Arq;
                getElement(A, r, q) = getElement(A, q, r) = s * Arp + c * Arq;
            }
        }
        
        for (int i = 0; i < 3; i++) {
            double Vip = getElement(V, i, p);
            double Viq = getElement(V, i, q);
            getElement(V, i, p) = c * Vip - s * Viq;
            getElement(V, i, q) = s * Vip + c * Viq;
        }
    }
    
    rintintin_eigen result;
    for (int i = 0; i < 3; i++) {
        result.values[i] = getElement(A, i, i) / scale_factor;
        result.vectors[i].x = getElement(V, 0, i);
        result.vectors[i].y = getElement(V, 1, i);
        result.vectors[i].z = getElement(V, 2, i);
    }
        
    makeConsistent(&result);    
    return rintintin_sort_eigen(&result);
}

static rintintin_vec4 quat_multiply(rintintin_vec4 const* a, rintintin_vec4 const* b)
{
   rintintin_vec4 result;
   
   // Quaternion multiplication: (w1 + x1*i + y1*j + z1*k) * (w2 + x2*i + y2*j + z2*k)
   result.w = a->w * b->w - a->x * b->x - a->y * b->y - a->z * b->z;
   result.x = a->w * b->x + a->x * b->w + a->y * b->z - a->z * b->y;
   result.y = a->w * b->y - a->x * b->z + a->y * b->w + a->z * b->x;
   result.z = a->w * b->z + a->x * b->y - a->y * b->x + a->z * b->w;
   
   return result;
}

rintintin_error_code rintintin_estimate_shapes(rintintin_inertia_estimation * dst, rintintin_metrics const* src, uint64_t no_items)
{
	if(dst == 0 || src == 0)
		return RINTINTIN_ERROR_NULL_POINTER;

	// 90° rotations as quaternions (axis, angle)
	const rintintin_vec4 rot_around_x = {.x=0.7071, .y=0, .z=0, .w=0.7071}; // 90° around X
	const rintintin_vec4 rot_around_y = {.x=0, .y=0.7071, .z=0, .w=0.7071}; // 90° around Y
	const rintintin_vec4 rot_around_z = {.x=0, .y=0.0, .z=0.7071, .w=0.7071}; // 90° around Z

    const double pi = 3.14159265358979323846;
    const double sphere_factor = 0.62035049089; // ratio for converting unit RADIUS sphere to unit VOLUME sphere
    (void)sphere_factor;
    (void)rot_around_x;
    (void)rot_around_y;
    (void)rot_around_z;
    
    for (uint64_t i = 0; i < no_items; i++) {
        // Compute eigendecomposition of the inertia tensor
        rintintin_symmetric_mat3 inertia = src[i].inertia;
        rintintin_eigen eigen = rintintin_compute_eigen(&inertia);
        
        // Set translation to the centroid
        dst[i].translation = src[i].centroid;
        
        // Set rotation based on principal axes
        dst[i].rotation = rintintin_compute_rotation_quat(&eigen);
        
        // always positive for a valid inertia tensor anyway... 
        double lambda_small = fabs(eigen.values[0]); // smallest eigenvalue
        double lambda_mid = fabs(eigen.values[1]);   // middle eigenvalue  
        double lambda_large = fabs(eigen.values[2]); // largest eigenvalue
        double mass = src[i].volume; // treating volume as mass (density = 1)
        
        double best_error = mass*1e6;
        (void)best_error;
        
        rintintin_vec3 best_scale = {0.01, 0.01, 0.01};
        int best_axis = 2;
        (void)best_axis;
        
        if(mass == 0)
        {
            dst[i].scale = (rintintin_vec3){0.0001, 0.0001, 0.0001};
            continue;            
        }
        
        // Avoid division by zero
        if (mass < 1e-30 || lambda_large < 1e-30) {
            dst[i].scale = (rintintin_vec3){0.01, 0.01, 0.01};
            continue;
        }
   
        // Test ELLIPSOID
        {
			double large_axis = 5.0 * (lambda_mid + lambda_large - lambda_small) / (2.0 * mass);
			double mid_axis = 5.0 * (lambda_small + lambda_large - lambda_mid) / (2.0 * mass);
			double small_axis = 5.0 * (lambda_small + lambda_mid - lambda_large) / (2.0 * mass);
			
			int update_error = (large_axis > 0 && mid_axis > 0 && small_axis > 0);
			
            large_axis = sqrt(fabs(large_axis)) * sphere_factor;
            mid_axis = sqrt(fabs(mid_axis)) * sphere_factor;
            small_axis = sqrt(fabs(small_axis)) * sphere_factor;
            
			double predicted_volume = (4.0/3.0) * pi * large_axis * mid_axis * small_axis;
			double error = fabs(predicted_volume - mass) / mass;
		
			if(update_error)
				best_error = error;
				
			best_scale = (rintintin_vec3){large_axis, mid_axis, small_axis}; // x=minor, z=major   
        }
#if 0
        // Test CUBE
		// Test BOX (rectangular cuboid)
		{
			// For a box with dimensions a, b, c:
			// I_xx = (1/12)*m*(b² + c²)
			// I_yy = (1/12)*m*(a² + c²)  
			// I_zz = (1/12)*m*(a² + b²)
			// 
			// Solving for dimensions:
			// a² = (6/m) * (I_yy + I_zz - I_xx)
			// b² = (6/m) * (I_xx + I_zz - I_yy)
			// c² = (6/m) * (I_xx + I_yy - I_zz)
			
			double a_squared = (6.0 / mass) * (lambda_mid + lambda_large - lambda_small);
			double b_squared = (6.0 / mass) * (lambda_small + lambda_large - lambda_mid);
			double c_squared = (6.0 / mass) * (lambda_small + lambda_mid - lambda_large);
			
			if (a_squared > 0 && b_squared > 0 && c_squared > 0) {
				double a = sqrt(a_squared);
				double b = sqrt(b_squared);
				double c = sqrt(c_squared);
				
				double predicted_volume = a * b * c;
			//	if(predicted_volume < mass)
				{
					double error = fabs(predicted_volume - mass) / mass;
					
					if (error < best_error) {
						best_error = error;
						best_type = RINTINTIN_UNIT_CUBE; // You'll need to define this constant
						best_scale = (rintintin_vec3){a, b, c};
					}
				}
			}
		}
#endif

#if 0
		// Test CYLINDER (try each axis as the cylinder axis)
		{
			for (int axis = 0; axis < 3; axis++) {
				double axis_eigenvalue, perp_eigenvalue1, perp_eigenvalue2;
				
				if (axis == 0) { // X-axis is cylinder axis
					axis_eigenvalue = lambda_small;
					perp_eigenvalue1 = lambda_mid;
					perp_eigenvalue2 = lambda_large;
				} else if (axis == 1) { // Y-axis is cylinder axis  
					axis_eigenvalue = lambda_mid;
					perp_eigenvalue1 = lambda_small;
					perp_eigenvalue2 = lambda_large;
				} else { // Z-axis is cylinder axis
					axis_eigenvalue = lambda_large;
					perp_eigenvalue1 = lambda_small;
					perp_eigenvalue2 = lambda_mid;
				}
				
				// For cylinder: I_perp = (1/12)*m*(3*r² + h²), I_axis = (1/2)*m*r²
				double r_squared = 2.0 * axis_eigenvalue / mass;
				double avg_perp = (perp_eigenvalue1 + perp_eigenvalue2) / 2.0;
				double h_squared = 12.0 * (avg_perp - axis_eigenvalue/2.0) / mass;
				
				if (r_squared > 0 && h_squared > 0) {
					double r = sqrt(r_squared);
					double h = sqrt(h_squared);
					double predicted_volume = pi * r * r * h;
					double error = fabs(predicted_volume - mass) / mass;
					
					if (predicted_volume < mass && error < best_error) {
						best_error = error;
						best_type = RINTINTIN_UNIT_CYLINDER;
						best_axis = axis;
						
						best_scale = (rintintin_vec3){r, r, h};
					}
				}
			}
		}
#endif
#if 0
        // Test CONE (try each axis as the cone axis)
		{
			for (int axis = 0; axis < 3; axis++) {
				double axis_eigenvalue, perp_eigenvalue1, perp_eigenvalue2;
				
				if (axis == 0) { // X-axis is cone axis
					axis_eigenvalue = lambda_small;
					perp_eigenvalue1 = lambda_mid;
					perp_eigenvalue2 = lambda_large;
				} else if (axis == 1) { // Y-axis is cone axis
					axis_eigenvalue = lambda_mid;
					perp_eigenvalue1 = lambda_small;
					perp_eigenvalue2 = lambda_large;
				} else { // Z-axis is cone axis
					axis_eigenvalue = lambda_large;
					perp_eigenvalue1 = lambda_small;
					perp_eigenvalue2 = lambda_mid;
				}
				
				// For cone: I_perp = (3/80)*m*(4*r² + h²), I_axis = (3/10)*m*r²
				double r_squared = 10.0 * axis_eigenvalue / (3.0 * mass);
				double avg_perp = (perp_eigenvalue1 + perp_eigenvalue2) / 2.0;
				double h_squared = 80.0 * (avg_perp - 3.0*axis_eigenvalue/20.0) / (3.0 * mass);
				
				if (r_squared > 0 && h_squared > 0) {
					double r = sqrt(r_squared);
					double h = sqrt(h_squared);
					double predicted_volume = (1.0/3.0) * pi * r * r * h;
					double error = fabs(predicted_volume - mass) / mass;
					
					if (error < best_error) {
						best_error = error;
						best_type = RINTINTIN_UNIT_CONE;
						best_scale = (rintintin_vec3){r, h, r};
						best_axis = axis;
					}
				}
			}
		}
		
#endif		
		if (best_axis == 0) { // Z axis needs to be rotated to face X direction
			dst[i].rotation = quat_multiply(&dst[i].rotation, &rot_around_y);
		} 
		else if (best_axis == 1) { // Z-axis needs to align with Y  
			dst[i].rotation = quat_multiply(&dst[i].rotation, &rot_around_x);
		}


        dst[i].scale = best_scale;
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
