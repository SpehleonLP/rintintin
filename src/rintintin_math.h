#include "../include/rintintin.h"
#include "rintintin_sqrt.h"

#if USING_NAMESPACE_RINTINTIN == 1
typedef rintintin_vec3 dvec3;
typedef rintintin_vec4 dvec4;
#define vec3_add rintintin_vec3_add
#define vec3_sub rintintin_vec3_sub
#define vec3_mul rintintin_vec3_mul
#define dot rintintin_vec3_dot
#define cross rintintin_vec3_cross

#define vec3_scale rintintin_vec3_scale
#define vec3_len rintintin_vec3_scale
#define vec3_normalize rintintin_vec3_normalize

#define vec4_normalize rintintin_vec3_normalize
#endif

static inline rintintin_vec3 rintintin_vec3_add(rintintin_vec3 const* a, rintintin_vec3 const* b)
{
	return (rintintin_vec3){a->x + b->x, a->y + b->y, a->z + b->z};
}

static inline rintintin_vec3 rintintin_vec3_sub(rintintin_vec3 const* a, rintintin_vec3 const* b)
{
	return (rintintin_vec3){a->x - b->x, a->y - b->y, a->z - b->z};
}

static inline rintintin_vec3 rintintin_vec3_mul(rintintin_vec3 const* a, rintintin_vec3 const* b)
{
	return (rintintin_vec3){
		.x = a->x*b->x,
		.y = a->y*b->y,
		.z = a->z*b->z,
	};
}

static inline double rintintin_vec3_dot(rintintin_vec3 const* a, rintintin_vec3 const* b)
{
	return a->x*b->x + a->y*b->y + a->z*b->z;
}

static inline rintintin_vec3 rintintin_vec3_cross(rintintin_vec3 const* a, rintintin_vec3 const* b)
{
	rintintin_vec3 r = {
		a->y*b->z - a->z*b->y,
		a->z*b->x - a->x*b->z,
		a->x*b->y - a->y*b->x
	};

	return r;
}

static inline rintintin_vec3 rintintin_vec3_scale(rintintin_vec3 const* a, double b)
{
	return (rintintin_vec3){
		.x = a->x*b,
		.y = a->y*b,
		.z = a->z*b,
	};
}

static inline double rintintin_vec3_len(rintintin_vec3 const* in)
{
	double length2 = (in->x*in->x + in->y*in->y + in->z*in->z);
	return sqrt(length2);
}

static inline rintintin_vec3 rintintin_vec3_normalize(rintintin_vec3 const* in)
{
	double length2 = (in->x*in->x + in->y*in->y + in->z*in->z);
	length2 = length2 > 0? 1.0 / sqrt(length2) : 0.0;
	rintintin_vec3 r = {in->x * length2, in->y * length2, in->z * length2};
	return r;
}

static inline rintintin_vec4 rintintin_vec4_normalize(rintintin_vec4 const* in)
{
	double length2 = (in->x*in->x + in->y*in->y + in->z*in->z + in->w*in->w);
	length2 = length2 > 0? 1.0 / sqrt(length2) : 0.0;
	rintintin_vec4 r;
	r.x = in->x * length2;
	r.y = in->y * length2;
	r.z = in->z * length2;
	r.w = in->w * length2;;
	return r;
}

static inline rintintin_vec3 vmin(rintintin_vec3 const* a, rintintin_vec3 const* b) 
{ return (rintintin_vec3){ (a->x < b->x? a->x : b->x), (a->y < b->y? a->y : b->y), (a->z < b->z? a->z : b->z) }; }
static inline rintintin_vec3 vmax(rintintin_vec3 const* a, rintintin_vec3 const* b) 
{ return (rintintin_vec3){ (a->x > b->x? a->x : b->x), (a->y > b->y? a->y : b->y), (a->z > b->z? a->z : b->z) }; }


///---- SMAT3
 
#if USING_NAMESPACE_RINTINTIN == 1
typedef rintintin_symmetric_mat3 smat3;

#define smat3_add rintintin_smat3_add
#define smat3_sub rintintin_smat3_sub
#define parallel_axis rintintin_compute_parallel_axis
#define smat3_square rintintin_compute_square_sym
#endif

static inline rintintin_symmetric_mat3 rintintin_smat3_add(rintintin_symmetric_mat3 const* a, rintintin_symmetric_mat3 const* b)
{
	return  (rintintin_symmetric_mat3){
		.xx = a->xx + b->xx,
		.yy = a->yy + b->yy,
		.zz = a->zz + b->zz,
		.xy = a->xy + b->xy,
		.xz = a->xz + b->xz,
		.yz = a->yz + b->yz,
	};
}

static inline rintintin_symmetric_mat3 rintintin_smat3_sub(rintintin_symmetric_mat3 const* a, rintintin_symmetric_mat3 const* b)
{
	return  (rintintin_symmetric_mat3){
		.xx = a->xx - b->xx,
		.yy = a->yy - b->yy,
		.zz = a->zz - b->zz,
		.xy = a->xy - b->xy,
		.xz = a->xz - b->xz,
		.yz = a->yz - b->yz,
	};
}

static inline rintintin_symmetric_mat3 rintintin_compute_parallel_axis(double mass, rintintin_vec3 const* o)
{ 
	return  (rintintin_symmetric_mat3){
		.xx = mass*o->x*o->x,
		.yy = mass*o->y*o->y,
		.zz = mass*o->z*o->z,
		.xy = -mass*o->x*o->y,
		.xz = -mass*o->x*o->z,
		.yz = -mass*o->y*o->z,
	};
}

static inline rintintin_symmetric_mat3 rintintin_compute_square_sym(rintintin_symmetric_mat3 const* in, rintintin_vec3 const * p)
{
	return (rintintin_symmetric_mat3)
	{
		.xx = in->xx * p->x*p->x,
		.yy = in->yy * p->y*p->y,
		.zz = in->zz * p->z*p->z,
		.xy = in->xy * p->x*p->y,
		.xz = in->xz * p->x*p->z,
		.yz = in->yz * p->y*p->z,
	};
}


/// ---- MAT3

#if USING_NAMESPACE_RINTINTIN == 1
typedef rintintin_mat3x3 dmat3;

#define mat3_mul_vec3 rintintin_mat3x3_mul_vec
#define mat3T_mul_vec3 rintintin_mat3x3_T_mul_vec

#define mat3_add rintintin_mat3x3_add
#define mat3_sub rintintin_mat3x3_sub
#define mat3_mul rintintin_mat3x3_mul
#define mat3_mul_T rintintin_mat3x3_mul_T

#define mat3_scale_eq rintintin_mat3x3_scale_eq
#define inverse rintintin_compute_inverse_3x3
#endif


// Column-oriented multiplication: treats vector as column vector, matrix * [x; y; z]
// Result = matrix * [x; y; z] (standard matrix-vector multiplication)
static inline rintintin_vec3 rintintin_mat3x3_mul_vec(const rintintin_mat3x3* mat, const rintintin_vec3* vec) {
    rintintin_vec3 result;
    
    result.x = rintintin_vec3_dot(&mat->m.x, vec);
    result.y = rintintin_vec3_dot(&mat->m.y, vec);
    result.z = rintintin_vec3_dot(&mat->m.z, vec);
    
    return result;
}

// Row-oriented multiplication: treats vector as row vector [x y z] * matrix
// Result = [x y z] * [dx dy dz] = x*dx + y*dy + z*dz (linear combination of matrix columns)
static inline rintintin_vec3 rintintin_mat3x3_T_mul_vec(const rintintin_mat3x3* mat, const rintintin_vec3* vec) {
    rintintin_vec3 result;
    
    result.x = vec->x * mat->dbl[0][0] + vec->y * mat->dbl[1][0] + vec->z * mat->dbl[2][0];
    result.y = vec->x * mat->dbl[0][1] + vec->y * mat->dbl[1][1] + vec->z * mat->dbl[2][1];
    result.z = vec->x * mat->dbl[0][2] + vec->y * mat->dbl[1][2] + vec->z * mat->dbl[2][2];
    
    return result;
}

static inline rintintin_mat3x3 rintintin_mat3x3_add(const rintintin_mat3x3* a, const rintintin_mat3x3* b)
{
	rintintin_mat3x3 c;
	
	for(int i = 0; i < 9; ++i)
		(&(c.dbl[0][0]))[i] = (&(a->dbl[0][0]))[i] + (&(b->dbl[0][0]))[i];

	return c;
}

static inline rintintin_mat3x3 rintintin_mat3x3_mul(const rintintin_mat3x3* a, const rintintin_mat3x3* b)
{
   rintintin_mat3x3 result;
   
   // Standard matrix multiplication: C[i][j] = sum over k of A[i][k] * B[k][j]
   for (int i = 0; i < 3; i++) {
       for (int j = 0; j < 3; j++) {
           result.dbl[i][j] = a->dbl[i][0] * b->dbl[0][j] +
                              a->dbl[i][1] * b->dbl[1][j] +
                              a->dbl[i][2] * b->dbl[2][j];
       }
   }
   
   return result;   
}


static inline rintintin_mat3x3 rintintin_mat3x3_mul_T(const rintintin_mat3x3* mat)
{
	rintintin_mat3x3 r;
		
	// First term: 2 * J^T * J (Gauss-Newton approximation)
	r.m.x.x = (mat->m.x.x * mat->m.x.x + mat->m.y.x * mat->m.y.x + mat->m.z.x * mat->m.z.x);
	r.m.x.y = (mat->m.x.x * mat->m.x.y + mat->m.y.x * mat->m.y.y + mat->m.z.x * mat->m.z.y);
	r.m.x.z = (mat->m.x.x * mat->m.x.z + mat->m.y.x * mat->m.y.z + mat->m.z.x * mat->m.z.z);
	r.m.y.y = (mat->m.x.y * mat->m.x.y + mat->m.y.y * mat->m.y.y + mat->m.z.y * mat->m.z.y);
	r.m.y.z = (mat->m.x.y * mat->m.x.z + mat->m.y.y * mat->m.y.z + mat->m.z.y * mat->m.z.z);
	r.m.z.z = (mat->m.x.z * mat->m.x.z + mat->m.y.z * mat->m.y.z + mat->m.z.z * mat->m.z.z);
	r.m.y.x = r.m.x.y;
	r.m.z.x = r.m.x.z;
	r.m.z.y = r.m.y.z;

	return r;
}

static inline  rintintin_mat3x3 rintintin_mat3x3_transpose(const rintintin_mat3x3* mat) {
	rintintin_mat3x3 r;
	
	r.dbl[0][0] = mat->dbl[0][0];
	r.dbl[0][1] = mat->dbl[1][0];
	r.dbl[0][2] = mat->dbl[2][0];
	r.dbl[1][0] = mat->dbl[0][1];
	r.dbl[1][1] = mat->dbl[1][1];
	r.dbl[1][2] = mat->dbl[2][1];
	r.dbl[2][0] = mat->dbl[0][2];
	r.dbl[2][1] = mat->dbl[1][2];
	r.dbl[2][2] = mat->dbl[2][2];
	
    return r;
}



/*
rintintin_mat3x3 rintintin_hessain_mul_vec3(rintintin_hessian *const a, rintintin_vec3*const b)
{
	rintintin_mat3x3 H;

	// Second term: 2 * Σ(F_i * H_i) - add constraint-weighted second derivatives
	// a->m.x = Hessian of constraint fx, a->m.y = Hessian of constraint fy, a->m.z = Hessian of constraint fz
	H.m.x.x = (b->x * a->m.x.m.x.x + b->y * a->m.y.m.x.x + b->z * a->m.z.m.x.x);
	H.m.x.y = (b->x * a->m.x.m.x.y + b->y * a->m.y.m.x.y + b->z * a->m.z.m.x.y);
	H.m.x.z = (b->x * a->m.x.m.x.z + b->y * a->m.y.m.x.z + b->z * a->m.z.m.x.z);
	H.m.y.x = (b->x * a->m.x.m.y.x + b->y * a->m.y.m.y.x + b->z * a->m.z.m.y.x);
	H.m.y.y = (b->x * a->m.x.m.y.y + b->y * a->m.y.m.y.y + b->z * a->m.z.m.y.y);
	H.m.y.z = (b->x * a->m.x.m.y.z + b->y * a->m.y.m.y.z + b->z * a->m.z.m.y.z);
	H.m.z.x = (b->x * a->m.x.m.z.x + b->y * a->m.y.m.z.x + b->z * a->m.z.m.z.x);
	H.m.z.y = (b->x * a->m.x.m.z.y + b->y * a->m.y.m.z.y + b->z * a->m.z.m.z.y);
	H.m.z.z = (b->x * a->m.x.m.z.z + b->y * a->m.y.m.z.z + b->z * a->m.z.m.z.z);

	return H;
}

rintintin_mat3x3 rintintin_mat3x3_add(const rintintin_mat3x3* a, const rintintin_mat3x3* b) 
{
	rintintin_mat3x3 r;
	
	for(int i = 0; i < 9; ++i)
	{
		(&r.m.x.x)[i] = (&a->m.x.x)[i] + (&b->m.x.x)[i];
	}
	
	return r;
}



typedef union rintintin_hessian {
	struct {
		rintintin_mat3x3 x;
		rintintin_mat3x3 y; 
		rintintin_mat3x3 z;
	} m;
	double dbl[3][3][3];
	rintintin_mat3x3 layers[3];
	rintintin_vec3 h[3][3];
} rintintin_hessian;

void rintintin_mat3x3_scale_eq(rintintin_mat3x3 * a, double b);
rintintin_mat3x3 rintintin_hessain_mul_vec3(rintintin_hessian *const a, rintintin_vec3*const b);
*/
