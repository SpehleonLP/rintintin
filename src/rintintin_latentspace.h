#ifndef RINTINTIN_LATENTSPACE_H
#define RINTINTIN_LATENTSPACE_H
#include "../include/rintintin.h"

enum
{
	THIN_MASS = 0,
	THICK_MASS = 3,
};


struct rintintin_input
{
	rintintin_vec3 p0, p1, p2;
	double d0, d1, d2;
};

typedef struct rintintin_tensor rintintin_tensor;
typedef struct rintintin_latent_space rintintin_latent_space;
typedef struct rintintin_covariance rintintin_covariance;

/**
 * @brief Axis-aligned bounding box
 */
typedef struct rintintin_aabb {
    rintintin_vec3 min;
    rintintin_vec3 max;
} rintintin_aabb;

struct rintintin_latent_space
{
	double thick[58];
	double thin[10];
};

// Predefined structs

typedef struct rintintin_o1 {
   double x_xx;
   double x_xy;
   double x_xz;
   double y_xy;
   double y_yy;
   double y_yz;
   double z_xz;
   double z_yz;
   double z_zz;
} rintintin_o1;

typedef struct rintintin_sym_jacobian
{
    rintintin_vec3 xx;
    rintintin_vec3 yy;
    rintintin_vec3 zz;
    rintintin_vec3 xy;
    rintintin_vec3 xz;
    rintintin_vec3 yz;
} rintintin_sym_jacobian;


typedef struct rintintin_coupling {
    rintintin_vec3 x;
    rintintin_vec3 y;
    rintintin_vec3 z;
} rintintin_coupling;

typedef struct rintintin_linear {
    rintintin_coupling coupling;
    rintintin_vec3 o;
    rintintin_sym_jacobian j;
} rintintin_linear;


typedef struct rintintin_quadratic {
#define RINTINTIN_QUAD_STRUCT(C, ...) struct rintintin_##C { double __VA_ARGS__; } C;

	RINTINTIN_QUAD_STRUCT(xx, xx, xy, xz)
	RINTINTIN_QUAD_STRUCT(yy, yy, xy, yz)
	RINTINTIN_QUAD_STRUCT(zz, zz, xz, yz)
	
	RINTINTIN_QUAD_STRUCT(xy, xx, yy, xy, xz, yz)
	RINTINTIN_QUAD_STRUCT(xz, xx, zz, xy, xz, yz)
	RINTINTIN_QUAD_STRUCT(yz, yy, zz, xy, yz, xz)
#undef RINTINTIN_QUAD_STRUCT
} rintintin_quadratic;

typedef struct rintintin_cubic {
    rintintin_o1 o1;
    rintintin_vec3 o2;
    rintintin_sym_jacobian h;
} rintintin_cubic;

struct rintintin_tensor {
    rintintin_symmetric_mat3 constant;
    rintintin_linear linear;
    rintintin_quadratic quadratic;
    rintintin_cubic cubic;
    double mass_o;
};

typedef struct rintintin_solidified
{
    rintintin_symmetric_mat3 constant;
    rintintin_vec3 linear_o;
    double mass_o;
} rintintin_solidified;

typedef struct rintintin_liquified
{
	struct cubic_liquid
	{
		rintintin_o1 o1;
		rintintin_vec3 o2;	
	} cubic;

	struct linear_liquid
	{
		rintintin_coupling coupling;
		rintintin_vec3 o;
	} linear;
	
    double mass_o;
} rintintin_liquified;

void rintintin_thin_construct(double * coeff, struct rintintin_input const* src, double thickness);
void rintintin_thick_construct(double * coeff, struct rintintin_input const* src);
rintintin_vec3 rintintin_v_func_from_coeff(double const* coeff);
rintintin_tensor rintintin_tensor_from_coeff(double const* coeff);
     
rintintin_liquified rintintin_liquid_from_coeff(double const* coeff);

rintintin_vec3 rintintin_evaluate_center(rintintin_liquified const * it, rintintin_vec3 const* p);
int rintintin_evaluate_constraints(rintintin_tensor const * it, rintintin_vec3 const* p, rintintin_vec3 * out);

#if USING_NAMESPACE_RINTINTIN
typedef rintintin_latent_space latent_t;
typedef rintintin_tensor tensor_t;
typedef rintintin_liquified liquid_t;
typedef rintintin_solidified solid_t;


#endif


#endif // RINTINTIN_LATENTSPACE_H
