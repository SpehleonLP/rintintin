#ifndef RINTINTIN_CONSTRAINT_SOLVER_H
#define RINTINTIN_CONSTRAINT_SOLVER_H
#include "../include/rintintin.h"
#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_latent_methods.h"
#include "rintintin_latentspace.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct rintintin_constrained_config rintintin_constrained_config;
typedef struct rintintin_constrained_result rintintin_constrained_result;

struct rintintin_constrained_config
{
	const char * name;
	int idx;
	int parent;
	
	int max_iterations;
	bool verbose;
	bool has_children;
	
	rintintin_vec3 const* centroid;
	rintintin_vec3 const* joint;	
	
	float alpha, lambda;
};

struct rintintin_constrained_result
{
	rintintin_vec3 solution;
	double residual;
	bool converged;
	bool improved;
	int8_t error_code;
	int8_t calc;
};

struct rintintin_tensor;

static inline rintintin_constrained_config rintintin_constrained_default_config(rintintin_vec3 const* centroid, rintintin_vec3 const* joint)
{
	rintintin_constrained_config r = {0};
	r.name = 0L;
	r.idx = 0;
	r.parent = 0;
	r.verbose = false;
	r.has_children=false;
	r.centroid = centroid;
	r.joint = joint;
	r.max_iterations = 10;
	r.lambda = 1e-6f;
	r.alpha = 1e-2f;
	return r;
}
 
rintintin_constrained_result rintintin_solve_center(rintintin_liquified const* tensor, rintintin_vec3 const* guess, rintintin_constrained_config const* config);
rintintin_constrained_result rintintin_walk_to_definite_matrix(rintintin_tensor const* tensor, rintintin_vec3 const* centroid, rintintin_constrained_config const* config);
		
#ifdef __cplusplus
}
#endif

#if USING_NAMESPACE_RINTINTIN
typedef rintintin_constrained_result constrained_result;
typedef rintintin_constrained_config config_t;
#endif

#endif // RINTINTIN_CONSTRAINT_SOLVER_H
