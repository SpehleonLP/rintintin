
#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_eigen.h"
#include "rintintin_sqrt.h"
#include "rintintin_math.h"
#include "rintintin_constraint_solver.h"
#include "rintintin_latent_methods.h"

static inline dvec3 mix(dvec3 const* a, dvec3 const* b, double alpha)
{
	return (dvec3){
		a->x * (1.0 - alpha) + b->x * alpha, 
		a->y * (1.0 - alpha) + b->y * alpha, 
		a->z * (1.0 - alpha) + b->z * alpha, 	
	};
}

#define abs(x) ((x) < 0? -(x) : (x))

double frobenius_norm(double (* m)[3]) {
    double sum = 0;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            sum += m[i][j] * m[i][j];
        }
    }
    
    return rintintin_sqrt(sum);
}

#define cross rintintin_vec3_cross

dmat3 get_matrix(dvec3 const* in)
{
	double len2 = dot(in, in);

	if(len2 == 0)
		return (dmat3){.m={
			.x={1, 0, 0},
			.y={0, 1, 0},
			.z={0, 0, 1}
		}};
			
	
	dvec3 norm = vec3_scale(in, 1.0/sqrt(len2));
	
	if(abs(norm.x) < 0.747)
	{
		dvec3 norm_y = cross(&((dvec3){1.0, 0.0, 0.0}), &norm);
		len2 = dot(&norm_y, &norm_y);
		norm_y = vec3_scale(&norm_y, 1.0/sqrt(len2));
			
		return (dmat3){.m={
			.x=cross(&norm_y, &norm),
			.y=norm_y,
			.z=norm
		}};
	}	
	else
	{
		dvec3 norm_x = cross(&((dvec3){0.0, 1.0, 0.0}), &norm);
		len2 = dot(&norm_x, &norm_x);
		norm_x = vec3_scale(&norm_x, 1.0/sqrt(len2));
			
		return (dmat3){.m={
			.x=norm_x,
			.y=cross(&norm_x, &norm),
			.z=norm
		}};
	}	
}

dmat3 get_matrix_in_plane(dvec3 const* src, dvec3 const* dst)
{
	dvec3 delta = vec3_sub(dst, src);
	double scale = sqrt(dot(&delta, &delta));
	dmat3 m1 = get_matrix(&delta);
	return (dmat3){.m={
		.x=vec3_scale(&m1.m.z, scale),
		.y=vec3_scale(&m1.m.y, scale),
		.z=vec3_scale(&m1.m.x, scale),
	}};
}


constrained_result rintintin_walk_to_definite_matrix(tensor_t const* tensor, dvec3 const* centroid, config_t const* config);


dvec3 get_candidate_solution(rintintin_liquified const* it);


static inline bool check_classification(tensor_t const* tensor, dvec3 const* centroid, int clac)
{
	enum 
	{
		SIGN_MASK = RINTINTIN_NEGATIVE|RINTINTIN_POSITIVE,
		DEFINITIE_MASK = RINTINTIN_DEFINITE|RINTINTIN_SEMI_DEFINITE
	};
	
	if((clac & DEFINITIE_MASK) == 0)
		return false;
		
	double mass = tensor->mass_o + dot(&tensor->cubic.o2, centroid);
	return (clac & SIGN_MASK) == (mass < 0? RINTINTIN_NEGATIVE : RINTINTIN_POSITIVE);
}

constrained_result rintintin_walk_to_definite_matrix(tensor_t const* tensor, dvec3 const* centroid, config_t const* config)
{	
	dvec3 constraints;
	dvec3 o2 = tensor->cubic.o2;
	double o2_len2 = o2.x*o2.x + o2.y*o2.y + o2.z*o2.z;
	double mass = tensor->mass_o + dot(&tensor->cubic.o2, centroid);
	
	int clac = rintintin_evaluate_constraints(tensor, centroid, &constraints); 
	
	if(o2_len2 == 0 || check_classification(tensor, centroid, clac))
	{
		return (rintintin_constrained_result){
			.solution = *centroid,
			.residual = constraints.x + constraints.y + constraints.z,
			.converged = (clac != 0),
			.improved = false,
			.error_code = o2_len2? 0 : -2,
			.calc=(char)clac 
		};		
	}
	
	// Extract o2 coefficients (mass gradient)
	dvec3 dir = tensor->cubic.o2;
	dir = vec3_scale(&dir, (mass < 0? -1.0 : 1.0) / sqrt(o2_len2));
	
	double max = 1.0;
	double min = 0;
		
	int iterations = 0;
	dvec3 p = *centroid;
	for(; iterations < config->max_iterations; ++iterations, min = max, max *= 2)
	{
		p = vec3_scale(&dir, max);
		p = vec3_add(&p, centroid);
		
		clac = rintintin_evaluate_constraints(tensor, &p, 0L);
	
		if(check_classification(tensor, &p, clac))
		{
			break;
		}
	}
	
	if(iterations >= config->max_iterations)
	{
		clac = rintintin_evaluate_constraints(tensor, &p, &constraints);
		
		return (rintintin_constrained_result){
			.solution = p,
			.residual = constraints.x + constraints.y + constraints.z,
			.converged = false,
			.improved = true,
			.error_code = -1,
			.calc=(char)clac 
		};	
	}
	
	double mid;
	while(min + (double)config->alpha < max)
	{
		mid = (min + max) * 0.5;
		p = vec3_scale(&dir, mid);
		p = vec3_add(&p, centroid);
				
		clac = rintintin_evaluate_constraints(tensor, &p, &constraints);
		
		if(clac)
		{
			max = mid;
		}
		else
		{
			min = mid;
		}
	}
	
	if(max != mid)
	{
		clac = rintintin_evaluate_constraints(tensor, &p, &constraints);
	}
		
	return (rintintin_constrained_result){
		.solution = p,
		.residual = constraints.x + constraints.y + constraints.z,
		.converged = true,
		.improved = true,
		.error_code = -1,
		.calc=(char)clac  
	};	
}


static error_code flip_w(liquid_t const* tensor, rintintin_vec3 * p)
{
	double mass = dot(&tensor->cubic.o2, p) + tensor->mass_o;
	double o2_norm_sq = dot(&tensor->cubic.o2, &tensor->cubic.o2);
	if(o2_norm_sq < 1e-36) return -2;
	
	double t = (-mass*2) / o2_norm_sq; 
	
	dvec3 o2 = vec3_scale(&tensor->cubic.o2, t);
	*p = vec3_add(p, &o2);
	
	return 1;
}


static rintintin_constrained_result rintintin_attractor_solver(rintintin_liquified const* tensor, rintintin_vec3 const* guess, rintintin_constrained_config const* config, int);

rintintin_constrained_result rintintin_solve_center(rintintin_liquified const* tensor, rintintin_vec3 const* guess, rintintin_constrained_config const* config)
{		
	dvec3 p = *guess;
		
	rintintin_constrained_result try0 = rintintin_attractor_solver(tensor, &p, config, 0);
	if(try0.converged)
	{
		return try0;
	}
	
	dvec3 n_p = try0.solution;
	flip_w(tensor, &n_p);
	
	rintintin_constrained_result try1 = rintintin_attractor_solver(tensor, &n_p, config, 1);
	if(try1.converged)
		return try1;
		
	dvec4 planes[3];
	rintintin_liquified const* it = tensor;
	planes[0] = (dvec4){
		2*it->cubic.o1.x_xx - 4*it->cubic.o2.x,
		(-it->cubic.o1.x_xy + 2*it->cubic.o2.y),
		(-it->cubic.o1.x_xz + 2*it->cubic.o2.z),	
		 -it->linear.coupling.x.x  + 2*it->mass_o,
	};
	planes[1] = (dvec4){
		(-it->cubic.o1.y_xy + 2*it->cubic.o2.x),
		(2*it->cubic.o1.y_yy - 4*it->cubic.o2.y),
		(-it->cubic.o1.y_yz +2*it->cubic.o2.z),	
		 -it->linear.coupling.y.y + 2*it->mass_o,	
	};
	planes[2] = (dvec4){
		(-it->cubic.o1.z_xz + 2*it->cubic.o2.x),
		(-it->cubic.o1.z_yz + 2*it->cubic.o2.y),
		(2*it->cubic.o1.z_zz - 4*it->cubic.o2.z),	
		 -it->linear.coupling.z.z + 2*it->mass_o	
	};
	
	rintintin_constrained_result best = (try0.residual < try1.residual)? try0 : try1;
	
	for(int i = 0; i < 3; ++i)
	{
		double m = dot((dvec3*)&planes[i].x, &try0.solution) + planes[i].w;
		double norm_sq = dot((dvec3*)&planes[i].x, (dvec3*)&planes[i].x);
		if(norm_sq < 1e-36) continue;
	
		for(int j = -1; j < 2; ++j)
		{
			double t = (-m + j) / norm_sq; 
		
			dvec3 o2 = vec3_scale((dvec3*)&planes[i].x, t);
			p = vec3_add(&try0.solution, &o2);
			
			rintintin_constrained_result try = rintintin_attractor_solver(tensor, &p, config, 0);
			if(try.converged)
			{
				return try0;
			}
			
			if(try.residual < best.residual)
			{
				best = try;
			}
		}	
	}
	
	
	
	return best;
}

static rintintin_constrained_result rintintin_attractor_solver(rintintin_liquified const* tensor, rintintin_vec3 const* guess, rintintin_constrained_config const* config, int flip_around)
{
	int iteration;
	const double epsilon = 5e-8;
	
	dvec3 p = *guess;
	dvec3 n_p = p;
	dvec3 c = rintintin_get_center(tensor, &n_p);	
	dvec3 dir = (dvec3){(n_p.x - c.x), (n_p.y - c.y), (n_p.z - c.z)};
	double residual = dot(&dir, &dir);
	double og_mass = tensor->mass_o + dot(&p, &tensor->cubic.o2);
	bool improved = false;
	
	if(residual < epsilon)
		return (rintintin_constrained_result){
			.solution = p,
			.residual = residual,
			.converged = true,
			.improved = false,
			.error_code = 0,  // denominator 0.
			.calc=0
		};	
	
	residual += 0.01;
		
	for(iteration = 0; iteration < config->max_iterations; ++iteration)
	{
		c = rintintin_get_center(tensor, &n_p);
		double mass = tensor->mass_o + dot(&c, &tensor->cubic.o2);
		dvec3 new_dir = (dvec3){(n_p.x - c.x), (n_p.y - c.y), (n_p.z - c.z)};
		double new_residual = dot(&new_dir, &new_dir);
		
		if(flip_around && mass * og_mass < 0)
		{
			flip_w(tensor, &c);
			continue;
		}
		
		if(new_residual >= residual)
		{
			n_p = c;
			continue;
		}
		
	// if log convergence then don't count towards iteration counter	
		if((new_residual / residual) < 0.5)
			--iteration;
			
		p = n_p;
		n_p = c;
		improved = true;
		residual = new_residual;
		
		if(residual < epsilon)
		{
			return (rintintin_constrained_result){
				.solution = p,
				.residual = residual,
				.converged = true,
				.improved = true,
				.error_code = 0,  // denominator 0.
				.calc=0
			};	
		}
	}
	
	return (rintintin_constrained_result){
		.solution = p,
		.residual = residual,
		.converged = false,
		.improved = improved,
		.error_code = 0,  // denominator 0.
		.calc=0
	};	
}
