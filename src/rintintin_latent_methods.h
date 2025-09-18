#ifndef rintintin_TENSOR_METHODS_H
#define rintintin_TENSOR_METHODS_H
#include "rintintin_latentspace.h"
#include "rintintin_type_hinting.h"

// this is the common header of all the *internal* rintintin parts.

rintintin_tensor rintintin_tensor_add(rintintin_tensor const* a, rintintin_tensor const* b);
rintintin_tensor rintintin_tensor_sub(rintintin_tensor const* a, rintintin_tensor const* b);

void rintintin_tensor_add_eq(rintintin_tensor * a, rintintin_tensor const* b);
void rintintin_tensor_sub_eq(rintintin_tensor * a, rintintin_tensor const* b);
void rintintin_liquid_add_eq(rintintin_liquified * a, rintintin_liquified const* b);
void rintintin_liquid_sub_eq(rintintin_liquified * a, rintintin_liquified const* b);

rintintin_solidified rintintin_solidify(rintintin_tensor const* m, rintintin_vec3 const* p);
void rintintin_solid_add_eq(rintintin_solidified * a, rintintin_solidified const* b);
void rintintin_tensor_sub_eq_solid(rintintin_tensor * a, rintintin_solidified const* b);

void rintintin_add_thin_shell(rintintin_tensor * dst, double const* thin);
rintintin_metrics rintintin_compute_metrics_from_solid(rintintin_solidified const* a);
void rintintin_mat3x3_add_eq(rintintin_mat3x3 *RESTRICT a, rintintin_mat3x3 const*RESTRICT b);
rintintin_vec3 rintintin_evaluate_centroid(rintintin_tensor const * it, rintintin_vec3 const* p);
rintintin_tensor rintintin_tensor_heuristic(rintintin_tensor const* base, rintintin_tensor const* to_remove, rintintin_vec3 const* p, int iteration);

/* these are the constraints i need all of them to be > 0 for a specific p!
 * 
 * meaning we can't use a log solver b/c they start out as negative and that means undefined as a log! 
 * constraint.w = dot({vec3 p, 1}, v_func)
 * general form of the other three is something like 
 * 
 * constraints = (o₂·p + mass_o) * (C + L·p + Q·p² + R·p³) - (L'·p + C'·p²)² 
*/

rintintin_vec3 rintintin_get_center(rintintin_liquified const * it, rintintin_vec3 const* p);
rintintin_liquified rintintin_liquid_less(rintintin_liquified const * it, rintintin_vec3 const* p);

#if USING_NAMESPACE_RINTINTIN
typedef rintintin_error_code error_code;
typedef rintintin_symmetric_mat3 sdmat3;
#define evaluate_jacobian rintintin_evaluate_jacobian

#define solve_constrained rintintin_solve_constrained
#define liquid_add_eq rintintin_liquid_add_eq
#define liquid_sub_eq rintintin_liquid_sub_eq
#endif

#endif // rintintin_TENSOR_METHODS_H
