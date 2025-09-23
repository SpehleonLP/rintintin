
#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_solver.h"
#include "rintintin_constraint_solver.h"
#include "rintintin_math.h"
#include "rintintin_type_hinting.h"

rintintin_mat3x3 rintintin_mat3x3_multiply(const rintintin_mat3x3* a, const rintintin_mat3x3* b);
dvec3 rintintin_mat3x3_mul_vec(const rintintin_mat3x3* mat, const rintintin_vec3* vec) ;
dvec3 rintintin_mat3x3_T_mul_vec(const rintintin_mat3x3* mat, const rintintin_vec3* vec) ;

#ifdef RINTINTIN_VERBOSE 
#include <stdio.h>
#define LOG(...) fprintf(stderr, __VA_ARGS__);
#include <assert.h>
#else
#define LOG(...)
#endif

// don't include *anything*
#ifdef RINTINTIN_USE_STRING
#include <string.h>
#else
void rintintin_memset(void * ptr, int value, uint64_t num);
#define memset rintintin_memset
#endif

#define pXYZ(x) ((rintintin_vec3*)&(x))
#define rXYZ(x) *((rintintin_vec3*)&(x))
	
rintintin_error_code rintintin_compute_child_table(struct rintintin_scratch_space * scratch, int32_t const* parents, uint32_t * counters)
{
	memset(counters, 0, sizeof(uint32_t)*scratch->no_joints);
	uint32_t * children = counters + scratch->no_joints;
	
	int32_t root = -1;
	
	for(uint64_t i = 0; i < scratch->no_joints; ++i)
	{
		int32_t p = parents[i];
		
		if(p < 0)
		{
			if(root >= 0)
				return RINTINTIN_ERROR_TOO_MANY_ROOTS;
			root = (int32_t)i;
			continue;
		}
		
		if((uint32_t)p >= scratch->no_joints)
			return RINTINTIN_ERROR_OUT_OF_BOUNDS;
			
		children[((uint32_t)p)*scratch->no_joints + counters[p]++] = (uint32_t)i;
	}
	
	return root;
}

#define NEXT_INSTRUCTION() goto *dispatch_table[*(uint8_t*)bc]

#define ALIGN(it) it = ((void*)((((intptr_t)(it)) + 63) & (intptr_t)~63))

#define ALLOC(type, name, length) \
	type * name = (type*)(tmp_begin); tmp_begin = (&name[length]); if(tmp_begin > tmp_end) return RINTINTIN_ERROR_OUT_OF_BOUNDS; 
#define M_ALLOC(type, name, length) \
	memo.name = (type*)(tmp_begin); tmp_begin = (&memo.name[length]); if(tmp_begin > tmp_end) return RINTINTIN_ERROR_OUT_OF_BOUNDS; 
	
struct Memo
{
	dvec3 * centroid;
	dvec3 * allow;
	dvec3 * joints;
	uint32_t * counters;
	uint32_t * children;
	tensor_t * subtree;
	void * tmp_begin;
	void * tmp_end;
	int root;
};

static int rintintin_solve_blame(rintintin_command * cmd, struct Memo * memo);
static int rintintin_build_subtree(rintintin_command * cmd, struct Memo *memo);

rintintin_error_code rintintin_solve(rintintin_command * cmd, struct rintintin_metrics * dst, void * tmp_begin, void * tmp_end)
{
	struct rintintin_scratch_space* scratch = (struct rintintin_scratch_space*)cmd->scratch_space;
	uint32_t N =  scratch->no_joints;
	
	struct Memo memo;
	
	assert((void*)(dst+N) <= (void*)tmp_begin);
		
	M_ALLOC(dvec3, centroid, N)
	M_ALLOC(dvec3, allow, N)
	memo.joints = 0L;
	M_ALLOC(uint32_t, counters, N)
	M_ALLOC(uint32_t, children, N*N)
	M_ALLOC(tensor_t, subtree, N)
	memo.tmp_begin = tmp_begin;
	memo.tmp_end = tmp_end;
	
	void * child_end = &memo.children[N*N];
	if((void*)memo.subtree < child_end)
		return RINTINTIN_ERROR_INVALID_INDEX;
	
	memo.root = rintintin_compute_child_table(scratch, cmd->skin.parents, memo.counters);
	if(memo.root < 0) return memo.root;	

	int root = rintintin_solve_blame(cmd, &memo);
	if(root < 0) return root;
		
	int ec = rintintin_build_subtree(cmd, &memo);
	if(ec < 0) return ec;
		
// repurpose memory	
	memo.joints = memo.allow;
	memo.allow = 0L;
			
	rintintin_constrained_config config = rintintin_constrained_default_config(memo.centroid, cmd->skin.joint_translation_mesh_space);
	config.verbose = false;
	
// finish	
	for(int32_t idx = (int32_t)N-1; idx >= 0; --idx)
	{
		uint32_t j = scratch->joints? scratch->joints[idx] : (uint32_t)(idx);

		rintintin_metrics     * acc_j = &dst[j];
		assert((void*)(acc_j+1) <= (void*)memo.centroid);
		rintintin_tensor const* subtree_j = &memo.subtree[j];
		
		config.name = cmd->skin.bone_names? cmd->skin.bone_names[j] : "<unnamed>";
		config.idx = idx;
		config.parent = cmd->skin.parents[j];
		rintintin_constrained_result r = rintintin_walk_to_definite_matrix(subtree_j, &cmd->skin.joint_translation_mesh_space[j], &config);
	
		rintintin_solidified as_solid = rintintin_solidify(subtree_j, &r.solution);

		acc_j->volume = as_solid.mass_o;
		acc_j->inertia = as_solid.constant;
		
		acc_j->centroid.x = as_solid.linear_o.x * 0.5;
		acc_j->centroid.y = as_solid.linear_o.y * 0.5;
		acc_j->centroid.z = as_solid.linear_o.z * 0.5;
	}
	
	return RINTINTIN_SUCCESS;
}


static int rintintin_solve_blame(rintintin_command * cmd, struct Memo * memo)
{
	struct rintintin_scratch_space* scratch = (struct rintintin_scratch_space*)cmd->scratch_space;
	uint32_t N =  scratch->no_joints;
	
	void * tmp_begin = memo->tmp_begin;
	void * tmp_end = memo->tmp_end;
	
	dvec3 * centroid = memo->centroid;	
	ALLOC(liquid_t, liquid, N)
	
// initialize centroids
	if(cmd->skin.joint_translation_mesh_space)
		rintintin_memcpy(centroid, cmd->skin.joint_translation_mesh_space, N * sizeof(dvec3));
	else
	{
		for(uint32_t idx = 0; idx < N; ++idx)
		{
			dvec3 sum = vec3_add(&scratch->aabbs[idx].max, &scratch->aabbs[idx].min);
			centroid[idx] = vec3_scale(&sum, 0.5);		
		}	
	}
	
// initialize allow	
	for(uint32_t i = 0; i < N; ++i)
		memo->allow[i] = (dvec3){1.0, 1.0, 1.0};
		
// compute centroids per local mesh.
	rintintin_constrained_config config = rintintin_constrained_default_config(memo->centroid, cmd->skin.joint_translation_mesh_space);
	config.verbose = false;
	
	{
		void * savestate = tmp_begin;
		ALLOC(liquid_t, subtree, N)
		ALLOC(dvec4, subtree_centroid, N)
		memset(subtree, 0, sizeof(*subtree)*N);
		memset(subtree_centroid, 0, sizeof(*subtree_centroid)*N);
		tmp_begin = savestate;
		
		for(uint32_t idx = 0; idx < N; ++idx)
		{
			uint32_t j = scratch->joints? scratch->joints[(N-1)-idx] : (uint32_t)((N-1)-idx);
			int32_t p = cmd->skin.parents[j];
			
			tensor_t solid = rintintin_tensor_from_coeff(scratch->latent[j].thick);
			liquid[j] = (liquid_t){
				.cubic={.o1=solid.cubic.o1, .o2=solid.cubic.o2},
				.linear={.coupling=solid.linear.coupling, .o=solid.linear.o},
				.mass_o=solid.mass_o
			};
			
			liquid_add_eq(&subtree[j], &liquid[j]);
			if(p >= 0)	liquid_add_eq(&subtree[p], &subtree[j]);
			
			config.name = cmd->skin.bone_names? cmd->skin.bone_names[j] : "<unnamed>";
			config.idx = (int32_t)idx;
			config.parent = p;
			rintintin_constrained_result result = rintintin_solve_center(&liquid[j], &centroid[j], &config);
		
			if(!result.improved && !result.converged)
				continue;
			
			if(subtree_centroid[j].w != 0)
			{
				int break_point = 0;
				++break_point;
			}
						
			centroid[j] = result.solution;	
			continue;
			double mass = subtree[j].mass_o + dot(&subtree[j].cubic.o2, &centroid[j]);
			double invMass = mass - subtree_centroid[j].w;
			invMass = invMass? 1.0 / invMass : 0.0;
			
			if(p >= 0)
			{								
				subtree_centroid[p] = (dvec4){
					.x=subtree_centroid[p].x + centroid[j].x * mass,
					.y=subtree_centroid[p].y + centroid[j].y * mass,
					.z=subtree_centroid[p].z + centroid[j].z * mass,
					.w=subtree_centroid[p].w + mass,
				};
			}
			
			if(invMass)
			{
				centroid[j] = (dvec3){
					.x=(centroid[j].x * mass - subtree_centroid[j].x) * invMass,
					.y=(centroid[j].y * mass - subtree_centroid[j].y) * invMass,
					.z=(centroid[j].z * mass - subtree_centroid[j].z) * invMass,
				};
			}
		}
	}
	
	return RINTINTIN_SUCCESS;	
}

static int rintintin_build_subtree(rintintin_command * cmd, struct Memo * memo)
{
	uint32_t idx;
		
	struct rintintin_scratch_space* scratch = (struct rintintin_scratch_space*)cmd->scratch_space;
	uint32_t N =  scratch->no_joints;
	
	void * tmp_begin = memo->tmp_begin;
	void * tmp_end = memo->tmp_end;
	ALLOC(tensor_t, i_func, N)
	ALLOC(tensor_t, error, N)
		
	for(idx = 0; idx < N; ++idx)
	{
		latent_t * latent_j = &scratch->latent[idx];		
		i_func[idx] = rintintin_tensor_from_coeff(latent_j->thick);
		
		if(latent_j->thin[THIN_MASS] != 0.0)
		{
			rintintin_add_thin_shell(&i_func[idx] , latent_j->thin);
		}
	}	
	
	{
		for(idx = 0; idx < N; ++idx) { memo->subtree[idx] = i_func[idx]; }
	
		for(idx = 0; idx < N; ++idx)
		{
			uint32_t j = scratch->joints? scratch->joints[(N-1)-idx] : (uint32_t)((N-1)-idx);
			int32_t p = cmd->skin.parents[j];
			
			if(p >= 0)
				memo->subtree[p] = rintintin_tensor_add(&memo->subtree[p], &memo->subtree[j]); 
		}
	}


	return RINTINTIN_SUCCESS;
}

