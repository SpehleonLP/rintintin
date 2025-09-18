#include "rintintin_scratch.h"
#include "../include/rintintin.h"
#include "rintintin_eigen.h"
#include "rintintin_solver.h"
#include "rintintin_latent_methods.h"
#include "rintintin_type_hinting.h"
#include "rintintin_mesh.h"
#include <assert.h>
#include <stddef.h>

#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_math.h"

/// ISO C restricts enumerator value to range of int before C2X
#define MAGIC_NUMBER 0xDEADBEEF

#ifdef RINTINTIN_USE_STRING
#include <string.h>
#else
// don't include *anything*
#define NEED_TO_MAKE_STRING
void rintintin_memset(void * ptr, int value, uint64_t num);
#define memset rintintin_memset
#endif

uint64_t rintintin_get_scratch_space_size(unsigned int no_joints, unsigned int max_threads)
{
	uint64_t N = no_joints;
	
	// Base scratch space structure
	uint64_t base_size = sizeof(struct rintintin_scratch_space);
	
	// Joints array (with padding to even boundary)
	uint64_t joints_size = sizeof(uint32_t) * ((N + 1) & 0xFFFFFFFE);
	
	// Latent space (cache line aligned)
	uint64_t flat_latent_size = N * sizeof(rintintin_latent_space);
	flat_latent_size = (flat_latent_size + 63) & ~63UL;
	
	// Latent space (cache line aligned)
	uint64_t flat_aabb_size = N * sizeof(dvec3);
	flat_aabb_size = (flat_aabb_size + 63) & ~63UL;
	
	// Fixed solver data section
	uint64_t solver_size = sizeof(dvec3) * N;
	uint64_t accumulator_size = sizeof(rintintin_metrics) * N;
	uint64_t counters_size = sizeof(uint32_t) * N;
	uint64_t children_size = sizeof(uint32_t) * N * N;
	uint64_t fixed_solver_total = solver_size + accumulator_size + counters_size + children_size;
	
	// Cache line align after solver data
	fixed_solver_total = (fixed_solver_total + 63) & ~63UL;
	
	// Dynamic workspace sizes (these use the tmp pointer after solver data)
	uint64_t flat_vec4_size = N * sizeof(rintintin_vec4);
	flat_vec4_size = (flat_vec4_size + 63) & ~63UL;  // Fix: was using flat_latent_size
	
	uint64_t flat_tensor_size = N * sizeof(rintintin_tensor);
	flat_tensor_size = (flat_tensor_size + 63) & ~63UL;
	
	uint64_t flat_solidified_size = N * sizeof(rintintin_solidified);
	flat_solidified_size = (flat_solidified_size + 63) & ~63UL;
	
	// Case 1: 4 vec4 arrays (v_func, work_space, subtree, to_distribute)
	uint64_t case1_size = flat_vec4_size * 4;
	
	// Case 2: tensor array + solidified array
	uint64_t case2_size = flat_tensor_size + flat_solidified_size;
	
	// Other potential workspace sizes (keeping your original calculations but fixing names)
	uint64_t preprocess_mesh_size = flat_latent_size * max_threads;
	uint64_t preprocessed_mesh_size = flat_latent_size + flat_vec4_size * 2;
	uint64_t repair_size = flat_latent_size + flat_vec4_size * 4 + N * (1 + N) * sizeof(uint32_t);
	uint64_t process_mesh_size = flat_tensor_size * max_threads + flat_latent_size;
	uint64_t finish_size = flat_latent_size + flat_tensor_size * 2 + sizeof(rintintin_metrics) * N * 2;
	
	// Find maximum workspace size needed
	uint64_t max_workspace = case1_size;
	if (case2_size > max_workspace) max_workspace = case2_size;
	if (preprocess_mesh_size > max_workspace) max_workspace = preprocess_mesh_size;
	if (preprocessed_mesh_size > max_workspace) max_workspace = preprocessed_mesh_size;
	if (repair_size > max_workspace) max_workspace = repair_size;
	if (process_mesh_size > max_workspace) max_workspace = process_mesh_size;
	if (finish_size > max_workspace) max_workspace = finish_size;
	
	// Total size calculation
	uint64_t total_size = flat_tensor_size + base_size + flat_aabb_size + joints_size + flat_latent_size + fixed_solver_total + max_workspace;
	
	// Add padding for safety
	total_size += 128;
	
	return total_size;
}

static rintintin_error_code rintintin_OPENING_CHECKS(rintintin_command * cmd)
{
	if(cmd == 0L
	|| cmd->scratch_space == 0L
	|| cmd->skin.parents == 0L
	//|| cmd->skin.joint_translation_mesh_space == 0L
	)
		return RINTINTIN_ERROR_NULL_POINTER;
		
	if(cmd->scratch_space_byte_length < rintintin_get_scratch_space_size(cmd->skin.no_joints, cmd->max_threads))
		return RINTINTIN_ERROR_OUT_OF_BOUNDS;

	return RINTINTIN_SUCCESS;
}

static rintintin_error_code rintintin_STAGE_CHECKS(rintintin_command * cmd, struct rintintin_scratch_space * scratch, uint32_t stage_id)
{
	if(scratch->magic != MAGIC_NUMBER) return RINTINTIN_ERROR_NOT_PREPARED;
	if(scratch->stage_id > stage_id)	return RINTINTIN_ERROR_OUT_OF_ORDER;
	if (scratch->max_threads !=(cmd->max_threads? cmd->max_threads : 1)
	|| scratch->no_joints !=	cmd->skin.no_joints)
		return RINTINTIN_ERROR_CONFIG_CHANGED;
		
	return RINTINTIN_SUCCESS;	
}

#define OPENING_CHECKS\
	{ int ec = rintintin_OPENING_CHECKS(cmd); if(ec != RINTINTIN_SUCCESS) return ec; };\
	struct rintintin_scratch_space* scratch = (struct rintintin_scratch_space*)cmd->scratch_space;
	
#define STAGE_CHECKS(x)\
	OPENING_CHECKS\
	{ int ec = rintintin_STAGE_CHECKS(cmd, scratch, x); if(ec != RINTINTIN_SUCCESS) return ec; };\
		
#define ROUND_TO_64(x)  (((int64_t)(x) + (int64_t)63) & ~(int64_t)63);

rintintin_error_code rintintin_begin(rintintin_command * cmd)
{
	OPENING_CHECKS
	
	scratch->max_threads = cmd->max_threads > 0? cmd->max_threads : 1;
	scratch->no_joints = cmd->skin.no_joints;
	scratch->stage_id = STAGE_PREPARED;
	scratch->last_stage_threads_used = 1;
	scratch->magic = MAGIC_NUMBER;
	scratch->byte_length = cmd->scratch_space_byte_length;
			
	scratch->flat_latent_size = (uint32_t)ROUND_TO_64(scratch->no_joints * sizeof(latent_t));
	scratch->flat_tensor_size = (uint32_t)ROUND_TO_64(scratch->no_joints * sizeof(tensor_t));
	
	scratch->aabbs  = (rintintin_aabb*)(scratch+1);
	scratch->joints = (uint32_t*)(scratch->aabbs+scratch->no_joints); 

	intptr_t ptr = (intptr_t)(&scratch->joints[scratch->no_joints]);
	
// round to cache line	
	ptr = ROUND_TO_64(ptr);
	scratch->latent = (latent_t*)ptr;
	
	scratch->workspace = (uint8_t*)scratch->latent + scratch->flat_latent_size;
			
	int8_t need_topological_sort = 0;

	for(uint64_t i = 0; i < scratch->no_joints; ++i)
	{
		scratch->joints[i] = (uint32_t)i;
		
		scratch->aabbs[i].min = (dvec3){2e200, 2e200, 2e200};
		scratch->aabbs[i].max = (dvec3){-2e200, -2e200, -2e200};
		
		if(cmd->skin.parents[i] == (int)(i))
			return RINTINTIN_ERROR_INVALID_INDEX;
			
		if(cmd->skin.parents[i] > (int)(i))
			need_topological_sort = 1;	
	}
		
	if(need_topological_sort)
	{
		uint32_t * counters = (uint32_t*)(scratch->latent);
		uint32_t * children = counters + scratch->no_joints;
		int32_t root = rintintin_compute_child_table(scratch, cmd->skin.parents, counters);
		
		if(root < 0)
			return root;
		
		uint32_t index = 0;
		uint32_t stack_size = 1;
		uint32_t * stack = children + scratch->no_joints*scratch->no_joints;
		stack[0] = (uint32_t)root;
		
		while(stack_size)
		{
			uint32_t top = stack[--stack_size];
			scratch->joints[index] = top;
			
			for(uint64_t i = 0; i < counters[top]; ++i)
			{
				stack[stack_size++] = children[top*scratch->no_joints + i];
			}
		}	
	}
		
	memset(scratch->latent, 0, scratch->flat_latent_size * scratch->max_threads);	

// DEBUG
//	cmd->skin.joint_translation_mesh_space[0] = cmd->skin.joint_translation_mesh_space[1];

	return RINTINTIN_SUCCESS;
}

struct stage_2_cb_data
{
	rintintin_latent_space* write;
	rintintin_aabb * aabbs;
	uint32_t no_joints;
	int is_thin_shell;
	double thickness;
};

static rintintin_error_code rintintin_tri_callback_stage_2(dvec3 const* verts, double const* weights, int32_t joint, void * user_data)
{
	struct stage_2_cb_data * scratch =  ((struct stage_2_cb_data*)user_data);
	
	if((uint32_t)joint >= scratch->no_joints)
		return RINTINTIN_ERROR_OUT_OF_BOUNDS;
		
	struct rintintin_input tri = 
    {
		.p0=verts[0],
		.p1=verts[1],
		.p2=verts[2],
		.d0=weights[0],
		.d1=weights[1],
		.d2=weights[2],	
    };
    
    if(scratch->is_thin_shell == 0)
		rintintin_thick_construct(scratch->write[joint].thick, &tri);
	else
		rintintin_thin_construct(scratch->write[joint].thin, &tri, scratch->thickness);

	if(scratch->aabbs)
	{
		for(int i = 0; i < 3; ++i)
		{
			scratch->aabbs[joint].min = vmin(&scratch->aabbs[joint].min, &verts[i]);
			scratch->aabbs[joint].max = vmax(&scratch->aabbs[joint].max, &verts[i]);
		}
	}
	
	return RINTINTIN_SUCCESS;
}


static rintintin_error_code rintintin_vert_callback_stage_2(struct rintintin_vertex const* vert, void * user_data)
{
	struct stage_2_cb_data * scratch =  ((struct stage_2_cb_data*)user_data);
	
	for(int i = 0; i < 4; ++i)
	{
		if((&vert->weight.x)[i] <= 0) continue;
		if((uint32_t)vert->joint[i] >= scratch->no_joints) continue;	
		
		int j = vert->joint[i];
		scratch->aabbs[j].min = vmin(&scratch->aabbs[j].min, (rintintin_vec3*)&vert->position);
		scratch->aabbs[j].max = vmax(&scratch->aabbs[j].max, (rintintin_vec3*)&vert->position);
	}

	return RINTINTIN_SUCCESS;
}

rintintin_error_code rintintin_non_manifold_repair(rintintin_command * cmd);

rintintin_error_code rintintin_read_mesh(rintintin_command * cmd, uint32_t thread_id, uint32_t no_threads_this_stage)
{
	STAGE_CHECKS(STAGE_CREATE_LATENT_SPACE)
			
	no_threads_this_stage = no_threads_this_stage > 1? no_threads_this_stage-1 : 1;
		
	if(scratch->stage_id == STAGE_PREPARED)
	{
		scratch->stage_id = STAGE_CREATE_LATENT_SPACE;
		scratch->last_stage_threads_used = no_threads_this_stage;
	}	
	
	
	struct stage_2_cb_data userdata;
	
	userdata.aabbs = 0L;
	userdata.write = (rintintin_latent_space*)((uint8_t*)scratch->latent + scratch->flat_latent_size * thread_id);
	userdata.no_joints = scratch->no_joints;
		
	// solve AABB
	if(thread_id == no_threads_this_stage)
	{
		userdata.aabbs = scratch->aabbs;
	
		if(thread_id != 0)
		{	
			for(uint64_t i = 0; i < cmd->no_meshes; ++i)
			{
				int ec = rintintin_visit_each_index(&cmd->meshes[i], 0, 1, rintintin_vert_callback_stage_2, &userdata);
				
				if(ec != RINTINTIN_SUCCESS)
					return ec;
			}
		
			return RINTINTIN_SUCCESS;
		}
	}
	
	for(uint64_t i = 0; i < cmd->no_meshes; ++i)
	{
		userdata.thickness = cmd->meshes[i].thickness;
		userdata.is_thin_shell = (cmd->meshes[i].surface_mode == RINTINTIN_SURFACE_THIN_SHELL);
		int ec = rintintin_process_mesh(&cmd->meshes[i], thread_id, no_threads_this_stage, rintintin_tri_callback_stage_2, &userdata);
		
		if(ec != RINTINTIN_SUCCESS)
			return ec;
	}
	
	return RINTINTIN_SUCCESS;
}

static void rintintin_sum_vec4_arrays(double * chunk, uint32_t chunk_size_bytes, uint32_t no_joints, uint32_t a, uint32_t b)
{
	uint64_t flat_latent_size = no_joints * chunk_size_bytes;
	flat_latent_size = (flat_latent_size + 63) & ~63UL;
	
// tell the compiler it can use SIMD here.
	double *RESTRICT A = (double*)((char*)(chunk) + flat_latent_size*a);
	double *RESTRICT B = (double*)((char*)(chunk) + flat_latent_size*b);
	
    ASSUME_ALIGNED(A, 64);
    ASSUME_ALIGNED(B, 64);
    
// note: actual size is padded so that its okay to write over the end with SIMD	
	uint32_t N = no_joints * (chunk_size_bytes / sizeof(double));
	
#if defined(__GNUC__) && !defined(__clang__)
    #pragma GCC ivdep
    #pragma GCC unroll 4
#elif defined(__clang__)
    #pragma clang loop vectorize(enable)
    #pragma clang loop unroll_count(4)
#elif defined(_MSC_VER)
    #pragma loop(ivdep)
#endif

	for(uint64_t i = 0u; i < N; ++i)
	{
		A[i] += B[i]; 
	}
}

#include <stdio.h>

rintintin_error_code rintintin_parallel_reduction(rintintin_command * cmd, uint32_t loop, uint32_t thread_id, uint32_t no_threads_this_stage)
{
	STAGE_CHECKS(STAGE_PARALLEL_REDUCTION)
		
	if(no_threads_this_stage < scratch->last_stage_threads_used)
		return RINTINTIN_ERROR_THREAD_COUNT_MISMATCH;

	if(scratch->stage_id == STAGE_CREATE_LATENT_SPACE)
		scratch->stage_id = STAGE_PARALLEL_REDUCTION;
		
	if(scratch->stage_id != STAGE_PARALLEL_REDUCTION)
	{
		return RINTINTIN_ERROR_OUT_OF_ORDER;
	}	

	// Each iteration reduces by half
	uint32_t stride = 1 << loop;  // 1, 2, 4, 8, ...
	uint32_t step = stride * 2;   // 2, 4, 8, 16, ...
        
	// Only some threads participate in each iteration
	uint32_t partner_id = thread_id + stride;
	
	if((thread_id % step != 0) || (partner_id) >= scratch->last_stage_threads_used)
		return RINTINTIN_SUCCESS;
	
	fprintf(stderr, "loop %d: copying from %i to %d\n", loop, partner_id, thread_id);
	
	// Sum partner's data into this thread's data
	rintintin_sum_vec4_arrays(
		(double*)scratch->latent,  
		sizeof(rintintin_latent_space),
		scratch->no_joints,        
		thread_id,                 // destination (a)
		partner_id                 // source (b)
	);

	// done with the partner reset it to 0.
	uint8_t * partner_offset = (uint8_t*)scratch->latent + scratch->flat_latent_size*partner_id;
	if(partner_offset + scratch->flat_latent_size > (uint8_t*)cmd->scratch_space+cmd->scratch_space_byte_length)
		return RINTINTIN_ERROR_OUT_OF_BOUNDS;
		
	memset(partner_offset, 0, scratch->flat_latent_size);
	
	// Check if reduction is complete
    // Done when stride >= half the thread count
    if(stride >= (scratch->last_stage_threads_used+1) / 2)
        return RINTINTIN_SUCCESS;  // Done
    else
        return 1;  // More work needed
}

// internal function
static void rintintin_stage_3_single_reduction(rintintin_command * cmd, struct rintintin_scratch_space * scratch)
{	
	if(scratch->stage_id == STAGE_CREATE_LATENT_SPACE)
	{
		scratch->stage_id = STAGE_PARALLEL_REDUCTION;
		
		for(uint64_t i = 1u; i < scratch->last_stage_threads_used; ++i)
		{
			rintintin_sum_vec4_arrays(
				(double*)scratch->latent,  // your chunk pointer
				sizeof(rintintin_latent_space),
				scratch->no_joints,       // no_joints
				0,                 // destination (a)
				(uint32_t)i        // source (b)
			);
		}
		
		// reset memory to 0, done with it all. 		
		void * begin = &scratch->latent[scratch->no_joints];
		memset(begin, 0, scratch->flat_latent_size*(scratch->last_stage_threads_used-1));
	}
	
	for(uint32_t i = 0; i < cmd->skin.no_joints; ++i)
	{
		if(scratch->aabbs[i].min.x > scratch->aabbs[i].max.x)
		{
			if( cmd->skin.joint_translation_mesh_space)
			{
				scratch->aabbs[i].min = 
				scratch->aabbs[i].max = cmd->skin.joint_translation_mesh_space[i];
			}
			else
			{
				scratch->aabbs[i].min = 
				scratch->aabbs[i].max = (rintintin_vec3){0.0, 0.0, 0.0};
			}
		}
	}
}

rintintin_error_code rintintin_end(rintintin_command * cmd)
{
	STAGE_CHECKS(STAGE_FINISH)
	
	if(cmd->results == 0L)
		return RINTINTIN_ERROR_NULL_POINTER;
		
	rintintin_stage_3_single_reduction(cmd, scratch);
	
	if(scratch->stage_id != STAGE_PARALLEL_REDUCTION)
		return RINTINTIN_ERROR_OUT_OF_ORDER;
	
	scratch->stage_id = STAGE_FINISH;
	
	rintintin_error_code ec;
	
	// set up solver;
	uint32_t N =  scratch->no_joints;
	void * begin = ((uint8_t*)(scratch->latent) + scratch->flat_latent_size);
	rintintin_metrics * accumulator = (rintintin_metrics*)(begin);
	
	// bump to cache line
	intptr_t after_solver_int = (intptr_t)(accumulator+N);
	after_solver_int = ROUND_TO_64(after_solver_int);
	void* after_solver = (void*)after_solver_int;
	// how is this failing?
	assert((void*)(accumulator+N) <= (void*)after_solver);
	ec = rintintin_solve(cmd, accumulator, after_solver, (uint8_t*)cmd->scratch_space + cmd->scratch_space_byte_length);
	if(ec < 0) return ec;
	
// ---------------------------------------------------
// FINAL OUTPUT
// ---------------------------------------------------		
	rintintin_metrics * final = cmd->results;
	for(int32_t idx = (int32_t)scratch->no_joints-1; idx >= 0; --idx)
	{
		uint32_t j = scratch->joints? scratch->joints[idx] : (uint32_t)(idx);
		int32_t p = cmd->skin.parents[j];
		
// for debugger view		
		const char * name = cmd->skin.bone_names? cmd->skin.bone_names[j] : "";
		(void)name;
		
		rintintin_metrics * acc_j = &accumulator[j];
		rintintin_metrics * final_j = &final[j];
		
		final_j->volume = acc_j->volume - final_j->volume;
				
		if(final_j->volume == 0)
			final_j->centroid = cmd->skin.joint_translation_mesh_space? 
				cmd->skin.joint_translation_mesh_space[j] 
				: (dvec3){0, 0, 0};
		else
		{
			double invVolume = 1.0 / final_j->volume;
			final_j->centroid.x = (acc_j->centroid.x - final_j->centroid.x) * invVolume;
			final_j->centroid.y = (acc_j->centroid.y - final_j->centroid.y) * invVolume;
			final_j->centroid.z = (acc_j->centroid.z - final_j->centroid.z) * invVolume;
		}
			
		final_j->inertia = smat3_sub(&acc_j->inertia, &final_j->inertia);	
		
		if(p >= 0)
		{
			rintintin_metrics * final_p = &final[p];
		
			final_p->volume += acc_j->volume;
			final_p->centroid.x += acc_j->centroid.x;
			final_p->centroid.y += acc_j->centroid.y;
			final_p->centroid.z += acc_j->centroid.z;
			final_p->inertia = smat3_add(&acc_j->inertia, &final_p->inertia);	
		}	
					
		// parallel axis theorem...
		// (we measured the inertia at the origin so move it to the centroid)				
		smat3 axis = parallel_axis(final_j->volume, &final_j->centroid);
		final_j->inertia = smat3_sub(&final_j->inertia, &axis);	
				
		// flip around
		double xx = final_j->inertia.xx;
		double yy = final_j->inertia.yy;
		double zz = final_j->inertia.zz;
		
		final_j->inertia.xx = yy + zz;
		final_j->inertia.yy = xx + zz;
		final_j->inertia.zz = xx + yy;
	}
	
	return RINTINTIN_SUCCESS;
}

#ifdef NEED_TO_MAKE_STRING

/* Memory functions */
void rintintin_memcpy(void* restrict dest, const void* restrict src, size_t n) {
    char* restrict d = (char* restrict)dest;
    const char* restrict s = (const char* restrict)src;
    
    #ifdef __x86_64__
    /* Use optimized word/dword copy for aligned data, then handle remainder */
    if (n >= 8 && ((uintptr_t)d & 7) == 0 && ((uintptr_t)s & 7) == 0) {
        size_t words = n / 8;
        __asm__ volatile (
            "rep movsq"
            : "=D" (d), "=S" (s), "=c" (words)
            : "0" (d), "1" (s), "2" (words)
            : "memory"
        );
        n &= 7;
    }
    __asm__ volatile (
        "rep movsb"
        : "=D" (d), "=S" (s), "=c" (n)
        : "0" (d), "1" (s), "2" (n)
        : "memory"
    );
    #else
    /* Optimized word copy for aligned data */
    if (n >= sizeof(size_t) && ((uintptr_t)d & (sizeof(size_t)-1)) == 0 && 
        ((uintptr_t)s & (sizeof(size_t)-1)) == 0) {
        size_t* wd = (size_t*)d;
        const size_t* ws = (const size_t*)s;
        size_t words = n / sizeof(size_t);
        
        while (words--) {
            *wd++ = *ws++;
        }
        
        d = (char*)wd;
        s = (const char*)ws;
        n &= (sizeof(size_t) - 1);
    }
    
    while (n--) {
        *d++ = *s++;
    }
    #endif
}


#if defined(__GNUC__) || defined(__clang__)
    #if defined(__x86_64__) || defined(__x86_64) || defined(__amd64__) || defined(__amd64)
		#define MADE_MEMSET
        // x86-64 version
        void rintintin_memset(void *ptr, int value, uint64_t num) {
            __asm__ volatile (
                "rep stosb"
                : "=D" (ptr), "=c" (num)
                : "0" (ptr), "1" (num), "a" (value)
                : "memory"
            );
        }
    #elif defined(__i386__) || defined(__i386) || defined(i386)
		#define MADE_MEMSET
        // x86-32 version (same code works)
        void rintintin_memset(void *ptr, int value, uint64_t num) {
            __asm__ volatile (
                "rep stosb"
                : "=D" (ptr), "=c" (num)
                : "0" (ptr), "1" (num), "a" (value)
                : "memory"
            );
        }
    #endif
#endif

#ifndef MADE_MEMSET

void rintintin_memset(void *ptr, int value, uint64_t num) {
    uint8_t *p = (uint8_t *)ptr;
    uint8_t val = (uint8_t)value;
    
    // Handle unaligned bytes at the beginning
    while (num > 0 && ((uintptr_t)p & (sizeof(uint64_t) - 1)) != 0) {
        *p++ = val;
        num--;
    }
    
    // Fill word-sized chunks
    if (num >= sizeof(uint64_t)) {
        uint64_t word_val = 0;
        for (int i = 0; i < sizeof(uint64_t); i++) {
            word_val = (word_val << 8) | val;
        }
        
        uint64_t *word_ptr = (uint64_t *)p;
        while (num >= sizeof(uint64_t)) {
            *word_ptr++ = word_val;
            num -= sizeof(uint64_t);
        }
        p = (uint8_t *)word_ptr;
    }
    
    // Handle remaining bytes
    while (num > 0) {
        *p++ = val;
        num--;
    }
}
#endif
#endif

// Error message strings array (index 0 = error code -1, index 1 = error code -2, etc.)
static const char* rintintin_error_messages[] = {
    "Invalid argument",                    // -1
    "A required pointer parameter was NULL",                    // -2
    "Index or size parameter was out of valid range",          // -3
    "Invalid data type specified",                              // -4
    "Invalid size parameter",                                   // -5
    "Invalid index parameter",                                  // -6
    "Thread count mismatch",                                    // -7
    "Out of order",                                            // -8
    "Not prepared",                                            // -9
    "Too many roots",                                          // -10
    "Configuration changed"                                     // -11
};

const char* rintintin_get_error_string(int ec)
{
	if(ec >= 0) return 0L;
    int32_t index = -ec - 1;  // Convert negative error code to array index
    int32_t max_index = sizeof(rintintin_error_messages) / sizeof(rintintin_error_messages[0]);
    
    if (index >= 0 && index < max_index) {
        return rintintin_error_messages[index];
    }
    return "Unknown error";
}
