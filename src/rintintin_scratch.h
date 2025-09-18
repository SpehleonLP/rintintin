#ifndef RINTINTIN_PREPARE_H
#define RINTINTIN_PREPARE_H
#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_latentspace.h"
#include "rintintin_latent_methods.h"


typedef enum 
{
	STAGE_PREPARED = 1,
	STAGE_CREATE_LATENT_SPACE,
	STAGE_PARALLEL_REDUCTION,
	STAGE_FINISH
} rintintin_stage;

struct rintintin_scratch_space
{
	uint64_t magic;
 
	uint32_t       * joints;
	rintintin_aabb * aabbs;
	rintintin_latent_space * latent;
	void * workspace;
			
	uint32_t max_threads;
	uint32_t no_joints;
	rintintin_stage stage_id;
	uint32_t last_stage_threads_used;
	uint32_t flat_latent_size;
	uint32_t flat_tensor_size;
	uint64_t byte_length;
};


void rintintin_memset(void *ptr, int value, uint64_t num);
void rintintin_memcpy(void *RESTRICT dst, void const*RESTRICT src, uint64_t num);


#endif // RINTINTIN_PREPARE_H
