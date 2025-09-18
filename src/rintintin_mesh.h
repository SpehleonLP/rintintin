#ifndef RINTINTIN_MESH_H
#define RINTINTIN_MESH_H
#include "../include/rintintin.h"


struct rintintin_vertex
{
	rintintin_vec4 position;
	rintintin_vec4 weight;
	int joint[4];
};

typedef rintintin_error_code rintintin_tri_callback(rintintin_vec3 const* verts, double const* weights, int32_t joint, void * user_data);
typedef rintintin_error_code rintintin_vert_callback(struct rintintin_vertex const* vert, void * user_data);
rintintin_error_code rintintin_process_mesh(rintintin_mesh const* src, uint32_t thread_id, uint32_t no_threads, rintintin_tri_callback * cb, void * userdata);
rintintin_error_code rintintin_visit_each_index(rintintin_mesh const* src, uint32_t thread_id, uint32_t no_threads, rintintin_vert_callback * cb, void * userdata);

#endif // RINTINTIN_MESH_H
