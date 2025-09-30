#include "rintintin_mesh.h"
#include "../include/rintintin.h"

typedef struct rintintin_1 rintintin_tensor_main;

typedef struct rintintin_tri
{
	double       position[3][3];
	float        weight[3][4];
	uint32_t joint[3][4];
} rintintin_tri;


// Helper function to read an index from the index buffer
static uint32_t read_index(const void* index_buffer, rintintin_type index_type,uint64_t index) {
    switch (index_type) {
        case RINTINTIN_TYPE_UNSIGNED_BYTE: {
            const uint8_t* buf = (const uint8_t*)index_buffer;
            return buf[index];
        }
        case RINTINTIN_TYPE_UNSIGNED_SHORT: {
            const uint16_t* buf = (const uint16_t*)index_buffer;
            return buf[index];
        }
        case RINTINTIN_TYPE_UNSIGNED_INT: {
            const uint32_t* buf = (const uint32_t*)index_buffer;
            return buf[index];
        }
        default:
            return 0;
    }
}

static rintintin_error_code rintintin_get_indices(uint64_t * vertex_indices, const rintintin_mesh* src, uint64_t tri_id) 
{
    switch (src->geometry_type) {
        case RINTINTIN_TRIANGLES: {
            // Simple case: each triangle uses 3 consecutive indices
           uint64_t base_index = tri_id * 3;
            if (base_index + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
            
            if (src->index_array_buffer) {
                // Indexed geometry
                vertex_indices[0] = read_index(src->index_array_buffer, src->index_type, base_index);
                vertex_indices[1] = read_index(src->index_array_buffer, src->index_type, base_index + 1);
                vertex_indices[2] = read_index(src->index_array_buffer, src->index_type, base_index + 2);
            } else {
                // Non-indexed geometry
                vertex_indices[0] = base_index;
                vertex_indices[1] = base_index + 1;
                vertex_indices[2] = base_index + 2;
            }
            break;
        }
        
        case RINTINTIN_TRIANGLE_STRIP: {
            // Triangle strip: triangles share vertices
            // Triangle 0: vertices 0,1,2
            // Triangle 1: vertices 1,2,3 (or 1,3,2 for correct winding)
            // Triangle 2: vertices 2,3,4
            // etc.
           uint64_t base_index = tri_id;
            if (base_index + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
            
            if (src->index_array_buffer) {
                vertex_indices[0] = read_index(src->index_array_buffer, src->index_type, base_index);
                vertex_indices[1] = read_index(src->index_array_buffer, src->index_type, base_index + 1);
                vertex_indices[2] = read_index(src->index_array_buffer, src->index_type, base_index + 2);
            } else {
                vertex_indices[0] = base_index;
                vertex_indices[1] = base_index + 1;
                vertex_indices[2] = base_index + 2;
            }
            
            // Swap vertices 1 and 2 for odd triangles to maintain consistent winding
            if (tri_id & 1) {
                uint64_t temp = vertex_indices[1];
                vertex_indices[1] = vertex_indices[2];
                vertex_indices[2] = temp;
            }
            break;
        }
        
        case RINTINTIN_TRIANGLE_FAN: {
            // Triangle fan: all triangles share the first vertex
            // Triangle 0: vertices 0,1,2
            // Triangle 1: vertices 0,2,3
            // Triangle 2: vertices 0,3,4
            // etc.
            if (tri_id + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
            
            if (src->index_array_buffer) {
                vertex_indices[0] = read_index(src->index_array_buffer, src->index_type, 0);
                vertex_indices[1] = read_index(src->index_array_buffer, src->index_type, tri_id + 1);
                vertex_indices[2] = read_index(src->index_array_buffer, src->index_type, tri_id + 2);
            } else {
                vertex_indices[0] = 0;
                vertex_indices[1] = tri_id + 1;
                vertex_indices[2] = tri_id + 2;
            }
            break;
        }
        
        default:
            return RINTINTIN_ERROR_INVALID_TYPE;
    }
    
    // Validate vertex indices
    for (int i = 0; i < 3; i++) {
        if (vertex_indices[i] >= src->no_verts) {
            return RINTINTIN_ERROR_OUT_OF_BOUNDS;
        }
    }
    
    return RINTINTIN_SUCCESS;
}

static rintintin_error_code rintintin_get_visit_list(uint64_t * vertex_indices, const rintintin_mesh* src, uint64_t tri_id)
{
	int length = 0;

    switch (src->geometry_type) {
        case RINTINTIN_TRIANGLES: {
            // Simple case: each triangle uses 3 consecutive indices
           uint64_t base_index = tri_id * 3;
            if (base_index + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
            
            
            if (src->index_array_buffer == 0L) {
                // Non-indexed geometry
                vertex_indices[0] = base_index;
                vertex_indices[1] = base_index + 1;
                vertex_indices[2] = base_index + 2;
                break;
            }    
            
            uint64_t tmp[3];
            uint64_t cpy[3] = {vertex_indices[0], vertex_indices[1], vertex_indices[2]};
		
			// Indexed geometry
			tmp[0] = read_index(src->index_array_buffer, src->index_type, base_index);
			tmp[1] = read_index(src->index_array_buffer, src->index_type, base_index + 1);
			tmp[2] = read_index(src->index_array_buffer, src->index_type, base_index + 2);
			
// likely to read things with shared vertes so try to remove redundancy. 
			length = 0;
			
			for(int i = 0; i < 3; ++i)
			{
				int found = 0;
			
				for(int j = 0; j < 3; ++j)
				{
					if(cpy[j] == tmp[i])
					{
						found = 1;
						break;
					}
				}
				
				if(!found)
					vertex_indices[length++] = tmp[i];
			}

            break;
        }
        
        case RINTINTIN_TRIANGLE_STRIP: {
            // Triangle strip: triangles share vertices
            // Triangle 0: vertices 0,1,2
            // Triangle 1: vertices 1,2,3 (or 1,3,2 for correct winding)
            // Triangle 2: vertices 2,3,4
            // etc.
            uint64_t base_index = tri_id;
            if (base_index + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
             
            length = 1;
            
            if (src->index_array_buffer) {
                vertex_indices[0] = read_index(src->index_array_buffer, src->index_type, base_index + 2);
            } else {
                vertex_indices[0] = base_index + 2;
            }
            
            if(tri_id == 0)
            {
				length = 3;
				
				if (src->index_array_buffer) {
					vertex_indices[1] = read_index(src->index_array_buffer, src->index_type, base_index + 1);
					vertex_indices[2] = read_index(src->index_array_buffer, src->index_type, base_index);
				} else {
					vertex_indices[1] = base_index + 1;
					vertex_indices[2] = base_index;
				}				
            }
            
            break;
        }
        
        case RINTINTIN_TRIANGLE_FAN: {
            // Triangle fan: all triangles share the first vertex
            // Triangle 0: vertices 0,1,2
            // Triangle 1: vertices 0,2,3
            // Triangle 2: vertices 0,3,4
            // etc.
            if (tri_id + 2 >= src->no_indices) {
                return RINTINTIN_ERROR_OUT_OF_BOUNDS;
            }
               
            length = 1;
            
            if (src->index_array_buffer) {
                vertex_indices[0] = read_index(src->index_array_buffer, src->index_type, tri_id + 2);
            } else {
                vertex_indices[0] = tri_id + 2;
            }
            
            if(tri_id == 0)
            {
				length = 3;
				
				if (src->index_array_buffer) {
					vertex_indices[1] = read_index(src->index_array_buffer, src->index_type, tri_id + 1);
					vertex_indices[2] = read_index(src->index_array_buffer, src->index_type, 0);
				} else {
					vertex_indices[1] = tri_id + 1;
					vertex_indices[2] = 0;
				}				
            }
         
            break;
        }
        
        default:
            return RINTINTIN_ERROR_INVALID_TYPE;
    }
    
    // Validate vertex indices
    for (int i = 0; i < length; i++) {
        if (vertex_indices[i] >= src->no_verts) {
            return RINTINTIN_ERROR_OUT_OF_BOUNDS;
        }
    }
    
    return length;
}

static rintintin_error_code rintintin_get_tri(rintintin_tri* dst, const rintintin_mesh* src, uint64_t tri_id) 
{ 			
    uint64_t vertex_indices[3];
    
    // Read attributes for each vertex
    rintintin_error_code err = rintintin_get_indices(vertex_indices, src, tri_id);
    
    if(err != RINTINTIN_SUCCESS)
		return RINTINTIN_SUCCESS;
    
    for (int vert = 0; vert < 3; vert++) 
    {
		double  tmp_pos[4];
		double  tmp_weights[4];
		int32_t tmp_joint[4];
   
		err = src->position(tmp_pos, (uint32_t)vertex_indices[vert], src->position_user_data);
		if(err != RINTINTIN_SUCCESS) return err;
		
		err = src->joints(tmp_joint, (uint32_t)vertex_indices[vert], src->joints_user_data);
		if(err != RINTINTIN_SUCCESS) return err;
		
		err = src->weights(tmp_weights, (uint32_t)vertex_indices[vert], src->weights_user_data);
		if(err != RINTINTIN_SUCCESS) return err;
		
// normalize weights.		
		double sum = tmp_weights[0] + tmp_weights[1] + tmp_weights[2] + tmp_weights[3];
		
		if(sum > 1.0)
		{
			double invSum = 1.0 / sum;
		
			for (int i = 0; i < 4; i++) {
				tmp_weights[i] = (float)(tmp_weights[i] * invSum);
			}
		}
		
		for (int i = 0; i < 4; i++)
		{ 
			if(i < 3)
				dst->position[vert][i] = tmp_pos[i]; 
				
			dst->joint[vert][i] = (uint32_t)tmp_joint[i]; 
			dst->weight[vert][i] = (float)tmp_weights[i]; 
		}  
    }
    
    return RINTINTIN_SUCCESS;
}

// Helper function to count triangles in a mesh
uint32_t rintintin_get_triangle_count(const rintintin_mesh* src) {
    if (!src) return 0;
    
    switch (src->geometry_type) {
        case RINTINTIN_TRIANGLES:
            return (uint32_t)(src->no_indices / 3);
            
        case RINTINTIN_TRIANGLE_STRIP:
            return (src->no_indices >= 3) ? (uint32_t)(src->no_indices - 2) : 0;
            
        case RINTINTIN_TRIANGLE_FAN:
            return (src->no_indices >= 3) ? (uint32_t)(src->no_indices - 2) : 0;
            
        default:
            return 0;
    }
}

rintintin_error_code rintintin_process_mesh(rintintin_mesh const* src, uint32_t thread_id, uint32_t no_threads, rintintin_tri_callback * cb, void * userdata)
{
	if(no_threads == 0)
	{
		thread_id = 0;
		no_threads = 1;
	}
	
	if(thread_id >= no_threads)
		return RINTINTIN_ERROR_INVALID_INDEX;

	if(src == 0
	|| src->position == 0 || src->position_user_data == 0
	|| src->joints == 0 || src->joints_user_data == 0
	|| src->weights == 0 || src->weights_user_data == 0) 
		return RINTINTIN_ERROR_NULL_POINTER;
		
	uint64_t N = rintintin_get_triangle_count(src);
	uint64_t begin = (thread_id * N) / no_threads;
	uint64_t end = ((thread_id+1) * N) / no_threads;
	
	rintintin_tri tri = {0};
	for(uint64_t i = begin; i < end; ++i)
	{
		rintintin_error_code ec = rintintin_get_tri(&tri, src, i);
		if(RINTINTIN_SUCCESS != ec)
			return ec;
	
		int va = 0, vb = 1;
		if(src->surface_mode == RINTINTIN_SURFACE_FLIPPED)
		{
			va = 1;
			vb = 0;
		}
			
		int w[3] = {va, vb, 2};
		rintintin_vec3 p[3];
		p[0].x = tri.position[va][0];
		p[0].y = tri.position[va][1];
		p[0].z = tri.position[va][2];
		
		p[1].x = tri.position[vb][0];
		p[1].y = tri.position[vb][1];
		p[1].z = tri.position[vb][2];
		
		p[2].x = tri.position[2][0];
		p[2].y = tri.position[2][1];
		p[2].z = tri.position[2][2];
		
		double weights[3] = {0, 0, 0};
		
		for(uint32_t v = 0; v < 3; ++v)
		{
			for(uint32_t j = 0; j < 3; ++j)
			{
				if(tri.weight[v][j] == 0)
					continue;
					
				weights[w[v]] = tri.weight[v][j];
				uint32_t joint = tri.joint[v][j];
				
				for(uint32_t v1 = v+1; v1 < 3; ++v1)
				{							
					for(uint32_t j1 = 0; j1 < 3; ++j1)
					{
						if(tri.joint[v1][j1] == joint)
						{
							weights[w[v1]] = tri.weight[v1][j1];
							// don't process same joint again.
							tri.weight[v1][j1] = 0;
							break;
						}
					}
				}
			
				ec = cb(p, weights, (int32_t)joint, userdata);
				// clear for next time.
				weights[0] = 0;
				weights[1] = 0;
				weights[2] = 0;
				
				if(ec != RINTINTIN_SUCCESS)
					return ec;
			}
		}
	}
		
	return RINTINTIN_SUCCESS;
}

rintintin_error_code rintintin_visit_each_index(rintintin_mesh const* src, uint32_t thread_id, uint32_t no_threads, rintintin_vert_callback * cb, void * userdata)
{		
	if(no_threads == 0)
	{
		thread_id = 0;
		no_threads = 1;
	}
	
	if(thread_id >= no_threads)
		return RINTINTIN_ERROR_INVALID_INDEX;
		
	if(src == 0
	|| src->position == 0 || src->position_user_data == 0
	|| src->joints == 0 || src->joints_user_data == 0
	|| src->weights == 0 || src->weights_user_data == 0) 
		return RINTINTIN_ERROR_NULL_POINTER;
		
	uint64_t N = rintintin_get_triangle_count(src);
	uint64_t begin = (thread_id * N) / no_threads;
	uint64_t end = ((thread_id+1) * N) / no_threads;
	
	for(uint64_t i = begin; i < end; ++i)
	{
		uint64_t vertex_indices[3] = {(uint64_t)(-1ll), (uint64_t)(-1ll), (uint64_t)(-1ll)};
		int length = rintintin_get_visit_list(vertex_indices, src, i);
		
		struct rintintin_vertex vertex = {0};
		rintintin_error_code err;
		
		if(length < 0)
			return length;
		
		for (int vert = 0; vert < length; vert++) 
		{		
			double * tmp_weights = &vertex.weight.x;
		
			err = src->position(&vertex.position.x, (uint32_t)vertex_indices[vert], src->position_user_data);
			if(err != RINTINTIN_SUCCESS) return err;
			
			err = src->joints(&vertex.joint[0], (uint32_t)vertex_indices[vert], src->joints_user_data);
			if(err != RINTINTIN_SUCCESS) return err;
			
			err = src->weights(&vertex.weight.x, (uint32_t)vertex_indices[vert], src->weights_user_data);
			if(err != RINTINTIN_SUCCESS) return err;
			
		// normalize weights.		
			double sum = tmp_weights[0] + tmp_weights[1] + tmp_weights[2] + tmp_weights[3];
			
			if(sum > 1.0)
			{
				double invSum = 1.0 / sum;
			
				for (int i = 0; i < 4; i++) {
					tmp_weights[i] = (float)(tmp_weights[i] * invSum);
				}
			}
			
			err = cb(&vertex, userdata);
			
			if(err != RINTINTIN_SUCCESS)
				return err;
		}
	}
	
	return RINTINTIN_SUCCESS;
}
