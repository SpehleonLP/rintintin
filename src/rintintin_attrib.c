#include "../include/rintintin.h"

typedef uint64_t size_t;

// Helper function to get size of a type in bytes
static size_t get_type_size(rintintin_type type) {
    switch (type) {
        case RINTINTIN_TYPE_BYTE:
        case RINTINTIN_TYPE_UNSIGNED_BYTE:
            return 1;
        case RINTINTIN_TYPE_SHORT:
        case RINTINTIN_TYPE_UNSIGNED_SHORT:
        case RINTINTIN_TYPE_HALF_FLOAT:
            return 2;
        case RINTINTIN_TYPE_INT:
        case RINTINTIN_TYPE_UNSIGNED_INT:
        case RINTINTIN_TYPE_FLOAT:
        case RINTINTIN_TYPE_FIXED:
        case RINTINTIN_TYPE_INT_2_10_10_10_REV:
        case RINTINTIN_TYPE_UNSIGNED_INT_2_10_10_10_REV:
        case RINTINTIN_TYPE_UNSIGNED_INT_10F_11F_11F_REV:
            return 4;
        case RINTINTIN_TYPE_DOUBLE:
            return 8;
        default:
            return 0;
    }
}

// Helper to convert half float to float
static float half_to_float(uint16_t h) {
    uint32_t sign = (h & 0x8000) << 16;
    uint32_t exponent = (h & 0x7C00) >> 10;
    uint32_t mantissa = h & 0x03FF;
    
    if (exponent == 0) {
        if (mantissa == 0) return *(float*)&sign; // ±0
        // Subnormal
        exponent = 1;
        while (!(mantissa & 0x0400)) {
            mantissa <<= 1;
            exponent--;
        }
        mantissa &= 0x03FF;
    } else if (exponent == 31) {
        // Inf or NaN
        exponent = 255;
    }
    
    exponent = exponent + 112;
    mantissa <<= 13;
    
    uint32_t result = sign | (exponent << 23) | mantissa;
    return *(float*)&result;
}

// Read floating-point attributes
rintintin_error_code rintintin_read_attrib_generic_f(double* dst, uint32_t index, void const* user_data) {
	rintintin_attrib * layout = (rintintin_attrib*)user_data;
	void * src = layout->src;
	unsigned long long byte_length = layout->byte_length;
	
    if (!dst || !src || !layout) {
        return RINTINTIN_ERROR_NULL_POINTER;
    }
    
    if (layout->size < 1 || layout->size > 4) {
        return RINTINTIN_ERROR_INVALID_SIZE;
    }
    
    size_t type_size = get_type_size(layout->type);
    if (type_size == 0) {
        return RINTINTIN_ERROR_INVALID_TYPE;
    }
    
    size_t stride = layout->stride;
    if (stride == 0) {
        stride = type_size * layout->size;
    }
    
    size_t offset = layout->offset + (index * stride);
    size_t required_size = offset + (type_size * layout->size);
    
    if (required_size > byte_length) {
        return RINTINTIN_ERROR_OUT_OF_BOUNDS;
    }
    
    uint8_t* base = (uint8_t*)src + offset;
    
    // Initialize output to 0
    for (int i = 0; i < 4; i++) { dst[i] = 0.0; }
    
	switch (layout->type) {
        default:
            return RINTINTIN_ERROR_INVALID_TYPE;
		case RINTINTIN_TYPE_FLOAT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((float*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_DOUBLE: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((double*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_HALF_FLOAT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = half_to_float(((uint16_t*)(base))[i]);
            }
        } break;
		case RINTINTIN_TYPE_BYTE: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((int8_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 127.0) : (double)eax;
            }
        } break;
		case RINTINTIN_TYPE_SHORT: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((int16_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 32767.0) : (double)eax;
            }
        } break;
		case RINTINTIN_TYPE_INT: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((int32_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 2147483647.0) : (double)eax;
            }
        } break;
        case RINTINTIN_TYPE_UNSIGNED_BYTE: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((uint8_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 255.0) : (double)eax;
            }
        } break;
		case RINTINTIN_TYPE_UNSIGNED_SHORT: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((uint16_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 65535.0) : (double)eax;
            }
        } break;
		case RINTINTIN_TYPE_UNSIGNED_INT: 
		{
			for (int i = 0; i < layout->size; i++) {
				int64_t eax = ((uint32_t*)(base))[i];
                dst[i] = layout->normalized? ((double)eax / 4294967295.0) : (double)eax;
            }
        } break;
	 // Special packed formats (only when size matches expected)
		case RINTINTIN_TYPE_INT_2_10_10_10_REV:
		case RINTINTIN_TYPE_UNSIGNED_INT_2_10_10_10_REV: {
			float tmp[4];			
			
			uint32_t packed = *(uint32_t*)base;
			if (layout->type == RINTINTIN_TYPE_INT_2_10_10_10_REV) {
				// Signed
				int32_t x = (int32_t)(packed << 22) >> 22;  // Sign extend
				int32_t y = (int32_t)(packed << 12) >> 22;
				int32_t z = (int32_t)(packed << 2) >> 22;
				int32_t w = (int32_t)(packed >> 30) & 0x3;
				
				tmp[0] = (float)(layout->normalized ? (x / 511.0) : x);
				tmp[1] = (float)(layout->normalized ? (y / 511.0) : y);
				tmp[2] = (float)(layout->normalized ? (z / 511.0) : z);
				tmp[3] = (float)(layout->normalized ? (w / 1.0) : w);
			} else {
				// Unsigned
				uint32_t x = packed & 0x3FF;
				uint32_t y = (packed >> 10) & 0x3FF;
				uint32_t z = (packed >> 20) & 0x3FF;
				uint32_t w = (packed >> 30) & 0x3;
				
				tmp[0] = (float)(layout->normalized ? (x / 1023.0) : x);
				tmp[1] = (float)(layout->normalized ? (y / 1023.0) : y);
				tmp[2] = (float)(layout->normalized ? (z / 1023.0) : z);
				tmp[3] = (float)(layout->normalized ? (w / 3.0) : w);
			}
				
			for (int i = 0; i < layout->size; i++) {
                dst[i] = tmp[i];
            }
		} break;
    }
    
    return RINTINTIN_SUCCESS;
}

// Read integer attributes
rintintin_error_code rintintin_read_attrib_generic_i(int32_t * dst, uint32_t index, void const* user_data) {
	rintintin_attrib * layout = (rintintin_attrib*)user_data;
	void * src = layout->src;
	unsigned long long byte_length = layout->byte_length;

    if (!dst || !src || !layout) {
        return RINTINTIN_ERROR_NULL_POINTER;
    }
    
    if (layout->size < 1 || layout->size > 4) {
        return RINTINTIN_ERROR_INVALID_SIZE;
    }
        
    size_t type_size = get_type_size(layout->type);
    if (type_size == 0) {
        return RINTINTIN_ERROR_INVALID_TYPE;
    }
    
    size_t stride = layout->stride;
    if (stride == 0) {
        stride = type_size * layout->size;
    }
    
    size_t offset = layout->offset + (index * stride);
    size_t required_size = offset + (type_size * layout->size);
    
    if (required_size > byte_length) {
        return RINTINTIN_ERROR_OUT_OF_BOUNDS;
    }
    
    uint8_t* base = (uint8_t*)src + offset;
    
    // Initialize output
    for (int i = 0; i < 4; i++) { dst[i] = 0; }
    
    // Initialize output to 0
    for (int i = 0; i < 4; i++) { dst[i] = 0.0; }
    
	switch (layout->type) {
        default:
            return RINTINTIN_ERROR_INVALID_TYPE;
		case RINTINTIN_TYPE_FLOAT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = (int32_t)((float*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_DOUBLE: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = (int32_t)((double*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_HALF_FLOAT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = (int32_t)half_to_float(((uint16_t*)(base))[i]);
            }
        } break;
		case RINTINTIN_TYPE_BYTE: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((int8_t*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_SHORT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((int16_t*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_INT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((int32_t*)(base))[i];
            }
        } break;
        case RINTINTIN_TYPE_UNSIGNED_BYTE: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((uint8_t*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_UNSIGNED_SHORT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = ((uint16_t*)(base))[i];
            }
        } break;
		case RINTINTIN_TYPE_UNSIGNED_INT: 
		{
			for (int i = 0; i < layout->size; i++) {
                dst[i] = (int32_t)((uint32_t*)(base))[i];
            }
        } break;
	 // Special packed formats (only when size matches expected)
		case RINTINTIN_TYPE_INT_2_10_10_10_REV:
		case RINTINTIN_TYPE_UNSIGNED_INT_2_10_10_10_REV: {
			int32_t tmp[4];			
			
			uint32_t packed = *(uint32_t*)base;
			if (layout->type == RINTINTIN_TYPE_INT_2_10_10_10_REV) {
				// Signed
				tmp[0] = (int32_t)(packed << 22) >> 22;  // Sign extend
				tmp[1] = (int32_t)(packed << 12) >> 22;
				tmp[2] = (int32_t)(packed << 2) >> 22;
				tmp[3] = (int32_t)(packed >> 30) & 0x3;
			} else {
				// Unsigned
				tmp[0] = (int32_t)(packed & 0x3FF);
				tmp[1] = (int32_t)((packed >> 10) & 0x3FF);
				tmp[2] = (int32_t)((packed >> 20) & 0x3FF);
				tmp[3] = (int32_t)((packed >> 30) & 0x3);
			}
				
			for (int i = 0; i < layout->size; i++) {
                dst[i] = tmp[i];
            }
		} break;
    }
    
    
    return RINTINTIN_SUCCESS;
}

