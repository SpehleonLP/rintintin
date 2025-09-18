/**
 * @file rintintin.h
 * @brief A zero-dependency library (only stdint!) for computing metrics of skinned meshes.
 * 
 * Rintintin processes rigged/skinned 3D meshes to compute volume, centroid, and inertia tensor
 * properties per joint in a skeletal hierarchy. The library is designed to be pure (no global state),
 * thread-safe, and dependency-free - it performs no memory allocation internally.
 * 
 */

#ifndef RINTINTIN_H
#define RINTINTIN_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _MSC_VER
#define RINTINTIN_API
#else
// Define export/import macros
#ifdef BUILDING_RINTINTIN_DLL
    #define RINTINTIN_API __declspec(dllexport)
#else
    #define RINTINTIN_API __declspec(dllimport)
#endif
#endif

/**
 * @brief Error codes returned by rintintin functions.
 */
typedef enum {
    RINTINTIN_SUCCESS = 0,                ///< Operation completed successfully
    RINTINTIN_ERROR_INVALID_ARGUMENT = -1,    ///< Invalid argument given
    RINTINTIN_ERROR_NULL_POINTER = -2,    ///< A required pointer parameter was NULL
    RINTINTIN_ERROR_OUT_OF_BOUNDS = -3,   ///< Index or size parameter was out of valid range
    RINTINTIN_ERROR_INVALID_TYPE = -4,    ///< Invalid data type specified
    RINTINTIN_ERROR_INVALID_SIZE = -5,    ///< Invalid size parameter
    RINTINTIN_ERROR_INVALID_INDEX = -6,    ///< Invalid index parameter
    RINTINTIN_ERROR_THREAD_COUNT_MISMATCH = -7,
    RINTINTIN_ERROR_OUT_OF_ORDER = -8,
    RINTINTIN_ERROR_NOT_PREPARED = -9,
    RINTINTIN_ERROR_TOO_MANY_ROOTS = -10,  ///< parents doesn't form a proper DAG 
    RINTINTIN_ERROR_CONFIG_CHANGED = -11,  
} rintintin_error_code;

typedef enum {
	RINTINTIN_SURFACE_NORMAL = 0,			///< glCullFace(GL_BACK)
	RINTINTIN_SURFACE_FLIPPED = 1,			///< glCullFace(GL_FRONT)
	RINTINTIN_SURFACE_THIN_SHELL = 2		///< glCullFace(GL_DISABLE)
} rintintin_surface_mode;

/**
 * @brief Convert the error code to a human readable string
 */
RINTINTIN_API const char* rintintin_get_error_string(int ec);

/**
 * @brief OpenGL-compatible data type constants.
 * 
 * These match OpenGL type constants for compatibility with existing vertex data formats.
 */
typedef enum {
    RINTINTIN_TYPE_BYTE           = 0x1400,  ///< Signed 8-bit integer
    RINTINTIN_TYPE_UNSIGNED_BYTE  = 0x1401,  ///< Unsigned 8-bit integer
    RINTINTIN_TYPE_SHORT          = 0x1402,  ///< Signed 16-bit integer
    RINTINTIN_TYPE_UNSIGNED_SHORT = 0x1403,  ///< Unsigned 16-bit integer
    RINTINTIN_TYPE_INT            = 0x1404,  ///< Signed 32-bit integer
    RINTINTIN_TYPE_UNSIGNED_INT   = 0x1405,  ///< Unsigned 32-bit integer
    RINTINTIN_TYPE_FLOAT          = 0x1406,  ///< 32-bit floating point
    RINTINTIN_TYPE_DOUBLE         = 0x140A,  ///< 64-bit floating point
    RINTINTIN_TYPE_HALF_FLOAT     = 0x140B,  ///< 16-bit floating point
    RINTINTIN_TYPE_FIXED          = 0x140C,  ///< Fixed-point number
    RINTINTIN_TYPE_INT_2_10_10_10_REV = 0x8D9F,            ///< Packed signed 2.10.10.10 format
    RINTINTIN_TYPE_UNSIGNED_INT_2_10_10_10_REV = 0x8368,   ///< Packed unsigned 2.10.10.10 format
    RINTINTIN_TYPE_UNSIGNED_INT_10F_11F_11F_REV = 0x8C3B   ///< Packed 10.11.11 floating point format
} rintintin_type;

/**
 * @brief Supported geometry primitive types.
 */
typedef enum 
{
    RINTINTIN_TRIANGLES,        ///< Independent triangles
    RINTINTIN_TRIANGLE_STRIP,   ///< Triangle strip
    RINTINTIN_TRIANGLE_FAN,     ///< Triangle fan
} rintintin_geometry_type;

/**
 * @brief Function pointer type for reading floating-point vertex attributes.
 * 
 * @param dst Output array of 4 doubles. If fewer than 4 components are available,
 *            remaining components should be zero-filled.
 * @param index Vertex index to read from
 * @param layout User data pointer (typically points to mesh data structure)
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note The dst array will always have 4 elements allocated. Implementations should
 *       zero-fill any unused components.
 */
typedef rintintin_error_code rintintin_read_attrib_f(double* dst, uint32_t index, void const* layout);

/**
 * @brief Function pointer type for reading integer vertex attributes.
 * 
 * @param dst Output array of 4 long long integers. If fewer than 4 components are available,
 *            remaining components should be zero-filled.
 * @param index Vertex index to read from  
 * @param layout User data pointer (typically points to mesh data structure)
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note The dst array will always have 4 elements allocated. Implementations should
 *       zero-fill any unused components.
 */
typedef rintintin_error_code rintintin_read_attrib_i(int32_t * dst, uint32_t index, void const* layout);

/**
 * @brief 3D vector structure.
 * 
 * @note order matches GLM {x, y, z} so you can cast to GLM's dvec3
 */
typedef struct rintintin_vec3 {
    double x;  ///< X component
    double y;  ///< Y component 
    double z;  ///< Z component
} rintintin_vec3;

/**
 * @brief 4D vector structure.
 * 
 * @note order matches GLM {x, y, z, w} so you can cast to GLM's dvec4 or dquat
 */ 
typedef struct rintintin_vec4 {
	double x;  ///< X component
	double y;  ///< Y component 
	double z;  ///< Z component
	double w;  ///< W component
} rintintin_vec4;

/**
 * @brief Symmetric 3x3 matrix stored in compressed form.
 * 
 * Represents a symmetric matrix using only the unique elements:
 * | xx  xy  xz |
 * | xy  yy  yz |
 * | xz  yz  zz |
 */
typedef struct rintintin_symmetric_mat3 {
    double xx;  ///< [0,0] element
    double yy;  ///< [1,1] element
    double zz;  ///< [2,2] element
    double xy;  ///< [0,1] and [1,0] elements
    double xz;  ///< [0,2] and [2,0] elements
    double yz;  ///< [1,2] and [2,1] elements
} rintintin_symmetric_mat3;

/**
 * @brief 4D vector structure.
 * 
 * @note data layout matches GLM so you can cast to GLM's dmat3
 */ 
typedef union rintintin_mat3x3 {
	struct {
		rintintin_vec3 x;
		rintintin_vec3 y; 
		rintintin_vec3 z;
	} m;
	double dbl[3][3];
	rintintin_vec3 vec[3];
} rintintin_mat3x3;

/**
 * @brief Computed properties for a joint.
 */
typedef struct rintintin_metrics {	 
    double volume;                      ///< Volume of geometry influenced by this joint
    rintintin_vec3 centroid;            ///< Center of mass in mesh space
    rintintin_symmetric_mat3 inertia;   ///< Inertia tensor about the centroid 
} rintintin_metrics;
	

/**
 * @brief Mesh data structure for processing skinned geometry.
 * 
 * This structure defines how to access vertex data for a skinned mesh. All function
 * pointers must be provided except for density, which is optional.
 */
typedef struct rintintin_mesh {
    rintintin_read_attrib_f * position;     ///< Function to read vertex positions (vec3)
    rintintin_read_attrib_i * joints;       ///< Function to read joint indices (ivec4)
    rintintin_read_attrib_f * weights;      ///< Function to read joint weights (vec4)
    
    void const* joints_user_data;           ///< User data passed to joints function
    void const* position_user_data;         ///< User data passed to position function
    void const* weights_user_data;          ///< User data passed to weights function
    
    void const* index_array_buffer;         ///< Index buffer for indexed geometry
    rintintin_type index_type;              ///< Data type of indices
    rintintin_geometry_type geometry_type;  ///< Type of geometry primitives
    uint64_t no_verts;						 ///< Number of vertices
    uint64_t no_indices;					 ///< Number of indices
    
	rintintin_surface_mode surface_mode;    ///< Surface interpretation mode
    float thickness;						 ///< how thick should it be interpreted as if a thin shell (populate this even when not in thin shell mode, sometimes used as fallback).
} rintintin_mesh;

/**
 * @brief Skin data structure for processing skinned geometry.
 * 
 * This structure defines how the armature goes together and must be defined. 
 */
typedef struct rintintin_skin {
	const char ** bone_names;						///< used for debugging, rintintin does not read from this. 
	
    rintintin_vec3 const* joint_translation_mesh_space;   ///< populate by inversing the inverse bind pose matrices and getting the translation. 
    int32_t const* parents;                       ///< Parent joint indices, -1 for root joints (length: no_joints)
    uint32_t no_joints;                           ///< Number of joints
} rintintin_skin;

/**
 * @brief Processing data structure for rintintin.
 *
 * @note output and scratch space arrays must be defined, rintintin performs no memory allocation internally.
 */
typedef struct rintintin_process_command
{
	/**
	 * @brief array of meshes to be processed
	 * 
	 * @note you can use this to include multiple primitives as one mesh, 
	 *       or one mesh multiple times with different veretx weights.
	 * 
	 * @warning Do not submit alpha cards, fins, or shells as they can cause artifacts.
	 */	
	rintintin_mesh * meshes;
	rintintin_skin skin;

	const char * name;
    rintintin_metrics * results;                      ///< Output array for computed metrics (length: no_joints)
    
    void * scratch_space;                             ///< Temporary working memory
    uint64_t scratch_space_byte_length;				///< Size of scratch space in bytes
    uint32_t max_threads;								///< number of threads that will be used.
    uint32_t no_meshes;									///< count of meshes in the mesh array.
} rintintin_command;

/**
 * @brief Get required size for mesh processing buffer.
 * 
 * @param no_joints Number of joints in the skeleton
 * @param max_threads Number of threads that will process the mesh
 * @return Size in bytes required for the processing buffer
 */
RINTINTIN_API uint64_t rintintin_get_scratch_space_size(uint32_t no_joints, uint32_t max_threads);

/**
 * @brief Prepare scratch space for compuation.
 * 
 * This is O(no_joints)
 * 
 * @param cmd Process command with all required fields filled
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note Run only once.
 * @warning Do not submit stages out of order.
 */
RINTINTIN_API rintintin_error_code rintintin_begin(rintintin_command * cmd);

/**
 * @brief Preprocess mesh geometry to get volume functions (multi-threaded).
 * 
 * This function processes the array of meshes and converts to latent space
 * 
 * This is O(no_tris) (multithreaded) followed by O(no_joints)
 *
 * ExampleUsage (C++):
 *
 * ThreadPool.run(no_threads, [&](int thread_id) { rintintin_preprocess_mesh(cmd, thread_id, no_threads); });
 *
 * @param cmd Processing command, 
 * @param thread_id Current thread ID (0 to no_threads-1)
 * @param no_threads Total number of threads processing this mesh
 * @return RINTINTIN_SUCCESS on success, error code on failure, 1 when (this thread is) finished
 * 
 * @warning wait for all threads to finish before continuing. 
 */
RINTINTIN_API rintintin_error_code rintintin_read_mesh(rintintin_command * cmd, uint32_t thread_id, uint32_t no_threads);

/**
 * @brief Do a parallel reduction (multi-threaded).
 * 
 * This step is optional (does not change the result)
 * The user is responsible for synchronizing threads. 
 * 
 * If you don't understand: skip this. 
 *
 * @param cmd Processing command, 
 * @param thread_id Current thread ID (0 to no_threads-1)
 * @param no_threads Total number of threads processing this mesh
 * @return RINTINTIN_SUCCESS on success, error code on failure, true if there is more work to do. 
 * 
 * @warning wait for all threads to finish before continuing. 
 */
RINTINTIN_API rintintin_error_code rintintin_parallel_reduction(rintintin_command * cmd, uint32_t loop, uint32_t thread_id, uint32_t no_threads);


/**
 * @brief Compute final mass properties from processed tensors.
 * 
 * Calculates volume, centroid, and inertia tensor for each joint based on
 * the processed mesh data and joint hierarchy.
 * 
 * This is O(no_joints)
 * 
 * Actually generates the output. 
 * Requires mesh to be processed.
 *
 * @param cmd Process command with all required fields filled
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note Run only once, does not check the mesh struct.
 */
RINTINTIN_API rintintin_error_code rintintin_end(rintintin_command * cmd);


/**
 * @brief Attribute layout structure for generic vertex data reading.
 * 
 * Similar to glVertexAttribPointer parameters, this structure describes
 * how to interpret packed vertex attribute data.
 */
typedef struct rintintin_attrib {
    void* src;                        ///< Pointer to attribute data
    uint64_t byte_length;   ///< Total size of data buffer
    rintintin_type type;              ///< Data type of components
    uint8_t size;               ///< Number of components (1-4)
    uint8_t normalized;         ///< Whether to normalize fixed-point data to [0,1] or [-1,1]
    uint16_t stride;            ///< Byte offset between consecutive attributes
    uint32_t offset;              ///< Initial byte offset in buffer
} rintintin_attrib;

/**
 * @brief Generic function for reading floating-point attributes from packed data.
 * 
 * Reads vertex attribute data based on the layout specification and converts
 * to double precision floating point.
 * 
 * @param dst Output array of 4 doubles
 * @param index Vertex index to read
 * @param layout Pointer to rintintin_attrib structure describing the data layout
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note Can be used as a rintintin_read_attrib_f function pointer.
 */
RINTINTIN_API rintintin_error_code rintintin_read_attrib_generic_f(double* dst, uint32_t index, void const* layout);

/**
 * @brief Generic function for reading integer attributes from packed data.
 * 
 * Reads vertex attribute data based on the layout specification and converts
 * to long long integers.
 * 
 * @param dst Output array of 4 long long integers
 * @param index Vertex index to read
 * @param layout Pointer to rintintin_attrib structure describing the data layout
 * @return RINTINTIN_SUCCESS on success, error code on failure
 * 
 * @note Can be used as a rintintin_read_attrib_i function pointer.
 */
RINTINTIN_API rintintin_error_code rintintin_read_attrib_generic_i(int32_t * dst, uint32_t index, void const* layout);

/// visualization stuff. 


/**
 * @brief Best-fit primitive estimation for an inertia tensor.
 * 
 * Contains the primitive type and its oriented bounding parameters.
 */
typedef struct rintintin_inertia_estimation
{
    rintintin_vec3 scale;         ///< Dimensions: x=largest axis, z=smallest axis
    rintintin_vec4 rotation;      ///< Orientation quaternion
    rintintin_vec3 translation;   ///< Translation vector
} rintintin_inertia_estimation;

/**
 * @brief Estimate best-fit ellipsoid shapes for a collection of inertia tensors.
 * 
 * Uses ellipsoids because they provide the best looking fit.
 * To adapt this for capsules or other shapes:
 * 
 * - Run rintintin_compute_eigen on the tensor
 * - Convert to a rotation (rintintin_compute_rotation_quat), use the centroid as translation
 * - The eigenvalues correspond to the principal moments of inertia
 * 
 * The rotation matrix will be such that:
 *   - X+ = largest eigenvalue
 *   - Z+ = smallest eigenvalue
 * So eigenvalues {0,1,2} correspond to {λ_large, λ_mid, λ_small}.
 * 
 * For an ellipsoid with semi-axes a,b,c and total mass M:
 *   Ix = (1/5) M (b² + c²)
 *   Iy = (1/5) M (a² + c²)
 *   Iz = (1/5) M (a² + b²)
 * 
 * Back-solving gives:
 *   a² = (5 / (2M)) * (Iy + Iz - Ix)
 * 
 * Example:
 *   double large_axis_sq = 5.0 * (λ_mid + λ_large - λ_small) / (2.0 * mass);
 *   double large_axis = sqrt(fabs(large_axis_sq));
 * 
 * With the axis lengths known, compute predicted mass:
 *   M_pred = (4/3) * π * density * a * b * c
 * 
 * Compare M_pred to actual mass:
 *   - If relative error is small → ellipsoid is a consistent fit
 *   - If not → the axes don’t correspond to a realizable ellipsoid
 * 
 * Always check that all squared axes are positive before sqrt.
 * 
 * @param dst Output array of shape estimations
 * @param src Input metrics containing inertia tensor data
 * @param no_items Number of items to process
 */
RINTINTIN_API rintintin_error_code rintintin_estimate_shapes(rintintin_inertia_estimation * dst, rintintin_metrics const* src, uint64_t no_items);


/// next few things are mostly internal stuff, exposed because it might be useful.

/**
 * @brief Transform decomposition into translation, rotation, and scaling components.
 * 
 * Applies to vertices in Scale → Rotate → Translate order:
 * `final_vertex = translation + scaling * (rotation * vertex)`
 * 
 * @note rotation is stored as a quaternion (x,y,z,w)
 */
typedef struct rintintin_transform
{
    rintintin_vec3 translation;  ///< Translation vector
    rintintin_vec3 scaling;      ///< Non-uniform scale factors
    rintintin_vec4 rotation;     ///< Rotation quaternion
} rintintin_transform;

/**
 * @brief Convert 3x3 rotation matrix to quaternion.
 * @param matrix 3x3 rotation matrix in row-major order
 * @note won't work if it also has scale.
 * @return Quaternion representation (x,y,z,w)
 */
RINTINTIN_API rintintin_vec4 rintintin_quat_from_3x3(rintintin_mat3x3 const* src);

/**
 * @brief Decompose 4x4 transformation matrix into SRT components.
 * 
 * Optimized for standard TRS matrices. Will give garbage results if the matrix
 * contains perspective projection or skew transformations.
 * 
 * @param matrix 4x4 transformation matrix in row-major order
 * @return Transform decomposition with Scale-Rotate-Translate application order
 * @note Applies as: translation + scaling * (rotation * vertex)
 */
RINTINTIN_API rintintin_transform rintintin_transform_from_4x4(double (*matrix)[4][4]);

/**
 * @brief Eigenvalue decomposition result for symmetric 3x3 matrices.
 * 
 * Contains eigenvectors and eigenvalues sorted in descending order
 * (greatest to least eigenvalue).
 */
typedef struct rintintin_eigen
{
    rintintin_vec3 vectors[3];  ///< Eigenvectors [0]=minor, [1]=middle, [2]=major axis
    double values[3];           ///< Eigenvalues [0]=smallest, [1]=middle, [2]=largest
} rintintin_eigen;

/**
 * @brief Compute eigenvalue decomposition of a symmetric 3x3 matrix.
 * @param matrix Symmetric 3x3 matrix (typically inertia tensor)
 * @return Eigen decomposition sorted by eigenvalue magnitude (smallest first)
 */
RINTINTIN_API rintintin_eigen rintintin_compute_eigen(const rintintin_symmetric_mat3 *matrix);

/**
 * @brief Compute 3x3 rotation matrix from eigendecomposition.
 * 
 * Creates rotation where X-axis aligns with major axis, Z-axis with minor axis.
 * 
 * @param dst Output 3x3 rotation matrix in row-major order
 * @param eigen Eigendecomposition with sorted eigenvalues/vectors
 */
RINTINTIN_API void rintintin_compute_rotation_3x3(rintintin_mat3x3 * dst, rintintin_eigen const* eigen);

/**
 * @brief Compute quaternion rotation from eigendecomposition.
 * 
 * Creates rotation where X-axis aligns with major axis, Z-axis with minor axis.
 * 
 * @param eigen Eigendecomposition with sorted eigenvalues/vectors  
 * @return Quaternion representing the principal axis alignment
 */
RINTINTIN_API rintintin_vec4 rintintin_compute_rotation_quat(rintintin_eigen const* eigen);


/**
 * @brief Compute 3x3 inverse matrix
 * 
 * @return -1 singular, 0 success
 */
RINTINTIN_API int rintintin_compute_inverse_3x3(rintintin_mat3x3 * dst, rintintin_mat3x3 const* src);


#ifdef __cplusplus
}
#endif

#endif
