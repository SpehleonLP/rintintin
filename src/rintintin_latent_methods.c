#define USING_NAMESPACE_RINTINTIN 1
#include "rintintin_latent_methods.h"
#include "rintintin_math.h"

#define pow2(x) (x)*(x)
#define pow3(x) (x)*(x)*(x)

static sdmat3 rintintin_compute_quadratic(rintintin_quadratic const* m, dvec3 const* p);

void rintintin_mat3x3_scale_eq(rintintin_mat3x3 * a, double b)
{
	for(int i = 0; i < 9; ++i)
		(&a->m.x.x)[i] *= b;
}

void rintintin_mat3x3_add_eq(rintintin_mat3x3 *RESTRICT a, rintintin_mat3x3 const*RESTRICT b)
{
	for(int i = 0; i < 9; ++i)
		(&a->m.x.x)[i] += (&b->m.x.x)[i];
}

static sdmat3 rintintin_compute_sym_jacobian(rintintin_sym_jacobian const* m, dvec3 const* p)
{
	sdmat3 r = {
		.xx=rintintin_vec3_dot(&m->xx, p),
		.yy=rintintin_vec3_dot(&m->yy, p),
		.zz=rintintin_vec3_dot(&m->zz, p),
		
		.xy=rintintin_vec3_dot(&m->xy, p),
		.xz=rintintin_vec3_dot(&m->xz, p),
		.yz=rintintin_vec3_dot(&m->yz, p),
	};
	
	return r;
}


static dvec3 rintintin_compute_coupling(rintintin_coupling const* m, dvec3 const* p)
{
	dvec3 r = {
		.x=rintintin_vec3_dot(&m->x, p),
		.y=rintintin_vec3_dot(&m->y, p),
		.z=rintintin_vec3_dot(&m->z, p),
	};
	
	return r;
}

static sdmat3 rintintin_compute_quadratic(rintintin_quadratic const* m, dvec3 const* p)
{ 
	sdmat3 r = {
		.xx = m->xx.xx * p->x*p->x + m->xx.xy * p->x*p->y + m->xx.xz * p->x*p->z,
		.yy = m->yy.yy * p->y*p->y + m->yy.xy * p->x*p->y + m->yy.yz * p->y*p->z,
		.zz = m->zz.zz * p->z*p->z + m->zz.xz * p->x*p->z + m->zz.yz * p->y*p->z,
		
		.xy = m->xy.xy * p->x*p->y + m->xy.xx * p->x*p->x + m->xy.xz * p->x*p->z + m->xy.yy * p->y*p->y + m->xy.yz * p->y*p->z,
		.xz = m->xz.xz * p->x*p->z + m->xz.xx * p->x*p->x + m->xz.xy * p->x*p->y + m->xz.yz * p->y*p->z + m->xz.zz * p->z*p->z,
		.yz = m->yz.yz * p->y*p->z + m->yz.xy * p->x*p->y + m->yz.xz * p->x*p->z + m->yz.yy * p->y*p->y + m->yz.zz * p->z*p->z,
	};

	return r;
}

static dvec3 rintintin_compute_cubic_o1(rintintin_o1 const* m, dvec3 const* p)
{
   double x_xx = m->x_xx * p->x * p->x;
   double x_xy = m->x_xy * p->x * p->y;
   double x_xz = m->x_xz * p->x * p->z;
   double y_xy = m->y_xy * p->x * p->y;
   double y_yy = m->y_yy * p->y * p->y;
   double y_yz = m->y_yz * p->y * p->z;
   double z_xz = m->z_xz * p->x * p->z;
   double z_yz = m->z_yz * p->y * p->z;
   double z_zz = m->z_zz * p->z * p->z;
	
	return (dvec3){
		.x=x_xx + x_xy + x_xz,
		.y=y_xy + y_yy + y_yz,
		.z=z_xz + z_yz + z_zz,
	};
}

/*
sdmat3 rintintin_evaluate(rintintin_tensor const* m, dvec3 const* o, dvec3 const* p)
{ 
	if(o->x == 0 && o->y == 0 && o->z == 0
	&& p->x == 0 && p->y == 0 && p->z == 0)
		return m->constant;
		
	sdmat3 r = m->constant;
	dvec3  l = m->linear.o;
	
	sdmat3 eax;	
	dvec3 tmp;
	
// compute p terms
// linear
	eax = rintintin_compute_sym_jacobian(&m->linear.j, p);
	r = rintintin_mat3_add(&eax, &r);

// quadratic
	eax = rintintin_compute_quadratic(&m->quadratic, p);
	r = rintintin_mat3_add(&eax, &r);

// cubic
	eax = rintintin_compute_sym_jacobian(&m->cubic.h, p);
	eax = rintintin_compute_square_sym(&eax, p);
	r = rintintin_mat3_add(&eax, &r);	

// compute o terms
	double mass = m->mass_o + dot(&m->cubic.o2, p); 
	tmp = rintintin_compute_cubic_o1(&m->cubic, p);
	l = vec3_add(&l, &tmp);
	
	tmp = rintintin_compute_coupling(&m->linear.coupling, p);
	l = vec3_add(&l, &tmp);

	eax = (sdmat3){
		.xx= mass * o->x*o->x + l.x * o->x,
		.yy= mass * o->y*o->y + l.y * o->y,
		.zz= mass * o->z*o->z + l.z * o->z,
		
		.xy= -mass * o->x*o->y - 0.5 * (l.x * o->y + l.y * o->x),
		.xz= -mass * o->x*o->z - 0.5 * (l.x * o->z + l.z * o->x),
		.yz= -mass * o->y*o->z - 0.5 * (l.y * o->z + l.z * o->y),
	};
	
	r = rintintin_mat3_add(&eax, &r);	
	
	return r;
}
*/


void rintintin_tensor_add_eq(rintintin_tensor *RESTRICT A, rintintin_tensor const*RESTRICT B)
{
	double * a = (double *)A;
	const double * b = (double const*)B;
	
	enum
	{
		SIZE = sizeof(*A) / sizeof(double)
	};
	
	for(unsigned int i = 0; i < SIZE; ++i)
		a[i] += b[i];
}

void rintintin_liquid_add_eq(rintintin_liquified *RESTRICT A, rintintin_liquified const*RESTRICT B)
{
	double * a = (double *)A;
	const double * b = (double const*)B;
	
	enum
	{
		SIZE = sizeof(*A) / sizeof(double)
	};
	
	for(unsigned int i = 0; i < SIZE; ++i)
		a[i] += b[i];
}

void rintintin_liquid_sub_eq(rintintin_liquified *RESTRICT A, rintintin_liquified const*RESTRICT B)
{
	double * a = (double *)A;
	const double * b = (double const*)B;
	
	enum
	{
		SIZE = sizeof(*A) / sizeof(double)
	};
	
	for(unsigned int i = 0; i < SIZE; ++i)
		a[i] += b[i];
}

void rintintin_tensor_sub_eq(rintintin_tensor *RESTRICT A, rintintin_tensor const*RESTRICT B)
{
	double * a = (double *)A;
	const double * b = (double const*)B;
	
	enum
	{
		SIZE = sizeof(*A) / sizeof(double)
	};
	
	for(unsigned int i = 0; i < SIZE; ++i)
		a[i] -= b[i];
}

rintintin_tensor rintintin_tensor_add(rintintin_tensor const* a, rintintin_tensor const* b)
{
	rintintin_tensor eax = *a;
	rintintin_tensor_add_eq(&eax, b);
	return eax;
}

rintintin_tensor rintintin_tensor_sub(rintintin_tensor const* a, rintintin_tensor const* b)
{
	rintintin_tensor eax = *a;
	rintintin_tensor_sub_eq(&eax, b);
	return eax;
}

rintintin_solidified rintintin_solidify(rintintin_tensor const* m, dvec3 const* p)
{
	if(p->x == 0 && p->y == 0 && p->z == 0)
		return (rintintin_solidified){
			.constant = m->constant,
			.linear_o = m->linear.o,
			.mass_o = m->mass_o,
		};	
		
	sdmat3 r = m->constant;
	dvec3  l = m->linear.o;
	sdmat3 eax;
	dvec3 tmp;
		
// compute p terms
// linear
	eax = rintintin_compute_sym_jacobian(&m->linear.j, p);
	r = smat3_add(&eax, &r);

// quadratic
	eax = rintintin_compute_quadratic(&m->quadratic, p);
	r = smat3_add(&eax, &r);

// cubic
	eax = rintintin_compute_sym_jacobian(&m->cubic.h, p);
	eax = rintintin_compute_square_sym(&eax, p);
	r = smat3_add(&eax, &r);
	
// compute o terms
	double mass = m->mass_o + dot(&m->cubic.o2, p); 
	tmp = rintintin_compute_cubic_o1(&m->cubic.o1, p);
	l = vec3_add(&l, &tmp);
	
	tmp = rintintin_compute_coupling(&m->linear.coupling, p);
	l = vec3_add(&l, &tmp);
	
	return (rintintin_solidified)
	{
		.constant = r,
		.linear_o = l,
		.mass_o = mass		
	};
}

void rintintin_solid_add_eq(rintintin_solidified * a, rintintin_solidified const* b)
{
	a->constant = smat3_add(&a->constant, &b->constant);
	a->linear_o = vec3_add(&a->linear_o, &b->linear_o);
	a->mass_o += b->mass_o;
}

void rintintin_tensor_sub_eq_solid(rintintin_tensor * a, rintintin_solidified const* b)
{
	a->constant = smat3_sub(&a->constant, &b->constant);
	a->linear.o = vec3_sub(&a->linear.o, &b->linear_o);
	a->mass_o -= b->mass_o;
}

rintintin_metrics rintintin_compute_metrics_from_solid(rintintin_solidified const* a)
{
	rintintin_metrics m;
	
	m.inertia = a->constant;	
	m.volume = a->mass_o;
	
	if(m.volume == 0.0)
		m.centroid = (rintintin_vec3){0, 0, 0};
	else
	{
		double invVolume = 1.0 / m.volume;
		m.centroid.x = 0.5 * a->linear_o.x * invVolume;
		m.centroid.y = 0.5 * a->linear_o.y * invVolume;
		m.centroid.z = 0.5 * a->linear_o.z * invVolume;
	}
			
	sdmat3 axis = parallel_axis(m.volume, &m.centroid);
	m.inertia = smat3_sub(&m.inertia, &axis);	

	return m;
}



rintintin_vec3 rintintin_evaluate_centroid(rintintin_tensor const * tensor, rintintin_vec3 const* p)
{
	rintintin_liquified liquid = {
		.cubic={.o1=tensor->cubic.o1, .o2=tensor->cubic.o2},
		.linear={.coupling=tensor->linear.coupling, .o=tensor->linear.o},
		.mass_o=tensor->mass_o,
	};
	
	return rintintin_get_center(&liquid, p);
}

rintintin_vec3 rintintin_get_center(rintintin_liquified const * m, rintintin_vec3 const* p)
{
	dvec3  l = m->linear.o;
	dvec3 tmp;
		
// compute p terms
// compute o terms
	double mass = m->mass_o + dot(&m->cubic.o2, p); 
	tmp = rintintin_compute_cubic_o1(&m->cubic.o1, p);
	l = vec3_add(&l, &tmp);
	
	tmp = rintintin_compute_coupling(&m->linear.coupling, p);
	l = vec3_add(&l, &tmp);
	
	double invMass = mass? 0.5 / mass : 0;
	
	return vec3_scale(&l, invMass);
}

rintintin_liquified rintintin_liquid_less(rintintin_liquified const * m, rintintin_vec3 const* p)
{
	dvec3  l = m->linear.o;
	dvec3 tmp;
		
// compute p terms
// compute o terms
	double mass = m->mass_o + dot(&m->cubic.o2, p); 
	tmp = rintintin_compute_cubic_o1(&m->cubic.o1, p);
	l = vec3_add(&l, &tmp);
	
	tmp = rintintin_compute_coupling(&m->linear.coupling, p);
	l = vec3_add(&l, &tmp);
	
	rintintin_liquified r = {0};
	r.mass_o = mass;
	r.linear.o = l;
	
	return r;
}
