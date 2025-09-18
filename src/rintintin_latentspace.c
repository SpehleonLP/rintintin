#include "rintintin_latentspace.h"
#include "rintintin_latent_methods.h"
#include "rintintin_eigen.h"
#include "rintintin_sqrt.h"

typedef struct rintintin_input rintintin_tri;
typedef rintintin_vec3 dvec3;
typedef rintintin_mat3x3 dmat3;

#ifdef __GNUC__
    #pragma GCC push_options
    #pragma GCC optimize("O3,fast-math,finite-math-only")
    #pragma GCC optimize("no-signed-zeros")
    #pragma GCC optimize("associative-math")
    #pragma GCC optimize("no-math-errno")
    #pragma GCC optimize("unsafe-math-optimizations")
#elif defined(__clang__)
    #pragma clang optimize on
    #pragma clang fp contract(fast)
#elif defined(_MSC_VER)
    #pragma optimize("", on)
    #pragma float_control(precise, off)
    #pragma float_control(except, off)
    #pragma fp_contract(on)
#endif


#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))


/*
    result[0] = 1.0*ox**2*coeff[0] + 1.0*ox*coeff[1] + 1.0*coeff[2];
    result[1] = 1.0*oy**2*coeff[0] + 1.0*oy*coeff[3] + 1.0*coeff[4];
    result[2] = 1.0*oz**2*coeff[0] + 1.0*oz*coeff[5] + 1.0*coeff[6];
    result[3] = -1.0*ox*oy*coeff[0] - 0.5*ox*coeff[3] - 0.5*oy*coeff[1] + 1.0*coeff[7];
    result[4] = -1.0*ox*oz*coeff[0] - 0.5*ox*coeff[5] - 0.5*oz*coeff[1] + 1.0*coeff[8];
    result[5] = -1.0*oy*oz*coeff[0] - 0.5*oy*coeff[5] - 0.5*oz*coeff[3] + 1.0*coeff[9];
*/

void rintintin_thin_construct(double * coeff, rintintin_tri const* src, double thickness)
{    
    double p0x = src->p0.x;
    double p0y = src->p0.y;
    double p0z = src->p0.z;
    double p1x = src->p1.x;
    double p1y = src->p1.y;
    double p1z = src->p1.z;
    double p2x = src->p2.x;
    double p2y = src->p2.y;
    double p2z = src->p2.z;
    double d0 = src->d0;
    double d1 = src->d1;
    double d2 = src->d2;

   // Common subexpression elimination
    const double x0 = 2*d0;
    const double x1 = d2 + x0;
    const double x2 = 2*p0x;
    const double x3 = p1x*x2;
    const double x4 = p0y*x3;
    const double x5 = p2x*x2;
    const double x6 = p0y*x5;
    const double x7 = p0z*x3;
    const double x8 = p0z*x5;
    const double x9 = p1y*p2y;
    const double x10 = p1z*p2z;
    const double x11 = p1y*p1z;
    const double x12 = p0y*p0z;
    const double x13 = 2*x12;
    const double x14 = 2*p2z;
    const double x15 = 2*p2y;
    const double x16 = p2y*p2z;
    const double x17 = 2*p0y;
    const double x18 = p1y*x17;
    const double x19 = p1x*p2x;
    const double x20 = p2y*x17;
    const double x21 = 2*p0z;
    const double x22 = p1z*x21;
    const double x23 = p2z*x21;
    const double x24 = 2*x19;
    const double x25 = pow2(p0x);
    const double x26 = pow2(p1y);
    const double x27 = pow2(p1z);
    const double x28 = pow2(p2y);
    const double x29 = pow2(p2z);
    const double x30 = pow2(p0y);
    const double x31 = pow2(p1x);
    const double x32 = pow2(p2x);
    const double x33 = pow2(p0z);
    const double x34 = 2*x25;
    const double x35 = 2*x30;
    const double x36 = 2*x33;
    const double x37 = thickness*sqrt(p1y*x12*x14 + p1y*x16*x21 - p1y*x4 + p1y*x6 + p1z*x12*x15 + p1z*x16*x17 - p1z*x7 + p1z*x8 + p2y*x11*x21 + p2y*x4 - p2y*x6 + p2z*x11*x17 + p2z*x7 - p2z*x8 - x10*x24 + x10*x3 - x10*x34 - x10*x35 + x10*x5 - x11*x13 - 2*x11*x16 - x13*x16 + x18*x19 - x18*x29 - x18*x32 + x19*x20 + x19*x22 + x19*x23 - x19*x35 - x19*x36 - x20*x27 - x20*x31 - x22*x28 - x22*x32 - x23*x26 - x23*x31 - x24*x9 + x25*x26 + x25*x27 + x25*x28 + x25*x29 + x26*x29 + x26*x32 + x26*x33 - x26*x5 + x27*x28 + x27*x30 + x27*x32 - x27*x5 - x28*x3 + x28*x31 + x28*x33 - x29*x3 + x29*x30 + x29*x31 + x3*x9 + x30*x31 + x30*x32 + x31*x33 + x32*x33 - x34*x9 - x36*x9 + x5*x9);
    const double x38 = 2*d1;
    const double x39 = x1 + x38;
    const double x40 = 2*d2;
    const double x41 = d1 + x40;
    const double x42 = x0 + x41;
    const double x43 = 6*d0 + x38 + x40;
    const double x44 = (1.0/60.0)*x37;
    const double x45 = d0 + d1 + d2;
    const double x46 = p1x*x45;
    const double x47 = 3*d1;
    const double x48 = x1 + x47;
    const double x49 = 3*d2;
    const double x50 = d1 + x0 + x49;
    const double x51 = 12*d0 + x47 + x49;
    const double x52 = 3*d0;
    const double x53 = d2 + x38 + x52;
    const double x54 = p1x*x53;
    const double x55 = x41 + x52;
    const double x56 = p2x*x55;
    const double x57 = (1.0/360.0)*x37;
    const double x58 = p1y*x45;
    const double x59 = p1y*x53;
    const double x60 = p2y*x55;
    const double x61 = p1z*x45;
    const double x62 = p1z*x53;
    const double x63 = p2z*x55;
    const double x64 = p2x*x45;
    const double x65 = p1x*x48;
    const double x66 = p2x*x50;
    const double x67 = p0x*x51;
    const double x68 = x59 + x60;
    const double x69 = x54 + x56;
    const double x70 = x62 + x63;
    
    coeff[0] += (1.0/24.0)*x37*(d1 + x1);
    coeff[1] += x44*(p0x*x43 + p1x*x39 + p2x*x42);
    coeff[2] += x57*(p0x*(2*x54 + 2*x56) + 2*p2x*x46 + x25*x51 + x31*x48 + x32*x50);
    coeff[3] += x44*(p0y*x43 + p1y*x39 + p2y*x42);
    coeff[4] += x57*(p0y*(2*x59 + 2*x60) + x15*x58 + x26*x48 + x28*x50 + x30*x51);
    coeff[5] += x44*(p0z*x43 + p1z*x39 + p2z*x42);
    coeff[6] += x57*(p0z*(2*x62 + 2*x63) + x14*x61 + x27*x48 + x29*x50 + x33*x51);
    coeff[7] += -x57*(p0x*x68 + p0y*x67 + p0y*x69 + p1y*x64 + p1y*x65 + p2y*x46 + p2y*x66);
    coeff[8] += -x57*(p0x*x70 + p0z*x67 + p0z*x69 + p1z*x64 + p1z*x65 + p2z*x46 + p2z*x66);
    coeff[9] += -x57*(p0y*x70 + p0z*x68 + p2y*x61 + p2z*x58 + x11*x48 + x12*x51 + x16*x50);
}

void rintintin_thick_construct(double * coeff, rintintin_tri const* src)
{    
    double p0x = src->p0.x;
    double p0y = src->p0.y;
    double p0z = src->p0.z;
    double p1x = src->p1.x;
    double p1y = src->p1.y;
    double p1z = src->p1.z;
    double p2x = src->p2.x;
    double p2y = src->p2.y;
    double p2z = src->p2.z;
    double d0 = src->d0;
    double d1 = src->d1;
    double d2 = src->d2;

   // Common subexpression elimination
    const double x0 = (1.0/18.0)*d0 + (1.0/18.0)*d1 + (1.0/18.0)*d2;
    const double x1 = -x0;
    const double x2 = p1y - p2y;
    const double x3 = p0x*x2;
    const double x4 = -p1x + p2x;
    const double x5 = p0y*x4;
    const double x6 = p1x*p2y;
    const double x7 = p1y*p2x;
    const double x8 = x6 - x7;
    const double x9 = x3 + x5 + x8;
    const double x10 = p1z - p2z;
    const double x11 = p0x*x10;
    const double x12 = p0z*x4;
    const double x13 = p1x*p2z;
    const double x14 = p1z*p2x;
    const double x15 = x13 - x14;
    const double x16 = x11 + x12 + x15;
    const double x17 = p0y*x10;
    const double x18 = -p0z*x2;
    const double x19 = p1y*p2z;
    const double x20 = p1z*p2y;
    const double x21 = x19 - x20;
    const double x22 = x17 + x18 + x21;
    const double x23 = p0x*x21;
    const double x24 = -p0y*x15;
    const double x25 = p0z*x8;
    const double x26 = (1.0/36.0)*d0 + (1.0/36.0)*d1 + (1.0/36.0)*d2;
    const double x27 = -x26;
    const double x28 = (1.0/180.0)*x3;
    const double x29 = (1.0/180.0)*x5;
    const double x30 = (1.0/180.0)*p1y*p2x - x28 - x29 - 1.0/180.0*x6;
    const double x31 = 4*d1;
    const double x32 = 4*d2;
    const double x33 = 7*d0;
    const double x34 = x32 + x33;
    const double x35 = x31 + x34;
    const double x36 = p0x*x35;
    const double x37 = 7*d1;
    const double x38 = 4*d0;
    const double x39 = x32 + x38;
    const double x40 = x37 + x39;
    const double x41 = p1x*x40;
    const double x42 = 7*d2;
    const double x43 = x31 + x42;
    const double x44 = x38 + x43;
    const double x45 = p2x*x44;
    const double x46 = x41 + x45;
    const double x47 = x36 + x46;
    const double x48 = (1.0/180.0)*x11;
    const double x49 = (1.0/180.0)*x12 + (1.0/180.0)*x13 - 1.0/180.0*x14 + x48;
    const double x50 = p0y*x35;
    const double x51 = x48*x50;
    const double x52 = (1.0/180.0)*x18;
    const double x53 = 2*d0 - d1 - d2;
    const double x54 = (1.0/180.0)*x23;
    const double x55 = 2*d2;
    const double x56 = d0 + d1 - x55;
    const double x57 = x56*x7;
    const double x58 = 2*d1;
    const double x59 = d0 + d2 - x58;
    const double x60 = x59*x6;
    const double x61 = x57 - x60;
    const double x62 = -x40;
    const double x63 = p1x*p1y;
    const double x64 = p2y*x45 + x62*x63;
    const double x65 = (1.0/180.0)*p0z;
    const double x66 = p2z*x45;
    const double x67 = x13*x59;
    const double x68 = x14*x56;
    const double x69 = x67 - x68;
    const double x70 = (1.0/180.0)*p0y;
    const double x71 = (1.0/180.0)*x25;
    const double x72 = (1.0/180.0)*x24 + x54 + x71;
    const double x73 = (1.0/180.0)*d0 + (1.0/180.0)*d1 + (1.0/180.0)*d2;
    const double x74 = -x73;
    const double x75 = (1.0/1080.0)*x3;
    const double x76 = (1.0/1080.0)*x5;
    const double x77 = (1.0/1080.0)*p1y*p2x - 1.0/1080.0*x6 - x75 - x76;
    const double x78 = 8*d0;
    const double x79 = 5*d1;
    const double x80 = 5*d2;
    const double x81 = x78 + x79 + x80;
    const double x82 = p0x*x81;
    const double x83 = 5*d0;
    const double x84 = 8*d1;
    const double x85 = x80 + x83 + x84;
    const double x86 = p1x*x85;
    const double x87 = 8*d2;
    const double x88 = x79 + x83 + x87;
    const double x89 = p2x*x88;
    const double x90 = x86 + x89;
    const double x91 = x82 + x90;
    const double x92 = (1.0/1080.0)*x11;
    const double x93 = (1.0/1080.0)*x12 + (1.0/1080.0)*x13 - 1.0/1080.0*x14 + x92;
    const double x94 = p0y*x81;
    const double x95 = x92*x94;
    const double x96 = (1.0/1080.0)*x18;
    const double x97 = (1.0/1080.0)*x23;
    const double x98 = x53*x97;
    const double x99 = (1.0/1080.0)*x21;
    const double x100 = -x85;
    const double x101 = p2y*x89 + x100*x63;
    const double x102 = x101 + x61;
    const double x103 = (1.0/1080.0)*p0z;
    const double x104 = p2z*x89;
    const double x105 = p1z*x86 - x104 + x69;
    const double x106 = (1.0/1080.0)*p0y;
    const double x107 = x34 + x37;
    const double x108 = p1x*x107;
    const double x109 = x33 + x43;
    const double x110 = p2x*x109;
    const double x111 = x108 + x110;
    const double x112 = pow2(p0x);
    const double x113 = 10*d0 + x31 + x32;
    const double x114 = x37 + x38 + x42;
    const double x115 = pow2(p1x);
    const double x116 = 10*d1 + x39;
    const double x117 = x115*x116;
    const double x118 = pow2(p2x);
    const double x119 = 2*d0 + x58 + x80;
    const double x120 = 2*x119;
    const double x121 = p1x*p2x*x114 + x117 + x118*x120;
    const double x122 = p0x*x111 + x112*x113 + x121;
    const double x123 = 2*p2z;
    const double x124 = 2*p1z - x123;
    const double x125 = x55 + x58 + x83;
    const double x126 = x112*x125;
    const double x127 = 2*p1y - 2*p2y;
    const double x128 = p1x + p2x;
    const double x129 = -x107;
    const double x130 = p2y*x110;
    const double x131 = x129*x63 + x130;
    const double x132 = x131 + x61;
    const double x133 = p0x*x103;
    const double x134 = p2z*x110;
    const double x135 = p1z*x108 - x134;
    const double x136 = x135 + x69;
    const double x137 = p0x*x106;
    const double x138 = x115*x59;
    const double x139 = x118*x56;
    const double x140 = x118*x119;
    const double x141 = -x116;
    const double x142 = p1y*x115;
    const double x143 = 2*p2y;
    const double x144 = (1.0/1080.0)*x25;
    const double x145 = x144 + (1.0/1080.0)*x24 + x97;
    const double x146 = (1.0/180.0)*p1z*p2y - 1.0/180.0*x17 - 1.0/180.0*x19 - x52;
    const double x147 = p1y*x40;
    const double x148 = p2y*x44;
    const double x149 = x147 + x148;
    const double x150 = x149 + x50;
    const double x151 = p0z*x35;
    const double x152 = x151*x29;
    const double x153 = x15*x53;
    const double x154 = x56*x6;
    const double x155 = x59*x7;
    const double x156 = x154 - x155;
    const double x157 = x19*x59;
    const double x158 = x20*x56;
    const double x159 = x157 - x158;
    const double x160 = p1z*x147 - p2z*x148;
    const double x161 = (1.0/180.0)*p0x;
    const double x162 = (1.0/1080.0)*p1z*p2y - 1.0/1080.0*x17 - 1.0/1080.0*x19 - x96;
    const double x163 = p1y*x85;
    const double x164 = p2y*x88;
    const double x165 = x163 + x164;
    const double x166 = x165 + x94;
    const double x167 = p0z*x81;
    const double x168 = x167*x76;
    const double x169 = x106*x153;
    const double x170 = (1.0/1080.0)*x15;
    const double x171 = x101 + x156;
    const double x172 = p1z*x163 - p2z*x164;
    const double x173 = x159 + x172;
    const double x174 = (1.0/1080.0)*p0x;
    const double x175 = p1y*x107;
    const double x176 = p2y*x109;
    const double x177 = x175 + x176;
    const double x178 = pow2(p0y);
    const double x179 = pow2(p1y);
    const double x180 = x116*x179;
    const double x181 = pow2(p2y);
    const double x182 = p1y*p2y*x114 + x120*x181 + x180;
    const double x183 = p0y*x177 + x113*x178 + x182;
    const double x184 = x125*x178;
    const double x185 = -2*p1x + 2*p2x;
    const double x186 = p1y + p2y;
    const double x187 = x131 + x156;
    const double x188 = p0z*x106;
    const double x189 = p2z*x176;
    const double x190 = p1z*x175 - x189;
    const double x191 = x159 + x190;
    const double x192 = x179*x59;
    const double x193 = x181*x56;
    const double x194 = x119*x181;
    const double x195 = p1x*x141;
    const double x196 = 2*p2x;
    const double x197 = p1z*x40 + p2z*x44;
    const double x198 = x151 + x197;
    const double x199 = p1x*p1z;
    const double x200 = x13*x56;
    const double x201 = x14*x59;
    const double x202 = x200 - x201;
    const double x203 = x20*x59;
    const double x204 = x19*x56;
    const double x205 = x203 - x204;
    const double x206 = p1z*x85 + p2z*x88;
    const double x207 = x167 + x206;
    const double x208 = x144*x53;
    const double x209 = (1.0/1080.0)*x8;
    const double x210 = x100*x199 + x104 + x202;
    const double x211 = x172 + x205;
    const double x212 = p1z*x107 + p2z*x109;
    const double x213 = pow2(p0z);
    const double x214 = pow2(p1z);
    const double x215 = x116*x214;
    const double x216 = pow2(p2z);
    const double x217 = p1z*p2z*x114 + x120*x216 + x215;
    const double x218 = p0z*x212 + x113*x213 + x217;
    const double x219 = x125*x213;
    const double x220 = p1z + p2z;
    const double x221 = x129*x199 + x134;
    const double x222 = x202 + x221;
    const double x223 = x190 + x205;
    const double x224 = x214*x59;
    const double x225 = x216*x56;
    const double x226 = x119*x216;
    const double x227 = (1.0/2160.0)*x3;
    const double x228 = (1.0/2160.0)*x5;
    const double x229 = x167*x228;
    const double x230 = (1.0/2160.0)*x25;
    const double x231 = x230*x53;
    const double x232 = (1.0/2160.0)*x8;
    const double x233 = (1.0/2160.0)*p0y;
    const double x234 = (1.0/2160.0)*p0x;
    const double x235 = p0x*x177;
    const double x236 = 20*d0 + x84 + x87;
    const double x237 = p0x*x236;
    const double x238 = 20*d1 + x78 + x87;
    const double x239 = x238*x63;
    const double x240 = 4*p2y;
    const double x241 = p2x*x119;
    const double x242 = x114*x6 + x114*x7 + x239 + x240*x241;
    const double x243 = p0y*x111 + p0y*x237 + x235 + x242;
    const double x244 = 4*p1y - 4*p2y;
    const double x245 = -x234*x244;
    const double x246 = p0y*p0z;
    const double x247 = x125*x246;
    const double x248 = 4*p2z;
    const double x249 = 4*p1z - x248;
    const double x250 = x234*x249;
    const double x251 = (1.0/2160.0)*x23;
    const double x252 = x251*x53;
    const double x253 = (1.0/2160.0)*x18;
    const double x254 = p0z*x233;
    const double x255 = (1.0/2160.0)*x178;
    const double x256 = -x83;
    const double x257 = -d2 - x256 - x31;
    const double x258 = -d1 - x256 - x32;
    const double x259 = p0y*x234;
    const double x260 = (1.0/2160.0)*x242;
    const double x261 = p2x*x56;
    const double x262 = p1x*x59;
    const double x263 = -x80;
    const double x264 = d0 + x263 + x31;
    const double x265 = p2y*x264;
    const double x266 = -x79;
    const double x267 = d0 + x266 + x32;
    const double x268 = p1y*x267;
    const double x269 = -x238;
    const double x270 = 4*p2x;
    const double x271 = (1.0/2160.0)*p0z;
    const double x272 = p2z*x60;
    const double x273 = p1z*x57;
    const double x274 = 3*x6;
    const double x275 = d0 - d2;
    const double x276 = p1z*x275;
    const double x277 = x274*x276;
    const double x278 = 3*x7;
    const double x279 = d0 - d1;
    const double x280 = p2z*x279;
    const double x281 = x278*x280;
    const double x282 = d2 + x266 + x38;
    const double x283 = p1y*x13;
    const double x284 = d1 + x263 + x38;
    const double x285 = p2y*x14;
    const double x286 = x241*x248;
    const double x287 = p2y*x286;
    const double x288 = p1z*x239 - x287;
    const double x289 = x282*x283 - x284*x285 + x288;
    const double x290 = p2z*x155;
    const double x291 = p1z*x154;
    const double x292 = x274*x280 - x276*x278;
    const double x293 = x115*x56;
    const double x294 = x118*x59;
    const double x295 = p2x*x264;
    const double x296 = p1x*x267;
    const double x297 = (1.0/2160.0)*x112;
    const double x298 = p0z*x234;
    const double x299 = x153*x233;
    const double x300 = p0z*x111;
    const double x301 = -4*p1x + 4*p2x;
    const double x302 = x234*x247;
    const double x303 = x228*x300 + x301*x302;
    const double x304 = -x230 - 1.0/2160.0*x24 - x251;
    const double x305 = (1.0/2160.0)*x11;
    const double x306 = x305*x94;
    const double x307 = (1.0/2160.0)*x15;
    const double x308 = x114*x13 + x114*x14 + x199*x238 + x286;
    const double x309 = p0x*x212 + p0z*x237 + x300 + x308;
    const double x310 = p1z*x282;
    const double x311 = x310*x6;
    const double x312 = p2y*x201;
    const double x313 = p1y*x200;
    const double x314 = p2z*x284;
    const double x315 = x314*x7;
    const double x316 = p1z*x269;
    const double x317 = p1z*x267;
    const double x318 = p2z*x264;
    const double x319 = x214*x238;
    const double x320 = x287 + x316*x63;
    const double x321 = 3*x275*x283 - 3*x279*x285 + x320;
    const double x322 = (1.0/2160.0)*x213;
    const double x323 = (1.0/2160.0)*x21;
    const double x324 = p1y*p1z;
    const double x325 = p0y*x212;
    const double x326 = x247*x250 + x305*x325;
    const double x327 = p0z*x177;
    const double x328 = p2y*x119*x248 + x114*x19 + x114*x20 + x238*x324;
    const double x329 = x236*x246 + x325 + x327 + x328;
    const double x330 = -x310*x7 + x314*x6;
    
    coeff[0] += x1*x9;
    coeff[1] += x0*x16;
    coeff[2] += x1*x22;
    coeff[3] += x0*(x23 + x24 + x25);
    coeff[4] += x27*x9;
    coeff[5] += x16*x26;
    coeff[6] += x22*x27;
    coeff[7] += x30*x47;
    coeff[8] += x47*x49;
    coeff[9] += -1.0/180.0*x21*x46 - x36*x52 - x51 - x53*x54 - x65*(x61 + x64) - x70*(p1z*x41 - x66 + x69);
    coeff[10] += x47*x72;
    coeff[11] += x74*x9;
    coeff[12] += x16*x73;
    coeff[13] += x22*x74;
    coeff[14] += x77*x91;
    coeff[15] += x91*x93;
    coeff[16] += -x102*x103 - x105*x106 - x82*x96 - x90*x99 - x95 - x98;
    coeff[17] += x122*x77;
    coeff[18] += x122*x93;
    coeff[19] += x103*x126*x127 - x103*(p1x*x57 + p1y*x139 - p2x*x60 - p2y*x138 + x140*x143 + x141*x142) - x106*x124*x126 - x106*(-p1x*x68 + p1z*x117 - p1z*x139 + p2x*x67 + p2z*x138 - x123*x140) - x112*x53*x99 - x121*x99 - x128*x98 - x132*x133 - x136*x137;
    coeff[20] += x122*x145;
    coeff[21] += x146*x150;
    coeff[22] += x150*x30;
    coeff[23] += (1.0/180.0)*x149*x15 + x152 + x153*x70 + x161*(x159 + x160) + x51 + x65*(x156 + x64);
    coeff[24] += x150*x72;
    coeff[25] += x162*x166;
    coeff[26] += x166*x77;
    coeff[27] += x103*x171 + x165*x170 + x168 + x169 + x173*x174 + x95;
    coeff[28] += x162*x183;
    coeff[29] += x183*x77;
    coeff[30] += x103*x184*x185 + x103*(p1x*x193 + p1y*x154 - p2x*x192 - p2y*x155 + x179*x195 + x194*x196) + x124*x174*x184 + x137*x191 + x169*x186 + x170*x178*x53 + x170*x182 + x174*(-p1y*x158 + p1z*x180 - p1z*x193 + p2y*x157 + p2z*x192 - x123*x194) + x187*x188;
    coeff[31] += x145*x183;
    coeff[32] += x198*x49;
    coeff[33] += x146*x198;
    coeff[34] += -x151*x28 - x152 - x161*(x160 + x205) - 1.0/180.0*x197*x8 - x53*x71 - x70*(x199*x62 + x202 + x66);
    coeff[35] += x198*x72;
    coeff[36] += x207*x93;
    coeff[37] += x162*x207;
    coeff[38] += -x106*x210 - x167*x75 - x168 - x174*x211 - x206*x209 - x208;
    coeff[39] += x218*x93;
    coeff[40] += x162*x218;
    coeff[41] += -x106*x185*x219 - x106*(p1x*x225 + p1z*x200 - p2x*x224 - p2z*x201 + x195*x214 + x196*x226) - x127*x174*x219 - x133*x223 - x174*(p1y*x215 - p1y*x225 - p1z*x204 + p2y*x224 + p2z*x203 - x143*x226) - x188*x222 - x208*x220 - x209*x213*x53 - x209*x217;
    coeff[42] += x145*x218;
    coeff[43] += -x167*x227 - x206*x232 - x210*x233 - x211*x234 - x229 - x231;
    coeff[44] += x243*(x227 + x228 + (1.0/2160.0)*x6 - 1.0/2160.0*x7);
    coeff[45] += x132*x254 + x136*x255 + x184*x250 + x186*x252 + x21*x260 + x233*(x272 - x273 - x277 + x281 + x289) + x235*x253 + x245*x247 + x259*(x19*x257 + x190 - x20*x258) + x271*(p1x*x179*x269 + x179*x261 - x181*x262 + x194*x270 + x265*x7 - x268*x6);
    coeff[46] += -x126*x233*x249 - x128*x299 - x15*x260 - x187*x298 - x191*x297 - x234*(x289 + x290 - x291 + x292) - x259*(x13*x257 + x135 - x14*x258) - x271*(-p1y*x294 + p2y*x293 + x140*x240 + x142*x269 + x295*x6 - x296*x7) - x303;
    coeff[47] += x243*x304;
    coeff[48] += x165*x307 + x171*x271 + x173*x234 + x229 + x299 + x306;
    coeff[49] += x309*((1.0/2160.0)*p1z*p2x - 1.0/2160.0*x12 - 1.0/2160.0*x13 - x305);
    coeff[50] += x126*x244*x271 + x128*x231 + x222*x259 + x223*x297 + x232*x308 + x233*(-p1z*x294 + p2z*x293 + x115*x316 + x13*x295 - x14*x296 + x140*x248) + x234*(x288 + x292 + x311 + x312 - x313 - x315) + x298*(p1y*x108 - x130 + x257*x6 - x258*x7) + x303;
    coeff[51] += x132*x322 + x136*x254 + x219*x245 + x220*x252 + x233*(p1x*x319 + x13*x317 - x14*x318 - x214*x261 + x216*x262 - x226*x270) + x271*(-x272 + x273 - x311 + x315 + x321) + x298*(x129*x324 + x189 + x19*x258 - x20*x257) + x308*x323 + x326;
    coeff[52] += x304*x309;
    coeff[53] += -x102*x271 - x105*x233 - x252 - x253*x82 - x306 - x323*x90;
    coeff[54] += x329*((1.0/2160.0)*x17 + (1.0/2160.0)*x19 - 1.0/2160.0*x20 + x253);
    coeff[55] += x184*x271*x301 + x186*x231 + x222*x255 + x223*x259 + x227*x327 + x232*x328 + x233*(x277 - x281 - x312 + x313 + x320 + x330) + x234*(p1z*x179*x238 + p1z*x181*x59 - p2z*x179*x56 - x19*x265 - x194*x248 + x20*x268) + x244*x302 + x254*(x131 - x257*x7 + x258*x6);
    coeff[56] += -x187*x322 - x191*x298 - x219*x233*x301 - x220*x299 - x234*(p1y*x216*x59 + p1y*x319 - p2y*x214*x56 + x19*x317 - x20*x318 - x226*x240) - x254*(x13*x258 - x14*x257 + x221) - x271*(-x290 + x291 + x321 + x330) - x307*x328 - x326;
    coeff[57] += x304*x329;
}

rintintin_vec3 rintintin_v_func_from_thick(double const* coeff)
{
	return (rintintin_vec3){ coeff[2] * (4.0 / 3.0), coeff[1] * (4.0 / 3.0), coeff[0] * (4.0 / 3.0)};
}


rintintin_vec4 rintintin_v_func(rintintin_tensor const* t)
{
	return (rintintin_vec4){t->cubic.o2.x, t->cubic.o2.y, t->cubic.o2.z, t->mass_o};
}

void rintintin_add_thin_shell(rintintin_tensor * dst, double const* coeff)
{    
    dst->constant.xx += coeff[2];
    dst->constant.yy += coeff[4];
    dst->constant.zz += coeff[6];
    
    dst->constant.xy += coeff[7];
    dst->constant.xz += coeff[8];
    dst->constant.yz += coeff[9];
    
	dst->linear.o.x += coeff[1];
	dst->linear.o.y += coeff[3];
	dst->linear.o.z += coeff[5];
	
    dst->mass_o += coeff[0];
}


rintintin_tensor rintintin_tensor_from_coeff(double const* coeff)
{
    rintintin_tensor o = {0};
        
   // Common subexpression elimination
    const double x0 = -coeff[12];
    const double x1 = -coeff[13];
    const double x2 = -coeff[11];
    const double x3 = -coeff[15]/2;
    const double x4 = -coeff[25]/2;
    const double x5 = -coeff[14]/2;
    const double x6 = -coeff[26]/2;
    const double x7 = -coeff[37]/2;
    const double x8 = -coeff[36]/2;
    o.constant.xx = coeff[20];
    o.constant.xy = coeff[47];
    o.constant.xz = coeff[52];
    o.constant.yy = coeff[31];
    o.constant.yz = coeff[57];
    o.constant.zz = coeff[42];
    o.cubic.h.xx.x = coeff[13];
    o.cubic.h.xx.y = coeff[12];
    o.cubic.h.xx.z = coeff[11];
    o.cubic.h.xy.x = x1;
    o.cubic.h.xy.y = x0;
    o.cubic.h.xy.z = x2;
    o.cubic.h.xz.x = x1;
    o.cubic.h.xz.y = x0;
    o.cubic.h.xz.z = x2;
    o.cubic.h.yy.x = coeff[13];
    o.cubic.h.yy.y = coeff[12];
    o.cubic.h.yy.z = coeff[11];
    o.cubic.h.yz.x = x1;
    o.cubic.h.yz.y = x0;
    o.cubic.h.yz.z = x2;
    o.cubic.h.zz.x = coeff[13];
    o.cubic.h.zz.y = coeff[12];
    o.cubic.h.zz.z = coeff[11];
    o.cubic.o1.x_xx = coeff[6];
    o.cubic.o1.x_xy = coeff[5];
    o.cubic.o1.x_xz = coeff[4];
    o.cubic.o1.y_xy = coeff[6];
    o.cubic.o1.y_yy = coeff[5];
    o.cubic.o1.y_yz = coeff[4];
    o.cubic.o1.z_xz = coeff[6];
    o.cubic.o1.z_yz = coeff[5];
    o.cubic.o1.z_zz = coeff[4];
    o.cubic.o2.x = coeff[2];
    o.cubic.o2.y = coeff[1];
    o.cubic.o2.z = coeff[0];
    o.linear.coupling.x.x = coeff[9];
    o.linear.coupling.x.y = coeff[8];
    o.linear.coupling.x.z = coeff[7];
    o.linear.coupling.y.x = coeff[21];
    o.linear.coupling.y.y = coeff[23];
    o.linear.coupling.y.z = coeff[22];
    o.linear.coupling.z.x = coeff[33];
    o.linear.coupling.z.y = coeff[32];
    o.linear.coupling.z.z = coeff[34];
    o.linear.j.xx.x = coeff[19];
    o.linear.j.xx.y = coeff[18];
    o.linear.j.xx.z = coeff[17];
    o.linear.j.xy.x = coeff[45];
    o.linear.j.xy.y = coeff[46];
    o.linear.j.xy.z = coeff[44];
    o.linear.j.xz.x = coeff[51];
    o.linear.j.xz.y = coeff[49];
    o.linear.j.xz.z = coeff[50];
    o.linear.j.yy.x = coeff[28];
    o.linear.j.yy.y = coeff[30];
    o.linear.j.yy.z = coeff[29];
    o.linear.j.yz.x = coeff[54];
    o.linear.j.yz.y = coeff[56];
    o.linear.j.yz.z = coeff[55];
    o.linear.j.zz.x = coeff[40];
    o.linear.j.zz.y = coeff[39];
    o.linear.j.zz.z = coeff[41];
    o.linear.o.x = coeff[10];
    o.linear.o.y = coeff[24];
    o.linear.o.z = coeff[35];
    o.mass_o = coeff[3];
    o.quadratic.xx.xx = coeff[16];
    o.quadratic.xx.xy = coeff[15];
    o.quadratic.xx.xz = coeff[14];
    o.quadratic.xy.xx = x4;
    o.quadratic.xy.xy = coeff[43];
    o.quadratic.xy.xz = x6;
    o.quadratic.xy.yy = x3;
    o.quadratic.xy.yz = x5;
    o.quadratic.xz.xx = x7;
    o.quadratic.xz.xy = x8;
    o.quadratic.xz.xz = coeff[48];
    o.quadratic.xz.yz = x3;
    o.quadratic.xz.zz = x5;
    o.quadratic.yy.xy = coeff[25];
    o.quadratic.yy.yy = coeff[27];
    o.quadratic.yy.yz = coeff[26];
    o.quadratic.yz.xy = x7;
    o.quadratic.yz.xz = x4;
    o.quadratic.yz.yy = x8;
    o.quadratic.yz.yz = coeff[53];
    o.quadratic.yz.zz = x6;
    o.quadratic.zz.xz = coeff[37];
    o.quadratic.zz.yz = coeff[36];
    o.quadratic.zz.zz = coeff[38];

    return o; 
}

rintintin_liquified rintintin_liquid_from_coeff(double const* coeff)
{
	rintintin_liquified o;
	
	  // Common subexpression elimination    
    o.cubic.o1.x_xx = coeff[6];
    o.cubic.o1.x_xy = coeff[5];
    o.cubic.o1.x_xz = coeff[4];
    o.cubic.o1.y_xy = coeff[6];
    o.cubic.o1.y_yy = coeff[5];
    o.cubic.o1.y_yz = coeff[4];
    o.cubic.o1.z_xz = coeff[6];
    o.cubic.o1.z_yz = coeff[5];
    o.cubic.o1.z_zz = coeff[4];
    o.cubic.o2.x = coeff[2];
    o.cubic.o2.y = coeff[1];
    o.cubic.o2.z = coeff[0];
    o.linear.coupling.x.x = coeff[9];
    o.linear.coupling.x.y = coeff[8];
    o.linear.coupling.x.z = coeff[7];
    o.linear.coupling.y.x = coeff[21];
    o.linear.coupling.y.y = coeff[23];
    o.linear.coupling.y.z = coeff[22];
    o.linear.coupling.z.x = coeff[33];
    o.linear.coupling.z.y = coeff[32];
    o.linear.coupling.z.z = coeff[34];
    o.linear.o.x = coeff[10];
    o.linear.o.y = coeff[24];
    o.linear.o.z = coeff[35];
    o.mass_o = coeff[3];
    
    return o;
}

//static inline double pow2(double a) { return a * a; }

int rintintin_evaluate_constraints(rintintin_tensor const * it, rintintin_vec3 const* p, rintintin_vec3 * out)
{
    rintintin_symmetric_mat3 constraint;
    
    // Apply SymPy's Common Subexpression Elimination
    double temp0 = p->x*p->y;
    double temp1 = p->x*p->z;
    double temp2 = pow2(p->x);
    double temp3 = it->cubic.o2.x*p->x + it->cubic.o2.y*p->y + it->cubic.o2.z*p->z + it->mass_o;
    double temp4 = 1.0/temp3;
    double temp5 = it->linear.coupling.x.x*p->x;
    double temp6 = it->linear.coupling.x.y*p->y;
    double temp7 = it->linear.coupling.x.z*p->z;
    double temp8 = it->cubic.o1.x_xy*temp0;
    double temp9 = it->cubic.o1.x_xz*temp1;
    double temp10 = it->cubic.o1.x_xx*temp2;
    double temp11 = p->y*p->z;
    double temp12 = pow2(p->y);
    double temp13 = it->linear.coupling.y.x*p->x;
    double temp14 = it->linear.coupling.y.y*p->y;
    double temp15 = it->linear.coupling.y.z*p->z;
    double temp16 = it->cubic.o1.y_xy*temp0;
    double temp17 = it->cubic.o1.y_yz*temp11;
    double temp18 = it->cubic.o1.y_yy*temp12;
    double temp19 = pow2(p->z);
    double temp20 = it->linear.coupling.z.x*p->x;
    double temp21 = it->linear.coupling.z.y*p->y;
    double temp22 = it->linear.coupling.z.z*p->z;
    double temp23 = it->cubic.o1.z_xz*temp1;
    double temp24 = it->cubic.o1.z_yz*temp11;
    double temp25 = it->cubic.o1.z_zz*temp19;
    double temp26 = (1.0/2.0)*it->linear.o.y + (1.0/2.0)*temp13 + (1.0/2.0)*temp14 + (1.0/2.0)*temp15 + (1.0/2.0)*temp16 + (1.0/2.0)*temp17 + (1.0/2.0)*temp18;
    double temp27 = pow2(temp4);
    double temp28 = -temp3;
    double temp29 = temp27*temp28*((1.0/2.0)*it->linear.o.x + (1.0/2.0)*temp10 + (1.0/2.0)*temp5 + (1.0/2.0)*temp6 + (1.0/2.0)*temp7 + (1.0/2.0)*temp8 + (1.0/2.0)*temp9);
    double temp30 = (1.0/2.0)*it->linear.o.z + (1.0/2.0)*temp20 + (1.0/2.0)*temp21 + (1.0/2.0)*temp22 + (1.0/2.0)*temp23 + (1.0/2.0)*temp24 + (1.0/2.0)*temp25;

    // Constraint assignments
    constraint.xx = it->constant.xx + it->linear.j.xx.x*p->x + it->linear.j.xx.y*p->y + it->linear.j.xx.z*p->z + it->quadratic.xx.xx*temp2 + it->quadratic.xx.xy*temp0 + it->quadratic.xx.xz*temp1 + temp2*(it->cubic.h.xx.x*p->x + it->cubic.h.xx.y*p->y + it->cubic.h.xx.z*p->z) - 0.25*temp4*pow2(it->linear.o.x + temp10 + temp5 + temp6 + temp7 + temp8 + temp9);
    constraint.yy = it->constant.yy + it->linear.j.yy.x*p->x + it->linear.j.yy.y*p->y + it->linear.j.yy.z*p->z + it->quadratic.yy.xy*temp0 + it->quadratic.yy.yy*temp12 + it->quadratic.yy.yz*temp11 + temp12*(it->cubic.h.yy.x*p->x + it->cubic.h.yy.y*p->y + it->cubic.h.yy.z*p->z) - 0.25*temp4*pow2(it->linear.o.y + temp13 + temp14 + temp15 + temp16 + temp17 + temp18);
    constraint.zz = it->constant.zz + it->linear.j.zz.x*p->x + it->linear.j.zz.y*p->y + it->linear.j.zz.z*p->z + it->quadratic.zz.xz*temp1 + it->quadratic.zz.yz*temp11 + it->quadratic.zz.zz*temp19 + temp19*(it->cubic.h.zz.x*p->x + it->cubic.h.zz.y*p->y + it->cubic.h.zz.z*p->z) - 0.25*temp4*pow2(it->linear.o.z + temp20 + temp21 + temp22 + temp23 + temp24 + temp25);
    constraint.xy = it->constant.xy + it->linear.j.xy.x*p->x + it->linear.j.xy.y*p->y + it->linear.j.xy.z*p->z + it->quadratic.xy.xx*temp2 + it->quadratic.xy.xy*temp0 + it->quadratic.xy.xz*temp1 + it->quadratic.xy.yy*temp12 + it->quadratic.xy.yz*temp11 + temp0*(it->cubic.h.xy.x*p->x + it->cubic.h.xy.y*p->y + it->cubic.h.xy.z*p->z) - temp26*temp29;
    constraint.xz = it->constant.xz + it->linear.j.xz.x*p->x + it->linear.j.xz.y*p->y + it->linear.j.xz.z*p->z + it->quadratic.xz.xx*temp2 + it->quadratic.xz.xy*temp0 + it->quadratic.xz.xz*temp1 + it->quadratic.xz.yz*temp11 + it->quadratic.xz.zz*temp19 + temp1*(it->cubic.h.xz.x*p->x + it->cubic.h.xz.y*p->y + it->cubic.h.xz.z*p->z) - temp29*temp30;
    constraint.yz = it->constant.yz + it->linear.j.yz.x*p->x + it->linear.j.yz.y*p->y + it->linear.j.yz.z*p->z + it->quadratic.yz.xy*temp0 + it->quadratic.yz.xz*temp1 + it->quadratic.yz.yy*temp12 + it->quadratic.yz.yz*temp11 + it->quadratic.yz.zz*temp19 + temp11*(it->cubic.h.yz.x*p->x + it->cubic.h.yz.y*p->y + it->cubic.h.yz.z*p->z) - temp26*temp27*temp28*temp30;

	if(out) { *out = (dvec3){constraint.xx, constraint.yy, constraint.zz}; };
	
	if(((constraint.xx >= 0 && constraint.yy >= 0 && constraint.zz >= 0)
	|| (constraint.xx <= 0 && constraint.yy <= 0 && constraint.zz >= 0)) == 0)
	{
		return 0;
	}
	
// flip around	
	double eax = constraint.xx;
	double eay = constraint.yy;
	constraint.xx = eay + constraint.zz;
	constraint.yy = eax + constraint.zz;
	constraint.zz = eax + eay;

    return rintintin_sym_mat3_classify(&constraint);
}

#ifdef __GNUC__
    #pragma GCC pop_options
#elif defined(__clang__)
    #pragma clang optimize off
#elif defined(_MSC_VER)
    #pragma optimize("", off)
    #pragma float_control(precise, on)
    #pragma float_control(except, on)
#endif
