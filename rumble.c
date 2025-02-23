#include <kos.h>
#include <stdlib.h>

//#include "map.h"
//#include "cheader.h"
#include "pal1.h"
#include "vrbg.h"

#include "beginner.h"

typedef float __attribute__((aligned(32))) Matrix[4][4];

// Structure to represent a 3D point
typedef struct {
    float x, y, z;
} Point3D;

// Compute cross product of two 3D vectors
Point3D crossProduct(Point3D a, Point3D b) {
    Point3D result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

// Compute dot product of two 3D vectors
float dotProduct(Point3D a, Point3D b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Subtract two points (vectors)
Point3D subtract(Point3D a, Point3D b) {
    Point3D result = {a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

// Project point P onto the plane defined by triangle (A, B, C)
Point3D projectPointOntoPlane(Point3D A, Point3D B, Point3D C, Point3D P) {
    Point3D AB = subtract(B, A);
    Point3D AC = subtract(C, A);
    Point3D normal = crossProduct(AB, AC);

    float D = -dotProduct(normal, A);
    float t = -(dotProduct(normal, P) + D) / dotProduct(normal, normal);

    Point3D projected = {
        P.x + t * normal.x,
        P.y + t * normal.y,
        P.z + t * normal.z
    };
    return projected;
}

// Check if a point P is inside triangle ABC using barycentric coordinates
bool isPointInsideTriangle(Point3D A, Point3D B, Point3D C, Point3D P) {
    Point3D v0 = subtract(C, A);
    Point3D v1 = subtract(B, A);
    Point3D v2 = subtract(P, A);

    float dot00 = dotProduct(v0, v0);
    float dot01 = dotProduct(v0, v1);
    float dot02 = dotProduct(v0, v2);
    float dot11 = dotProduct(v1, v1);
    float dot12 = dotProduct(v1, v2);

    float denom = dot00 * dot11 - dot01 * dot01;
    if (denom == 0.0f) return false;  // Degenerate triangle

    float invDenom = 1.0f / denom;
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}

// Main function to check if point is above and inside a 3D triangle
bool isPointAboveAndInsideTriangle(Point3D A, Point3D B, Point3D C, Point3D P) {
    // Project point P onto the plane of the triangle
    Point3D projected = projectPointOntoPlane(A, B, C, P);

    // Check if P is above the plane (compare original and projected y)
    if (P.y < projected.y) return false;

    // Check if the projected point lies inside the triangle
    return isPointInsideTriangle(A, B, C, projected);
}


Matrix R_ViewportMatrix;
Matrix R_ModelMatrix;
Matrix R_ProjectionMatrix;
Matrix R_RotX;
Matrix R_RotY;
Matrix R_Tran;

float camerax = 0.0f;
float cameray = 0.0f;
float cameraz = 0.0f;

float viewangle = 0.0f;
float viewpitch = 0.0f;

// must have
// `float t` and `float invt = 1.0f - t;`
// defined local to calling function
// this is just for cleaner code
#define lerp(a, b) (invt * (a) + t * (b))
// lerp two 32-bit colors
static uint32_t color_lerp(float ft, uint32_t c1, uint32_t c2) {
	uint8_t t = (ft * 255);
   	uint32_t maskRB = 0xFF00FF;  // Mask for Red & Blue channels
    uint32_t maskG  = 0x00FF00;  // Mask for Green channel
    uint32_t maskA  = 0xFF000000; // Mask for Alpha channel

    // Interpolate Red & Blue
    uint32_t rb = ((((c2 & maskRB) - (c1 & maskRB)) * t) >> 8) + (c1 & maskRB);
    
    // Interpolate Green
    uint32_t g  = ((((c2 & maskG) - (c1 & maskG)) * t) >> 8) + (c1 & maskG);

    // Interpolate Alpha
    uint32_t a  = ((((c2 & maskA) >> 24) - ((c1 & maskA) >> 24)) * t) >> 8;
    a = (a + (c1 >> 24)) << 24;  // Shift back into position

    return (a & maskA) | (rb & maskRB) | (g & maskG);
}

static void nearz_clip(float x0,float y0,float z0,float w0,float x1,float y1, float z1,float w1,
				float *xout, float *yout, float *zout, float *wout)
{
	const float d0 = w0 + z0;
	const float d1 = w1 + z1;

	float t = (fabsf(d0) * (1.0f / sqrtf((d1 - d0) * (d1 - d0)))) + 0.000001f;
	float invt = 1.0f - t;

	*wout = lerp(w0, w1);
	*xout = lerp(x0, x1);
	*yout = lerp(y0, y1);
	*zout = lerp(z0, z1);
}


static inline void perspdiv(float *x, float *y, float *z, float w)
{
	if (w == 0) w = 0.0001f;
	float invw = 1.0f / w;
	*x = (*x * invw);
	*y = (-*y *invw);

	if (w == 1.0f)
		*z = 1.0f / (1.0001f + (*z));
	else
		*z = invw;
}

// lifted from modern libultra
static inline void R_Ident(Matrix mf)
{
	int i, j;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			if (i == j)
				mf[i][j] = 1.0;
			else
				mf[i][j] = 0.0;
		}
	}
}
// the matrix setup in Doom 64 RE looks like it was "inlined"
// I broke it out into functions to better understand what it was doing
// and how viewproj was constructed

static inline void R_Translate(Matrix mf, float x, float y, float z)
{
	R_Ident(mf);
	mf[3][0] = x;
	mf[3][1] = y;
	mf[3][2] = z;
}
static inline void R_RotateX(Matrix mf, float in_sin, float in_cos)
{
	R_Ident(mf);
	mf[1][1] = in_cos;
	mf[1][2] = -in_sin;
	mf[2][1] = in_sin;
	mf[2][2] = in_cos;
}
static inline void R_RotateY(Matrix mf, float in_sin, float in_cos)
{
	R_Ident(mf);
	mf[0][0] = in_sin;
	mf[0][2] = -in_cos;
	mf[2][0] = in_cos;
	mf[2][2] = in_sin;
}

static inline void R_RotateZ(Matrix mf, float in_sin, float in_cos)
{
	R_Ident(mf);
	mf[0][0] = in_cos;
	mf[0][1] = -in_sin;
	mf[1][0] = in_sin;
	mf[1][1] = in_cos;
}


static inline void R_Frustum(Matrix mf, float l, float r, float b, float t, float n, float f, float scale)
{
	int i, j;
	R_Ident(mf);
	mf[0][0] = 2 * n / (r - l);
	mf[1][1] = 2 * n / (t - b);
	mf[2][0] = (r + l) / (r - l);
	mf[2][1] = (t + b) / (t - b);
	mf[2][2] = -(f + n) / (f - n);
	mf[2][3] = -1;
	mf[3][2] = -2 * f * n / (f - n);
	mf[3][3] = 0;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			mf[i][j] *= scale;
		}
	}
}

// I derived this from the glDC commit that moved screenspace transform into matrix
static inline void R_Viewport(Matrix mf, int x, int y, int width, int height) {
    mf[0][0] = -(float)width / 2.0f;
    mf[1][1] = (float)height / 2.0f;
    mf[2][2] = 1.0f;
    mf[3][3] = 1.0f;

    mf[3][0] = 640.0f - ((float)x + ((float)width / 2.0f));
    mf[3][1] = ((float)y + ((float)height / 2.0f));
}

#define get_color_argb1555(rrr, ggg, bbb, aaa)						\
	((uint16_t)(((aaa & 1) << 15) | (((rrr >> 3) & 0x1f) << 10) |	\
		    (((ggg >> 3) & 0x1f) << 5) | ((bbb >> 3) & 0x1f)))


KOS_INIT_FLAGS(INIT_DEFAULT);
pvr_init_params_t pvr_params = { { 0, 0, PVR_BINSIZE_16, 0, 0 },
				(2048576*2) / 2,
				1, // dma enabled
				0, // fsaa
				0, // 1 is autosort disabled
				2, // extra OPBs
				0, // Vertex buffer double-buffering enabled
 };

#define LINEWIDTH 1.5f
pvr_vertex_t __attribute__((aligned(32))) pvrlineverts[4];

void draw_pvr_line(vector_t *v1, vector_t *v2, int color)
{
	vector_t *ov1;
	vector_t *ov2;
	pvr_vertex_t *vert;
	float hlw_invmag;
	float dx,dy;
	float nx,ny;

	if (v1->x <= v2->x) {
		ov1 = v1;
		ov2 = v2;
	} else {
		ov1 = v2;
		ov2 = v1;
	}

	// https://devcry.heiho.net/html/2017/20170820-opengl-line-drawing.html
	dx = ov2->x - ov1->x;
	dy = ov1->y - ov2->y;//ov2->y - ov1->y;
	// I have been *told* that it is better for codegen to not call `frsqrt` 
	hlw_invmag = (1.0f / sqrtf((dx * dx) + (dy * dy))) * (LINEWIDTH * 0.5f);
	nx = dy * hlw_invmag;//-dy * hlw_invmag;
	ny = dx * hlw_invmag;

	vert = pvrlineverts;
	vert->flags = PVR_CMD_VERTEX;
	vert->x = ov1->x + nx;
	vert->y = ov1->y + ny;
	vert->z = ov1->z;
	vert++->argb = color;

	vert->flags = PVR_CMD_VERTEX;
	vert->x = ov1->x - nx;
	vert->y = ov1->y - ny;
	vert->z = ov2->z;
	vert++->argb = color;

	vert->flags = PVR_CMD_VERTEX;
	vert->x = ov2->x + nx;
	vert->y = ov2->y + ny;
	vert->z = ov1->z;
	vert++->argb = color;

	vert->flags = PVR_CMD_VERTEX_EOL;
	vert->x = ov2->x - nx;
	vert->y = ov2->y - ny;
	vert->z = ov2->z;
	vert->argb = color;

//	sq_fast_cpy(SQ_MASK_DEST(PVR_TA_INPUT), pvrlineverts, 4);
	pvr_list_prim(PVR_LIST_TR_POLY, pvrlineverts, 4 * sizeof(pvr_vertex_t));
}

void compute_norm(float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2, float *nx, float *ny, float *nz) {
		float u[3];
		u[0] = (-1.0f * x0) + (x1);
		u[1] = (-1.0f * y0) + (y1);
		u[2] = (-1.0f * z0) + (z1);

		float v[3];
		v[0] = (-1.0f * x0) + (x2);
		v[1] = (-1.0f * y0) + (y2);
		v[2] = (-1.0f * z0) + (z2);

		float uvi, uvj, uvk;
		uvi = u[1] * v[2] - v[1] * u[2];
		uvj = v[0] * u[2] - u[0] * v[2];
		uvk = u[0] * v[1] - v[0] * u[1];
		*nx = uvi;
		*ny = uvj;
		*nz = uvk;
}

uint8_t __attribute__((aligned(32))) tr_buf[2048576*2];

pvr_vertex_t __attribute__((aligned(32))) verts[4];
pvr_poly_cxt_t ccxt;
pvr_poly_hdr_t __attribute__((aligned(32))) chdr;

pvr_ptr_t tex;
static pvr_vertex_t skypic_verts[4];

float xangle = 0;
float yangle = 0;
float zangle = 0;

int main(int argc, char **argv)
{
	pvr_init(&pvr_params);

	pvr_set_vertbuf(PVR_LIST_TR_POLY, tr_buf, 2048576*2);

	pvr_ptr_t pvrsky = pvr_mem_malloc(256 * 256 * 2);
	pvr_poly_cxt_t pvrskycxt;
	pvr_poly_hdr_t pvrskyhdr;

	pvr_poly_cxt_txr(&pvrskycxt, PVR_LIST_TR_POLY, PVR_TXRFMT_ARGB1555 | PVR_TXRFMT_TWIDDLED, 256, 256, pvrsky, PVR_FILTER_BILINEAR);
	pvrskycxt.depth.write = PVR_DEPTHWRITE_DISABLE;
	pvrskycxt.txr.uv_flip = PVR_UVFLIP_NONE;
	pvr_poly_compile(&pvrskyhdr, &pvrskycxt);
	pvr_txr_load(vrbg__data, pvrsky, 256*256*2);
#if 0
        uint16_t *dst16 = (uint16_t *)pvrsky;
	uint16_t *src16 = (uint16_t *)vrbg__data;
	for (int i=0;i<131072;i++) {
//		dbgio_printf("src16[i] == %04x\n", src16[i]);
		dst16[i] = src16[i];
	}
#endif
	skypic_verts[0].flags = PVR_CMD_VERTEX;
	skypic_verts[0].x = 0;
	skypic_verts[0].y = 480;
	skypic_verts[0].z = 0.0000011f;
	skypic_verts[0].argb = 0xffeeeeee;

	skypic_verts[1].flags = PVR_CMD_VERTEX;
	skypic_verts[1].x = 0;
	skypic_verts[1].y = 0;
	skypic_verts[1].z = 0.0000011f;
	skypic_verts[1].argb = 0xffeeeeee;

	skypic_verts[2].flags = PVR_CMD_VERTEX;
	skypic_verts[2].x = 640;
	skypic_verts[2].y = 480;
	skypic_verts[2].z = 0.0000011f;
	skypic_verts[2].argb = 0xffeeeeee;

	skypic_verts[3].flags = PVR_CMD_VERTEX_EOL;
	skypic_verts[3].x = 640;
	skypic_verts[3].y = 0;
	skypic_verts[3].z = 0.0000011f;
	skypic_verts[3].argb = 0xffeeeeee;

	R_Frustum(R_ProjectionMatrix, -16.0f, 16.0f, -9.0f, 9.0f, 8.0f, 16384.0f, 1.0f);
	R_Viewport(R_ViewportMatrix, 0, 0, 640, 480);

	pvr_set_bg_color(0, 32.0f/255.0f, 1.0f);
	tex = pvr_mem_malloc(32*32);
	uint16_t *tex16 = (uint16_t *)tex;
	for (int i=0;i<32*32;i++) {
		tex16[i] = 0xFFFF;
	}
#define PVR_MIN_Z 0.000001f
	pvr_set_zclip(PVR_MIN_Z);

	pvr_poly_cxt_txr(&ccxt, PVR_LIST_TR_POLY, PVR_TXRFMT_ARGB1555 | PVR_TXRFMT_TWIDDLED, 32, 32, tex, PVR_FILTER_NONE);

	ccxt.gen.shading = PVR_SHADE_GOURAUD;
	ccxt.gen.culling = PVR_CULLING_NONE;

	pvr_poly_compile(&chdr, &ccxt);

	verts[0].flags = PVR_CMD_VERTEX;
	verts[1].flags = PVR_CMD_VERTEX;
	verts[2].flags = PVR_CMD_VERTEX_EOL;
	verts[3].flags = PVR_CMD_VERTEX_EOL;

	verts[0].u = 0; verts[0].v = 0;
	verts[1].u = 1; verts[1].v = 0;
	verts[2].u = 1; verts[2].v = 1;

	camerax = (models_bounding_boxes[0][1] + models_bounding_boxes[0][0]) / 2;
	cameraz = (models_bounding_boxes[0][3] + models_bounding_boxes[0][2]) / 2;

	float next_camera_y = -1.0f;
	cameray = 200;

	int last_segment = 0;

	float last_camerax = camerax;
	float last_cameraz = cameraz;

	while(1) {
		maple_device_t *controller;
		cont_state_t *cont;
		controller = maple_enum_type(0, MAPLE_FUNC_CONTROLLER);

		if (controller) {
			cont = maple_dev_status(controller);

			if ((cont->buttons & CONT_START) && cont->ltrig && cont->rtrig)
				exit(0);

			if (cont->buttons & CONT_Y) {
				cameraz -= (200 * sinf(viewangle));
				camerax += (200 * cosf(viewangle));
				//cameray += (100 * sinf(viewpitch));
			}

			if (cont->buttons & CONT_A) {
				cameraz += (200 * sinf(viewangle));
				camerax -= (200 * cosf(viewangle));
				//cameray -= (100 * sinf(viewpitch));
			}

			if (cont->buttons & CONT_X)  {
				//camerax -=50;
			}
			if (cont->buttons & CONT_B) {
				// camerax +=50;
			}

//			if (cont->buttons & CONT_A) lightz-=100;
//			if (cont->buttons & CONT_B) lightz+=100;
//			if (cont->buttons & CONT_X) {camerax = cameray = cameraz = 0;}
//			if (cont->joyy < 0) viewpitch -= 0.01f;
//			if (cont->joyy > 0) viewpitch += 0.01f;
			if (cont->joyx < 0) viewangle -= 0.03f;
			if (cont->joyx > 0) viewangle += 0.03f;
//			if (cont->ltrig) cameray -=20;
//			if (cont->rtrig) cameray +=20;
		}

		R_RotateX(R_RotX, sinf(viewpitch+/*0.9*/1.05), cosf(viewpitch+/*0.9*/1.05));
		R_RotateY(R_RotY, sinf(viewangle), cosf(viewangle));
		R_Translate(R_Tran, -camerax, -cameray, -cameraz);
		mat_load(&R_ViewportMatrix);
		mat_apply(&R_ProjectionMatrix);
		mat_apply(&R_RotX);
		mat_apply(&R_RotY);
		mat_apply(&R_Tran);

		if (next_camera_y > 0) {
			cameray = next_camera_y;
		}
		float drawangle = viewangle;
		while(drawangle < 0) drawangle += F_PI * 2.0f;
		while(drawangle > (F_PI * 2.0f)) drawangle += F_PI * 2.0f;
		float ang = 0 - ((int)((drawangle / (F_PI * 2.0f))*255.0f) & 255);
		float u0, v0, u1, v1;
		u0 = (float)ang / 256.0f;
		u1 = u0 + 1.0f;
		v0 = 1.0 / 256.0f; // 0.5f / 128.0f;
		v1 = 191.0f / 256.0f;

		skypic_verts[0].u = u0;
		skypic_verts[0].v = v1;
		skypic_verts[1].u = u0;
		skypic_verts[1].v = v0;
		skypic_verts[2].u = u1;
		skypic_verts[2].v = v1;
		skypic_verts[3].u = u1;
		skypic_verts[3].v = v0;

		pvr_wait_ready();
		pvr_scene_begin();
		pvr_list_prim(PVR_LIST_TR_POLY, &pvrskyhdr, sizeof(pvr_poly_hdr_t));
		pvr_list_prim(PVR_LIST_TR_POLY, skypic_verts, sizeof(pvr_vertex_t)*4);
		pvr_list_prim(PVR_LIST_TR_POLY, &chdr, sizeof(pvr_poly_hdr_t));

                int current_segment_number = -1;
		for (int j=0;j<NUM_MODELS;j++) {
			if (models_bounding_boxes[j][0] <= camerax && camerax <= models_bounding_boxes[j][1]) {
				if (models_bounding_boxes[j][2] <= cameraz && cameraz <= models_bounding_boxes[j][3]) {
					current_segment_number = j;
					break;
				}
			}
		}

		if (current_segment_number == -1) {
			camerax = last_camerax;//(models_bounding_boxes[last_segment][1] + models_bounding_boxes[last_segment]][0]) / 2;
			cameraz = last_cameraz;//(models_bounding_boxes[last_segment][3] + models_bounding_boxes[last_segment][2]) / 2;

		} else {
			last_segment = current_segment_number;
		}

		Point3D P = {camerax, cameray, cameraz};  // Point above the triangle

		dbgio_printf("camera is in segment %d\n", current_segment_number);

		int over_face_found = 0;
		next_camera_y = -1.0f;
		for (int j=0;j<NUM_MODELS;j++) {

		int NUM_FACES = face_counts[j];

		for (int i=0;i<NUM_FACES;i++) {
			int over_this_face = 0;

			int v0 = models_faces[j][i][0];
			int v1 = models_faces[j][i][1];
			int v2 = models_faces[j][i][2];

			float w0,w1,w2;

			float x0 = models_vertices[j][v2][0];
			float y0 = models_vertices[j][v2][1] * 0.4f;
			float z0 = models_vertices[j][v2][2];


			float x1 = models_vertices[j][v1][0];
			float y1 = models_vertices[j][v1][1] * 0.4f;
			float z1 = models_vertices[j][v1][2];

			float x2 = models_vertices[j][v0][0];
			float y2 = models_vertices[j][v0][1] * 0.4f;
			float z2 = models_vertices[j][v0][2];
//			if (!over_face_found) {
				Point3D A = {x0, y0, z0};
				Point3D B = {x1, y1, z1};
				Point3D C = {x2, y2, z2};
				over_this_face = isPointInsideTriangle(A, B, C, P);
				over_face_found = over_this_face;
				if (next_camera_y == -1.0f) {
					if (over_this_face) 
						next_camera_y = fmaxf(y0, fmaxf(y1, y2)) + 50;
				} else {
					if (over_this_face) {

					if ((fmaxf(y0, fmaxf(y1, y2)) + 50) < next_camera_y) {
						next_camera_y = fmaxf(y0, fmaxf(y1, y2)) + 50;

					}
					}
				}
//			}

			mat_trans_single3_nodivw(x0,y0,z0,w0);
			mat_trans_single3_nodivw(x1,y1,z1,w1);
			mat_trans_single3_nodivw(x2,y2,z2,w2);

			uint32_t fc,fc2;

			if (over_this_face) {
				fc = 0xffff0000;
				fc2 = 0xffff0000;
			} else {
			int ci = models_colors[j][i][0];
			float r = (float)colors[ci][0] / 255.0f;
			float g = (float)colors[ci][1] / 255.0f;
			float b = (float)colors[ci][2] / 255.0f;
			fc = PVR_PACK_COLOR((ci != 0),r,g,b);

			int ci2 = models_colors[j][i][1];
			float r2 = (float)colors[ci2][0] / 255.0f;
			float g2 = (float)colors[ci2][1] / 255.0f;
			float b2 = (float)colors[ci2][2] / 255.0f;
			fc2 = PVR_PACK_COLOR((ci2 != 0),r2,g2,b2);
			}

			uint32_t vismask = (z0 > -w0) | ((z1 > -w1) << 1) | ((z2 > -w2) << 2);
			float x3,y3,z3,w3;
			int usespare = 0;
			int nodraw = 0;
			if (vismask == 0) {
				nodraw = 1;
			} else if (vismask == 7) {
				goto sendit;
			} else {
				switch (vismask) {
				case 1:
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1,&x1,&y1,&z1,&w1);
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2,&x2,&y2,&z2,&w2);
				break;
				case 2:
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1,&x0,&y0,&z0,&w0);
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2,&x2,&y2,&z2,&w2);
				break;
				case 3:
				usespare = 1;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2,&x3,&y3,&z3,&w3);
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2,&x2,&y2,&z2,&w2);
				break;
				case 4:
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2,&x0,&y0,&z0,&w0);
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2,&x1,&y1,&z1,&w1);
				break;
				case 5:
				usespare = 1;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2,&x3,&y3,&z3,&w3);
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1,&x1,&y1,&z1,&w1);
				break;
				case 6:
				usespare = 1;
				x3 = x2; y3 = y2; z3 = z2; w3 = w2;
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2,&x2,&y2,&z2,&w2);
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1,&x0,&y0,&z0,&w0);
				break;
				}
			}

sendit:
			if (nodraw) goto endloop;
			if (!usespare) {
				verts[2].flags = PVR_CMD_VERTEX_EOL;
			} else {
				verts[2].flags = PVR_CMD_VERTEX;
			}

			perspdiv(&x0,&y0,&z0,w0);
			perspdiv(&x1,&y1,&z1,w1);
			perspdiv(&x2,&y2,&z2,w2);

		 	if (usespare) {
				perspdiv(&x3,&y3,&z3,w3);
			}

			verts[0].x = (x0);
			verts[0].y = (y0);
			verts[0].z = (z0);
			verts[0].argb = fc2;

			verts[1].x = (x1);
			verts[1].y = (y1);
			verts[1].z = (z1);
			verts[1].argb = fc;

			verts[2].x = (x2);
			verts[2].y = (y2);
			verts[2].z = (z2);
			verts[2].argb = fc2;

			if (usespare) {
				verts[3].x = (x3);
				verts[3].y = (y3);
				verts[3].z = (z3);
				verts[3].argb = fc;
			}

			pvr_list_prim(PVR_LIST_TR_POLY, &verts[0], sizeof(pvr_vertex_t) * (3 + usespare));

#if 0
			vector_t a;
			vector_t b;
			a.x = verts[0].x;
			a.y = verts[0].y;
			a.z = 5;
			b.x = verts[1].x;
			b.y = verts[1].y;
			b.z = 5;
			draw_pvr_line(&a,&b,0xffffffff);

			a.x = verts[1].x;
			a.y = verts[1].y;
			a.z = 5;
			b.x = verts[2].x;
			b.y = verts[2].y;
			b.z = 5;
			draw_pvr_line(&a,&b,0xffffffff);

			a.x = verts[2].x;
			a.y = verts[2].y;
			a.z = 5;
			b.x = verts[0].x;
			b.y = verts[0].y;
			b.z = 5;
			draw_pvr_line(&a,&b,0xffffffff);

#endif

endloop:
		} // for i num faces
		} // for j num models

		pvr_scene_finish();
	}
	return 0;
}


#if 0

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <kos.h>
#include <kos/init.h>

#include <dc/maple.h>
#include <dc/vmu_fb.h>
#include <dc/vmu_pkg.h>


KOS_INIT_FLAGS(INIT_DEFAULT);

static unsigned short vmu_icon_pal[] = {
0xF404, 0xF303, 0xF202, 0xF101, 0xFF62, 0xFE41, 0xFB42, 0xF732, 0xFC31, 0xF610, 0xF721, 0xF921, 0xF842, 0xFA21, 0xFA55, 0xF510, };

static unsigned char vmu_icon_img[] = {
0x23,0x32,0x21,0x22,0x2B,0x8B,0x31,0x22,0x22,0x23,0xAA,0xF3,0x33,0x22,0x23,0x21,
0x23,0x93,0x32,0x2F,0xBD,0xA3,0x33,0x33,0x33,0x33,0xFD,0x8A,0xF3,0x23,0x3F,0x31,
0x23,0xF7,0x73,0x2A,0xDD,0x33,0xCC,0xAA,0xAA,0x73,0x33,0x8D,0x93,0x3F,0x9F,0x31,
0x23,0xFF,0xA7,0x78,0xDA,0xCC,0xAF,0xFF,0xFF,0x9C,0x79,0xB8,0x8F,0xF9,0xAF,0x21,
0x22,0x3F,0xA7,0x9D,0xDA,0x99,0xFF,0x9F,0xFF,0x9F,0xF9,0x78,0x8A,0x99,0xF3,0x20,
0x21,0x23,0x99,0xA8,0xDF,0x9F,0xFF,0x33,0x33,0x99,0x99,0x98,0x56,0x7A,0xF2,0x00,
0x21,0x13,0xF9,0x9B,0xBF,0xF3,0x33,0x33,0x33,0x33,0x39,0xFB,0x8B,0x99,0x32,0x11,
0x21,0x12,0x37,0xA8,0xDF,0xFF,0x33,0x39,0xF3,0x33,0x97,0x98,0x8B,0x99,0x31,0x11,
0x22,0x12,0x33,0xCA,0x8A,0xFA,0x9F,0x3E,0xD3,0xFC,0x7A,0xA5,0xB9,0xA2,0x21,0x11,
0x22,0x23,0xE7,0xAD,0xAD,0xFF,0x99,0xFD,0xDA,0xC7,0x9A,0xDC,0xB7,0x79,0x21,0x11,
0x12,0x33,0xCA,0xF8,0x5A,0xF9,0x99,0xFD,0xD9,0xAA,0xAA,0x74,0x89,0x77,0x22,0x12,
0x12,0x3C,0x7A,0xFD,0x88,0xA9,0xAA,0xFB,0xDF,0xAB,0xAC,0x68,0x8A,0x77,0x92,0x22,
0x12,0x3E,0x7A,0x39,0xA5,0x8A,0xBB,0x9B,0xD9,0xBB,0xB8,0x4B,0xB2,0x77,0xA2,0x22,
0x12,0x3E,0x79,0x3F,0xA5,0xDB,0xBD,0x98,0xD9,0xDD,0xD6,0x4B,0x92,0xE7,0xA2,0x12,
0x22,0x3E,0xA9,0x33,0x98,0x8A,0xBB,0x7B,0x8F,0xB6,0xC8,0x4B,0x22,0xE7,0xB2,0x01,
0x23,0x3C,0x79,0x33,0x9A,0x58,0xDA,0xA8,0x5F,0x76,0x44,0x6A,0x22,0xE7,0xB2,0x01,
0x32,0x37,0x79,0x37,0xCF,0x54,0x5D,0xB5,0x5A,0x84,0x44,0xAA,0x62,0xE7,0x72,0x11,
0x22,0x37,0x77,0xC6,0xBC,0x98,0x44,0x85,0x48,0x44,0x48,0xA6,0xBE,0xC7,0xA2,0x01,
0x22,0x3E,0x77,0xBD,0x8C,0xFF,0xB4,0x45,0x54,0x45,0xFF,0xB5,0x66,0xE7,0xA2,0x10,
0x23,0xEE,0x7B,0xB8,0x8A,0x9A,0x9B,0x48,0x54,0x89,0xAF,0x85,0x8C,0xC7,0x7E,0x21,
0x3E,0xE7,0x7A,0xAA,0xAA,0xA8,0xAF,0xD8,0xD8,0xFD,0x8D,0x66,0x6B,0x77,0x7C,0xE2,
0xAB,0x9A,0x9B,0xB9,0xF9,0x9C,0xDD,0x8F,0xA5,0x88,0xA9,0x99,0x76,0xB9,0x99,0x9C,
0x33,0x33,0x3A,0xDD,0x63,0x3C,0xAB,0x95,0x4A,0xDA,0xF2,0x2C,0x66,0xA2,0x21,0x22,
0x11,0x11,0x22,0xA8,0x8C,0x33,0xC9,0xDD,0xD5,0xAB,0x22,0x74,0x5D,0x21,0x00,0x00,
0x01,0x11,0x01,0x2D,0x55,0xE7,0x89,0xF9,0xA9,0xA8,0x9C,0x45,0xD2,0x10,0x00,0x00,
0x00,0x10,0x00,0x12,0x85,0x58,0xB5,0x8A,0xBD,0x5A,0x45,0x48,0x21,0x01,0x11,0x00,
0x00,0x00,0x00,0x12,0x2D,0x54,0x8A,0x5B,0xD5,0x88,0x45,0xD3,0x21,0x11,0x11,0x00,
0x00,0x00,0x00,0x00,0x12,0x2D,0x8A,0x5D,0x85,0xB8,0xD3,0x32,0x22,0x21,0x11,0x11,
0x10,0x11,0x11,0x00,0x01,0x23,0x33,0xAD,0xDA,0x33,0x32,0x22,0x33,0x21,0x01,0x22,
0x11,0x11,0x11,0x00,0x10,0x12,0x12,0xF9,0x9F,0x31,0x12,0x22,0x22,0x11,0x01,0x22,
0x10,0x10,0x00,0x00,0x00,0x11,0x12,0x39,0x93,0x32,0x22,0x22,0x22,0x22,0x22,0x22,
0x10,0x00,0x00,0x00,0x00,0x00,0x11,0x3F,0xF3,0x22,0x32,0x22,0x22,0x22,0x32,0x22,
};

static vmu_pkg_t pkg;

typedef struct doom64_settings_s {
	int version;
	int HUDopacity;
	int SfxVolume;
	int MusVolume;
	int brightness;
	int enable_messages;
	int M_SENSITIVITY;
	int MotionBob;
	int Rumble;
	int VideoFilter;
	int Autorun;
	int runintroduction;
	int StoryText;
	int MapStats;
	int HUDmargin;
	int ColoredHUD;
	int Quality;
	int FpsUncap;
	int PlayDeadzone;
	int Interpolate;
	int VmuDisplay;
} doom64_settings_t;

static char full_fn[512];

static char *get_vmu_fn(maple_device_t *vmudev, char *fn) {
	if (fn)
		sprintf(full_fn, "/vmu/%c%d/%s", 'a'+vmudev->port, vmudev->unit, fn);
	else
		sprintf(full_fn, "/vmu/%c%d", 'a'+vmudev->port, vmudev->unit);

	return full_fn;
}

int I_SavePakSettings(doom64_settings_t *msettings)
{
	uint8 *pkg_out;
	ssize_t pkg_size;
	maple_device_t *vmudev = NULL;

	vmudev = maple_enum_type(0, MAPLE_FUNC_MEMCARD);
	if (!vmudev)
		return -1;

	file_t d = fs_open(get_vmu_fn(vmudev, "TEST_SETTINGS"), O_WRONLY | O_CREAT);
	if (!d)
		return -2;

	memset(&pkg, 0, sizeof(vmu_pkg_t));
	strcpy(pkg.desc_short,"TEST settings");
	strcpy(pkg.desc_long, "TEST settings data");
	strcpy(pkg.app_id, "TEST");
	pkg.icon_cnt = 1;
	pkg.icon_data = vmu_icon_img;
	memcpy(pkg.icon_pal, vmu_icon_pal, sizeof(vmu_icon_pal));
	pkg.data_len = 384;
	// doesn't matter, just not NULL
	pkg.data = vmu_icon_img;

	vmu_pkg_build(&pkg, &pkg_out, &pkg_size);

	if (!pkg_out || pkg_size <= 0) {
		fs_close(d);
		return -2;
	}

	memcpy(&pkg_out[640], msettings, sizeof(doom64_settings_t));

	ssize_t rv = fs_write(d, pkg_out, pkg_size);
	if (rv < 0) {
		fs_close(d);
		return -2;
	}
	ssize_t total = rv;
	while (total < pkg_size) {
		rv = fs_write(d, pkg_out + total, pkg_size - total);
		if (rv < 0) {
			fs_close(d);
			return -2;
		}
		total += rv;
	}

	fs_close(d);
	free(pkg_out);

	if (rv == pkg_size)
		return 0;
	else
		return -2;
}

int I_ReadPakSettings(doom64_settings_t *msettings)
{
	ssize_t size;
	maple_device_t *vmudev = NULL;
	uint8_t *data;

	vmudev = maple_enum_type(0, MAPLE_FUNC_MEMCARD);
	if (!vmudev)
		return -1;

	file_t d = fs_open(get_vmu_fn(vmudev, "TEST_SETTINGS"), O_RDONLY);
	if (!d)
		return -2;

	size = fs_total(d);
	data = calloc(1, size);

	if (!data) {
		fs_close(d);
		return -2;
	}

	// read version first
	ssize_t res = fs_read(d, data, size);
	if (res < 0) {
		fs_close(d);
		return -2;
	}
	ssize_t total = res;
	while (total < size) {
		res = fs_read(d, data + total, size - total);
		if (res < 0) {
			fs_close(d);
			return -2;
		}
		total += res;
	}

	fs_close(d);


	if (res != total) {
		free(data);
		return -2;
	}

	int save_size = sizeof(doom64_settings_t);
	memcpy(msettings, &data[640], save_size);

	free(data);

	return 0;
}

int main(int argc, char *argv[]) {

    pvr_init_defaults();

    doom64_settings_t ms;
    doom64_settings_t ms2;

    memset(&ms2, 0, sizeof(doom64_settings_t));
    uint8_t *ms8 = (uint8_t *)&ms;
    for (int i=0;i<sizeof(doom64_settings_t);i++) {
        ms8[i] = (uint8_t)(i % 16);
    }

    int rv = I_SavePakSettings(&ms);
    if (rv != 0) {
        dbgio_printf("error saving settings\n");
        return 0;
    }

    int rv2 = I_ReadPakSettings(&ms2);
    if (rv2 != 0) {
        dbgio_printf("error loading settings\n");
        return 0;
    }
    uint8_t *ms28 = (uint8_t *)&ms2;
    for (int i=0;i<sizeof(doom64_settings_t);i++) {
        if (ms28[i] != ms8[i]) {
            dbgio_printf("settings byte %d mismatch (save: %x, load %x)\n", i, ms8[i], ms28[i]);
            return 0;
        }
    }

    dbgio_printf("successful settings save and load\n");

    return 0;
}

#endif


