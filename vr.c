#include <kos.h>
#include <stdlib.h>

#include "vrbg.h"

//#include "pal1.h"
//#include "pal2.h"
#include "pal3.h"

//#include "beginner.h"
//#include "medium.h"
#include "expert.h"

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
    float dnn = dotProduct(normal, normal);
    if (dnn == 0.0f) dnn = 0.0001f;
    float t = -(dotProduct(normal, P) + D) / dnn;

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
Matrix R_RotZ;
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

static float xout,yout,zout,wout;

static void nearz_clip(float x0,float y0,float z0,float w0,float x1,float y1, float z1,float w1)
{
	const float d0 = w0 + z0;
	const float d1 = w1 + z1;

	float d1subd0 = d1 - d0;
	if (d1subd0 == 0.0f) d1subd0 = 0.0001f;

	float t = (fabsf(d0) * (1.0f / sqrtf(d1subd0 * d1subd0))) + 0.000001f;
	float invt = 1.0f - t;

	wout = lerp(w0, w1);
	xout = lerp(x0, x1);
	yout = lerp(y0, y1);
	zout = lerp(z0, z1);
}


static inline void perspdiv(float *x, float *y, float *z, float w)
{
	if (w == 0) w = 0.0001f;
	float invw = 1.0f / w;
	*x = (*x * invw);
	*y = (-*y *invw);

	if (w == 1.0f) {
		float zz = *z;
		if (zz == -1.0001f) zz = -1.0002f;
		*z = 1.0f / (1.0001f + zz);
	}
	else {
		*z = invw;
	}
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
pvr_init_params_t pvr_params = { { PVR_BINSIZE_16, 0, 0, 0, 0 },
				2048576,
				0, // dma disabled
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
	pvr_prim(pvrlineverts, 4 * sizeof(pvr_vertex_t));
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

//uint8_t __attribute__((aligned(32))) tr_buf[2048576*2];

pvr_vertex_t __attribute__((aligned(32))) verts[4];
pvr_poly_cxt_t ccxt;
pvr_poly_hdr_t __attribute__((aligned(32))) chdr;

pvr_ptr_t tex;
static pvr_vertex_t skypic_verts[4];

float xangle = 0;
float yangle = 0;
float zangle = 0;
pvr_dr_state_t dr_state;

float lastnx = -100000.0f;
float lastny = -100000.0f;
float lastnz = -100000.0f;

int main(int argc, char **argv)
{
	pvr_init(&pvr_params);

//	pvr_set_vertbuf(PVR_LIST_TR_POLY, tr_buf, 2048576*2);

	pvr_ptr_t pvrsky = pvr_mem_malloc(256 * 256 * 2);
	pvr_poly_cxt_t pvrskycxt;
	pvr_poly_hdr_t pvrskyhdr;

	pvr_poly_cxt_txr(&pvrskycxt, PVR_LIST_OP_POLY, PVR_TXRFMT_ARGB1555 | PVR_TXRFMT_TWIDDLED, 256, 256, pvrsky, PVR_FILTER_BILINEAR);
	pvrskycxt.depth.write = PVR_DEPTHWRITE_DISABLE;
	pvrskycxt.txr.uv_flip = PVR_UVFLIP_NONE;
	pvr_poly_compile(&pvrskyhdr, &pvrskycxt);
	pvr_txr_load(vrbg__data, pvrsky, 256*256*2);

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

	R_Frustum(R_ProjectionMatrix, -18.0f, 18.0f, -12.0f, 12.0f, 8.0f, 16384.0f, 1.0f);
	R_Viewport(R_ViewportMatrix, 0, 0, 640, 480);

	pvr_set_bg_color(0, 32.0f/255.0f, 1.0f);
	tex = pvr_mem_malloc(32*32);
	uint16_t *tex16 = (uint16_t *)tex;
	for (int i=0;i<32*32;i++) {
		tex16[i] = 0xFFFF;
	}
#define PVR_MIN_Z 0.000001f
	pvr_set_zclip(PVR_MIN_Z);

	pvr_poly_cxt_txr(&ccxt, PVR_LIST_OP_POLY, PVR_TXRFMT_ARGB1555 | PVR_TXRFMT_TWIDDLED, 32, 32, tex, PVR_FILTER_NONE);

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
				cameraz -= (50 * sinf(viewangle));
				camerax += (50 * cosf(viewangle));
				//cameray += (100 * sinf(viewpitch));
			}

			if (cont->buttons & CONT_A) {
				cameraz += (50 * sinf(viewangle));
				camerax -= (50 * cosf(viewangle));
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

		if (viewangle < 0) viewangle += F_PI * 2.0f;
		if (viewangle > F_PI * 2.0f) viewangle -= F_PI * 2.0f;

		int do_rotz = 0;
		if (lastnx != -100000.0f) {
			float rotangle = atan2f(lastny, sqrtf((lastnx*lastnx)+(lastnz*lastnz))) - viewangle;
			if (rotangle < 0) rotangle += F_PI * 2.0f;
			if (rotangle > F_PI * 2.0f) rotangle -= F_PI * 2.0f;

			dbgio_printf("nx %f ny %f nz %f\nviewangle %f\nrotangle %f\n", lastnx, lastny, lastnz, viewangle, rotangle);

			R_RotateZ(R_RotZ, sinf(rotangle), cosf(rotangle));
			do_rotz = 1;
		}

		R_RotateX(R_RotX, sinf(viewpitch+/*0.9*/1.1), cosf(viewpitch+/*0.9*/1.1));
		R_RotateY(R_RotY, sinf(viewangle), cosf(viewangle));
		R_Translate(R_Tran, -camerax, -cameray, -cameraz);
//		if (do_rotz)  {
//			mat_load(&R_RotZ);
//		mat_apply(&R_ViewportMatrix);
//		} else {
		mat_load(&R_ViewportMatrix);

//		}
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
		pvr_list_begin(PVR_LIST_OP_POLY);
		pvr_dr_init(&dr_state);

		pvr_prim(&pvrskyhdr, sizeof(pvr_poly_hdr_t));
		pvr_prim(skypic_verts, sizeof(pvr_vertex_t)*4);
		pvr_prim(&chdr, sizeof(pvr_poly_hdr_t));

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

		last_camerax = camerax;
		last_cameraz = cameraz;

//		dbgio_printf("camera is in segment %d\n", current_segment_number);

		int over_face_found = 0;
		next_camera_y = -1.0f;

//		float cnx = (lastnx != -100000.0f) ? lastnx : 0.0f;
//		float cny = (lastny != -100000.0f) ? lastny : 0.0f;
//		float cnz = (lastnz != -100000.0f) ? lastnz : 0.0f;

		for (int j=0;j<NUM_MODELS;j++) {

		int NUM_FACES = face_counts[j];

		for (int i=0;i<NUM_FACES;i++) {
			int over_this_face = 0;

			int v0 = models_faces[j][i][0];
			int v1 = models_faces[j][i][1];
			int v2 = models_faces[j][i][2];

//			cnx *= 20;
//			cny *= 20;
//			cnz *= 20;

			float w0,w1,w2;

			float x0 = models_vertices[j][v2][0];
			float y0 = (models_vertices[j][v2][1] * 0.4f);
			float z0 = models_vertices[j][v2][2];


			float x1 = models_vertices[j][v1][0];
			float y1 = (models_vertices[j][v1][1] * 0.4f);
			float z1 = models_vertices[j][v1][2];

			float x2 = models_vertices[j][v0][0];
			float y2 = (models_vertices[j][v0][1] * 0.4f);
			float z2 = models_vertices[j][v0][2];
//			if (!over_face_found) {
				Point3D A = {x0, y0, z0};
				Point3D B = {x1, y1, z1};
				Point3D C = {x2, y2, z2};
				over_this_face = isPointInsideTriangle(A, B, C, P);
				over_face_found = over_this_face;
				if (next_camera_y == -1.0f) {
					if (over_this_face) {
						next_camera_y = fmaxf(y0, fmaxf(y1, y2)) + 50;
						lastnx = models_normals[j][i][0];
						lastny = models_normals[j][i][1];
						lastnz = models_normals[j][i][2];
						if (lastny < 0) {
							lastnx = -lastnx;
							lastny = -lastny;
							lastnz = -lastnz;
						}
						if (lastnx == -0.0f) lastnx = 0.0f;
						if (lastny == -0.0f) lastny = 0.0f;
						if (lastnz == -0.0f) lastnz = 0.0f;
					}
				} else {
					if (over_this_face) {

					if ((fmaxf(y0, fmaxf(y1, y2)) + 50) < next_camera_y) {
						next_camera_y = fmaxf(y0, fmaxf(y1, y2)) + 50;

						lastnx = models_normals[j][i][0];
						lastny = models_normals[j][i][1];
						lastnz = models_normals[j][i][2];
						if (lastny < 0) {
							lastnx = -lastnx;
							lastny = -lastny;
							lastnz = -lastnz;
						}

						if (lastnx == -0.0f) lastnx = 0.0f;
						if (lastny == -0.0f) lastny = 0.0f;
						if (lastnz == -0.0f) lastnz = 0.0f;
					}
					}
				}
//			}

			mat_trans_single3_nodivw(x0,y0,z0,w0);
			mat_trans_single3_nodivw(x1,y1,z1,w1);
			mat_trans_single3_nodivw(x2,y2,z2,w2);

			uint32_t fc,fc2;

			{
//			if (over_this_face) {
//				fc = 0xffff0000;
//				fc2 = 0xffff0000;
//			} else {
			int ci = models_colors[j][i][0];
			float r = (float)colors[ci][0] / 255.0f;
			float g = (float)colors[ci][1] / 255.0f;
			float b = (float)colors[ci][2] / 255.0f;
			fc2 = PVR_PACK_COLOR((ci != 0),r,g,b);

			int ci2 = models_colors[j][i][1];
			float r2 = (float)colors[ci2][0] / 255.0f;
			float g2 = (float)colors[ci2][1] / 255.0f;
			float b2 = (float)colors[ci2][2] / 255.0f;
			fc = PVR_PACK_COLOR((ci2 != 0),r2,g2,b2);
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
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1);//,&x1,&y1,&z1,&w1);
				x1=xout;y1=yout;z1=zout;w1=wout;
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2);//,&x2,&y2,&z2,&w2);
				x2=xout;y2=yout;z2=zout;w2=wout;
				break;
				case 2:
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1);//&x0,&y0,&z0,&w0);
				x0=xout;y0=yout;z0=zout;w0=wout;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2);//,&x2,&y2,&z2,&w2);
				x2=xout;y2=yout;z2=zout;w2=wout;
				break;
				case 3:
				usespare = 1;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2);//,&x3,&y3,&z3,&w3);
				x3=xout;y3=yout;z3=zout;w3=wout;
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2);//,&x2,&y2,&z2,&w2);
				x2=xout;y2=yout;z2=zout;w2=wout;
				break;
				case 4:
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2);//,&x0,&y0,&z0,&w0);
				x0=xout;y0=yout;z0=zout;w0=wout;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2);//,&x1,&y1,&z1,&w1);
				x1=xout;y1=yout;z1=zout;w1=wout;
				break;
				case 5:
				usespare = 1;
				nearz_clip(x1,y1,z1,w1,x2,y2,z2,w2);//,&x3,&y3,&z3,&w3);
				x3=xout;y3=yout;z3=zout;w3=wout;
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1);//,&x1,&y1,&z1,&w1);
				x1=xout;y1=yout;z1=zout;w1=wout;
				break;
				case 6:
				usespare = 1;
				x3 = x2; y3 = y2; z3 = z2; w3 = w2;
				nearz_clip(x0,y0,z0,w0,x2,y2,z2,w2);//,&x2,&y2,&z2,&w2);
				x2=xout;y2=yout;z2=zout;w2=wout;
				nearz_clip(x0,y0,z0,w0,x1,y1,z1,w1);//,&x0,&y0,&z0,&w0);
				x0=xout;y0=yout;z0=zout;w0=wout;
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

			pvr_prim(&verts[0], sizeof(pvr_vertex_t) * (3 + usespare));

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
			continue;
		} // for i num faces
		} // for j num models
		pvr_list_finish();
		pvr_scene_finish();
	}
	return 0;
}
