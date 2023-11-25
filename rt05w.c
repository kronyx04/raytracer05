#include <stdio.h>
#include <SDL2/SDL.h>
#include <math.h>
#include <unistd.h>

#define CANVAS_WIDTH 960
#define CANVAS_HEIGHT 960
#define inf 9999999999

//SDL do not touch
SDL_Event event;
SDL_Window *window;
SDL_Renderer *renderer;

//structs
struct vector {
	float x;
	float y;
	float z;
};

struct rgb {
	int r;
	int g;
	int b;
};

struct sphere {
	struct vector center;
	float radius;
	struct rgb color;
	float specular; //shininess
	float reflective;
	char null; //0 for visible, 1 for not visible
};

struct light {
	char type;
	float intensity;
	struct vector position; //can also represent the direction of a directional light
};

// constants

const struct rgb red = {255, 0, 0};
const struct rgb green = {0, 255, 0};
const struct rgb blue = {0, 0, 255};
const struct rgb yellow = {255, 255, 0};
const struct rgb Aradia = {160, 1, 3};
const struct rgb Nepeta = {65, 102, 0};
const struct rgb Vriska = {0, 65, 130};
const struct rgb Sollux = {161, 161, 0};
const struct rgb pawl = {0, 147, 205};
const struct rgb black = {0, 0, 0};
const struct rgb white = {255, 255, 255};

const struct rgb BACKGROUND = white;

//the scene
const struct vector O = {0, 0, 0}; //Origin, the camera, rays come from here

const float VP[] = {100, 100, 100}; //the viewport, values are width, height and
									//distance from the camera respectively in cm
	//spheres
const struct sphere sphere1 = {{0, -25, 350}, 100, pawl, 500, 0.2, 1};
const struct sphere sphere2 = {{-50, 50, 325}, 50, pawl, 500, 0.2, 1};
const struct sphere sphere3 = {{50, 50, 325}, 50, pawl, 500, 0.2, 1};
const struct sphere sphere4 = {{-50, 100, 325}, 50, pawl, 500, 0.2, 1};
const struct sphere sphere5 = {{50, 100, 325}, 50, pawl, 500, 0.2, 1};

const struct sphere sphere6 = {{0, -100, 300}, 100, Aradia, 500, 0.2, 0};
const struct sphere sphere7 = {{200, 0, 400}, 100, Vriska, 500, 0.3, 0}; 
const struct sphere sphere8 = {{-200, 0, 400}, 100, Nepeta, 10, 0.4, 0};
const struct sphere sphere9 = {{0, -500100, 0}, 500000, Sollux, 1000, 0.2, 0}; 

const struct sphere spheres[] = {
	sphere1,
	sphere2,
	sphere3,
	sphere4,
	sphere5,
	sphere6,
	sphere7,
	sphere8,
	sphere9,
};
const int SPHERE_COUNT = 9;

	//lights
const struct light light1 = {'a', 0.2}; 
const struct light light2 = {'p', 0.6, {200, 100, 0}};
const struct light light3 = {'d', 0.2, {100, 400, 400}}; //the "position" 
																   //is the 
																   //direction
struct light lights[] = {
	light1,
	light2,
	light3
	};
const int LIGHT_COUNT = 3;

//#####functions#####

//putpixel simplifies the SDL ordeal, its mission is to set the color of a pixel in the canvas using cartesian coordinates
int PutPixel(int x, int y, struct rgb color) {
	int sx, sy;
	sx = CANVAS_WIDTH/2 + x;
	sy = CANVAS_HEIGHT/2 - y;
	SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
	SDL_RenderDrawPoint(renderer, sx, sy);
	return 0;
}

struct vector VectorSub(struct vector va, struct vector vb) {
	struct vector result;

	result.x = va.x - vb.x;
	result.y = va.y - vb.y;
	result.z = va.z - vb.z;

	return result;
};

struct vector VectorAdd(struct vector va, struct vector vb) {
	struct vector result;

	result.x = va.x + vb.x;
	result.y = va.y + vb.y;
	result.z = va.z + vb.z;

	return result;
};

struct vector VectorNumAdd(struct vector va, float b) {
	struct vector result;

	result.x = va.x + b;
	result.y = va.y + b;
	result.z = va.z + b;

	return result;
};

struct vector VectorMul(struct vector va, struct vector vb) {
	struct vector result;

	result.x = va.x * vb.x;
	result.y = va.y * vb.y;
	result.z = va.z * vb.z;

	return result;
};

struct vector VectorNumMul(struct vector va, float b) {
	struct vector result;

	result.x = va.x * b;
	result.y = va.y * b;
	result.z = va.z * b;

	return result;
};

struct rgb RgbNumMul(struct rgb a, float b) {
	struct rgb result;

	result.r = a.r * b;
	result.g = a.g * b;
	result.b = a.b * b;
	
	if (result.r > 255)
		result.r = 255;
	if (result.r < 0)
		result.r = 0;

	if (result.g > 255)
		result.g = 255;
	if (result.g < 0)
		result.g = 0;

	if (result.b > 255)
		result.b = 255;
	if (result.b < 0)
		result.b = 0;

	return result;
};

struct rgb RgbAdd(struct rgb a, struct rgb b) {
	struct rgb result;

	result.r = a.r + b.r;
	result.g = a.g + b.g;
	result.b = a.b + b.b;

	return result;
};

struct vector VectorDiv(struct vector va, struct vector vb) {
	struct vector result;

	result.x = va.x / vb.x;
	result.y = va.y / vb.y;
	result.z = va.z / vb.z;

	return result;
};

struct vector VectorNumDiv(struct vector va, float b) {
	struct vector result;

	result.x = va.x / b;
	result.y = va.y / b;
	result.z = va.z / b;

	return result;
};

struct vector VectorInv(struct vector v) {
	struct vector result;

	result.x = -1 * v.x;
	result.y = -1 * v.y;
	result.z = -1 * v.z;

	return result;
};

float Dot(struct vector va, struct vector vb) {
	return va.x*vb.x + va.y*vb.y + va.z*vb.z;
};

struct vector ReflectRay(struct vector R, struct vector N) {
	return VectorSub(VectorNumMul(VectorNumMul(N,2),Dot(N,R)), R);
};

float Length(struct vector v) {
	return sqrt(Dot(v, v));
};

struct vector CanvasToViewport(int x, int y) { //takes a canvas coordinate and spits out
	struct vector VPP;						   //a VP coordinate

	VPP.x = x * VP[0]/CANVAS_WIDTH;
	VPP.y = y * VP[1]/CANVAS_HEIGHT;
	VPP.z = VP[2];

	return VPP;
};

//ALGEBRA GRAAAAHHHHHHH
void IntersectRaySphere(float *t1, float *t2, struct vector O, struct vector D, struct sphere sphere) {
	float r = sphere.radius;
	struct vector CO = VectorSub(O, sphere.center);

	float a = Dot(D,D);
	float b = 2*Dot(CO, D);
	float c = Dot(CO, CO) - r*r;
	float discriminant = b*b - 4*a*c;

	if (discriminant < 0) {
		*t1 = inf;
		*t2 = inf;
		return;
	};

	*t1 = (-b + sqrt(discriminant)) / (2*a);
	*t2 = (-b - sqrt(discriminant)) / (2*a);
	return;
};

void ClosestIntersection(struct sphere *csphere, float *ct, struct vector O, struct vector D, float t_min, float t_max) {
	float closest_t = inf;
	struct sphere closest_sphere;
	closest_sphere.null = 1;
	float t[] = {0, 0}; 

	for (int ns=0 ; ns<=SPHERE_COUNT-1 ; ns++) {
		IntersectRaySphere(&t[0], &t[1], O, D, spheres[ns]);
		if (t[0] >= t_min && t[0] <= t_max && t[0] < closest_t && !spheres[ns].null) {
			closest_t = t[0];
			closest_sphere = spheres[ns];
		};
		if (t[1] >= t_min && t[1] <= t_max && t[1] < closest_t && !spheres[ns].null) {
			closest_t = t[1];
			closest_sphere = spheres[ns];
		};
	};
	*csphere = closest_sphere;
	*ct = closest_t;
	return;
};

float ComputeLighting(struct vector P, struct vector N, struct vector V, float s) {
	float i = 0;
	float n_dot_l;
	float r_dot_v;
	int nl;
	struct vector L;
	struct vector R;
	float t_min;
	float t_max;
	struct sphere shadow_sphere;
	float shadow_t;

	for (nl = 0; nl <= LIGHT_COUNT-1; nl++) {
		if (lights[nl].type == 'a') 
			i += lights[nl].intensity;
		else {
			if (lights[nl].type == 'p') {
				L = VectorSub(lights[nl].position, P);
				t_max = 1;
			} else {
				L = lights[nl].position;
				t_max = inf;
			};
			//Shadow check
			ClosestIntersection(&shadow_sphere, &shadow_t, P, L, 0.001, t_max);
			if (!shadow_sphere.null) {
				continue;
			};

			//Diffuse
			n_dot_l = Dot(N, L);
			if (n_dot_l > 0)
				i += lights[nl].intensity * n_dot_l/(Length(N) * Length(L));

			//Specular
			if (s != -1) {
				R = ReflectRay(L, N);
				r_dot_v = Dot(R, V);
				if (r_dot_v > 0) {
					i += lights[nl].intensity * pow(r_dot_v/(Length(R) * Length(V)), s);
				};
			};
		};
	};
	return i;
};

//TraceRay, the big one, self explanatory
struct rgb TraceRay(struct vector O, struct vector D, float t_min, float t_max, int recursion_depth) {
	struct sphere closest_sphere;
	float closest_t;
	float r;
	struct rgb local_color;
	struct rgb reflected_color;
	struct vector R;
	struct vector P;
	struct vector N;


	ClosestIntersection(&closest_sphere, &closest_t, O, D, t_min, t_max);

	if (closest_sphere.null) {
		return BACKGROUND;
	};

	//compute local color
	P = VectorAdd(VectorNumMul(D, closest_t), O);
	N = VectorSub(P, closest_sphere.center);	   // compute sphere normal
	N = VectorNumDiv(N, Length(N));							   // at intersection
															   
	local_color = RgbNumMul(closest_sphere.color, ComputeLighting(P, N, VectorInv(D), closest_sphere.specular));

	//if we hit the recursion limit ot the object is not reflective, we're done
	r = closest_sphere.reflective;
	if (recursion_depth <= 0 || r <= 0) {
		return local_color;
	};

	//compute the reflected color
	R = ReflectRay(VectorInv(D), N);
	reflected_color = TraceRay(P, R, 0.1, inf, recursion_depth - 1);
	return RgbAdd(RgbNumMul(local_color, 1-r), RgbNumMul(reflected_color, r));

};

//###### main function ######
int WinMain(int argc, char **argv) {

//window name
	const char* title = "bepis";

//initialize SDL, the window and the renderer
	SDL_Init(SDL_INIT_VIDEO);
	window = SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, CANVAS_WIDTH, CANVAS_HEIGHT, 0);
	renderer = SDL_CreateRenderer(window, -1, 0);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
	SDL_RenderClear(renderer);

//stuff
	int x;
	int y;
	int recursion_depth = 6;
	struct vector D;
	struct rgb color;

	for(x=-CANVAS_WIDTH/2 ; x<=CANVAS_WIDTH/2 ; x++) {
		for(y=CANVAS_HEIGHT/2 ; y>=-CANVAS_HEIGHT/2 ; y--) {
			D = CanvasToViewport(x, y);
			color = TraceRay(O, D, 1, inf, recursion_depth);
			PutPixel(x, y, color);
		};
	};
	
//loop to keep the window open and additional stuff for gracefully exiting
//the render call is also here
	SDL_RenderPresent(renderer);
	while (1) {
		if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
			break;
	};
	SDL_DestroyWindow(window);
	SDL_Quit();
	return 0;
}

