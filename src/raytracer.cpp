/*
  Brittaney Barber 
  7786450
  Ray Tracer Assignment 
  Winter 2021

*/

#include "raytracer.h"
#include "json2cpp.h"
#include "schema.h"

#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <cmath>
#include <math.h>

using namespace std;

const char *PATH = "scenes/";

double fov = 60;
colour3 background_colour(0, 0, 0);
const float epsilon = 2.2e-3;
float etaAir = 1.0;

json jScene;
Scene scene;

bool lightEffects(const point3& e, const point3& s, colour3& colour, int reflectionDepth, bool passThrough, bool pick);
bool intersection(Material& material, point3 &N, point3 &p, const point3& e, const point3& d, bool pick);
colour3 light(Material& material, const point3 &N, const point3 &p, const point3 e, bool pick);
colour3 reflect(point3 N, point3 p, point3 e, colour3 reflective, int reflectionDepth, bool pick);
bool refract(const point3& Vi, point3 N, const float& etaI, const float& etaR, point3& Vr, bool pick, bool passing);
/****************************************************************************/

// read in json file and convert it to a scene object 
// by default use scene c
void choose_scene(char const *fn) {
	if (fn == NULL) {
		std::cout << "Using default input file " << PATH << "c.json\n";
		fn = "c";
	}

	std::cout << "Loading scene " << fn << std::endl;
	
	std::string fname = PATH + std::string(fn) + ".json";
	std::fstream in(fname);
	if (!in.is_open()) {
		std::cout << "Unable to open scene file " << fname << std::endl;
		exit(EXIT_FAILURE);
	}
	
	in >> jScene;
	
	if (json_to_scene(jScene, scene) < 0) {
		cout << "Error in scene file " << fname << endl;
		exit(EXIT_FAILURE);
	}

	fov = scene.camera.field;
	cout << "Setting fov to " << fov << " degrees.\n";
	
	background_colour = scene.camera.background;
	cout << "Setting background colour to " << glm::to_string(background_colour) << endl;
}


// called in q1.c 
// returns a boolean indicating if there was an intersection and what colour the pixel should be
bool trace(const point3 &e, const point3 &s, colour3 &colour, bool pick) {
	return lightEffects(e, (s-e), colour, pick, 0, false);
}


// recursively finds the refelction, refraction, and transmission of a point 
bool lightEffects(const point3& e, const point3& d, colour3& colour, int reflectionDepth, bool passThrough, bool pick) {
	Material material;
	colour3 reflective = colour3(1.0);
	colour3 transmissive;
	float refraction;
	point3 N;
	point3 p;

	reflectionDepth++;

	if (reflectionDepth < 7) {

		if (intersection(material, N, p, e, d, pick)) {
			colour = light(material, N, p, e, pick);

			reflective = material.reflective;
			colour = colour + reflect(N, p, e, reflective, reflectionDepth, pick);

			transmissive = material.transmissive;
			colour3 nextColour;
			point3 nextRay = p + d;
			bool refractPass = false;
			bool internalReflection = true;

			float etaI = etaAir;
			refraction = material.refraction;

			if (passThrough)
				swap(etaI, refraction);
		
			if ((internalReflection = refract(d, N, etaI, refraction, nextRay, pick, passThrough))) 
				refractPass = !passThrough;

			if (internalReflection) {
				if (!lightEffects(p, nextRay, nextColour, reflectionDepth, refractPass, pick)) 
					nextColour = background_colour;
			}

			else {
				if (passThrough)
					N *= -1.0f;
				nextColour = reflect(N, p, e, reflective, reflectionDepth, pick);
			}

			colour = ((colour3(1.0) - transmissive) * colour) + (transmissive * nextColour);
			return true;
		}
	}
	return false;
}


// return a boolean indicating if our ray point intersects with an object in our scene 
bool intersection(Material &material, point3 &N, point3 &p, const point3& e, const point3& d, bool pick){
	float currT = 100000;  //set to a very high t value so we can set t to the lowest value out of all the objects
	bool hitObject = false;

	// traverse the objects
	for (int i = 0; i < scene.objects.size(); i++) {
		Object* obj = scene.objects[i];

		// every object in the scene will have a "type"
		if (obj->type == "sphere") {
			Sphere* sphere = (Sphere*)(obj);

			// get sphere position and radius and calculate d 
			point3 c = sphere->position;
			float r = sphere->radius;

			//calcualte discriminant
			float disc = pow(dot(d, (e - c)), 2) - dot(d, d) * (dot((e - c), (e - c)) - pow(r, 2));

			//check for an intersection
			if (disc >= 0) {
				//calculate t
				float t1 = (dot(-d, (e - c)) - sqrt(disc)) / dot(d, d);
				float t2 = (dot(-d, (e - c)) + sqrt(disc)) / dot(d, d);
				float t = min(t1, t2);

				if (t < epsilon)
					t = max(t1, t2);

				if (t < currT && t > epsilon) {
					currT = t;
					material = sphere->material;

					//get position of the point 
					p = e + (t * d);

					//get the surface normal 
					N = normalize(p - c);

					hitObject = true;

					if (pick) {
						cout << "You hit a sphere" << endl;
						cout << "First hit spot is at t = " << t << endl;
						cout << "Normal: " << glm::to_string(N) << endl;
					}
				}
			}
		}
		else if (obj->type == "plane") {
			Plane* plane = (Plane*)(obj);
			point3 n = normalize(plane->normal);
			point3 a = plane->position;

			//calcualte denominator 
			float denom = dot(n, d);

			//check for an intersection
			if (denom != 0) {
				//calculate t
				float t = dot(n, (a - e)) / denom;

				//verify hit point is visible to the eye (a positive t value)
				if (t < currT && t >= 0 && t > epsilon) {
					currT = t;
					material = plane->material;

					//get position of the point 
					p = e + (t * d);
					N = n;

					hitObject = true;

					if (pick) {
						cout << "You hit a plane" << endl;
						cout << "First hit spot is at t = " << t << endl;
						cout << "Normal: " << glm::to_string(N) << endl;
					}
				}
			}
		}
		else if (obj->type == "mesh") {
			Mesh* mesh = (Mesh*)(obj);

			//check each triangle 
			for (int i = 0; i < mesh->triangles.size(); i++) {
				point3 a = point3(mesh->triangles[i].vertices[0]);
				point3 b = point3(mesh->triangles[i].vertices[1]);
				point3 c = point3(mesh->triangles[i].vertices[2]);

				point3 n = normalize(cross(b - a, c - b));

				//check if we intercept plane and if so get the hit point
				float denom = dot(n, d);
				if (denom != 0) {
					float t = dot(n, (a - e)) / denom;
					if (t < currT && t >= 0 && t > epsilon) {
						point3 x = e + (t * d);

						//check if all cross product signs are positive 
						if (dot(cross(b - a, x - a), n) > 0 && dot(cross(c - b, x - b), n) > 0 && dot(cross(a - c, x - c), n) > 0) {
							currT = t;
							material = mesh->material;

							//get position of the point 
							p = x;
							N = n;

							hitObject = true;

							if (pick) {
								cout << "You hit a mesh" << endl;
								cout << "First hit spot is at t = " << t << endl;
								cout << "Normal: " << glm::to_string(N) << endl;

							}
						}
					}

				}
			}

		}
	}
	return hitObject;
}


// return the colour of out point based on the lights in our scene 
colour3 light( Material &material, const point3 &N, const point3 &p, const point3 e, bool pick ) {
	//iterate through lights for each object 
	Material temp;
	point3 ambient, lightEffect, diffuse, specular;
	point3 D, P, L, id, is, R;
	point3 V = normalize(e - p);

	//get material properties 
	point3 ka = material.ambient;
	point3 kd = material.diffuse;
	point3 ks = material.specular;
	float alpha = material.shininess;

	for (int i = 0; i < scene.lights.size(); i++) {
		Light* light = scene.lights[i];
		bool shadow = false;

		//get ambience 
		if (light->type == "ambient") {
			AmbientLight* a = (AmbientLight *)(light);
			ambient = a->color * ka;
		}
		
		//directional light 
		else if (light->type == "directional") {
			DirectionalLight* d = (DirectionalLight*)(light);
			D = d->direction;
			L = normalize(-D);

			shadow = intersection(temp, point3(0,0,0), point3(0,0,0), p, L, pick);
		}

		//point light 
		else if (light->type == "point") {
			PointLight* pl = (PointLight*)(light);
			P = pl->position;
			L = normalize(P - p);

			point3 pPrime;
			shadow = intersection(temp, point3(0, 0, 0), pPrime, p, L, pick);
		    shadow = (shadow && (distance(p, L) > distance(p, pPrime)));
		}

		// spot light 
		else if (light->type == "spot") {
			SpotLight* sl = (SpotLight*)(light);
			P = sl->position;
			L = normalize(P - p);
			D = normalize(sl->direction);

			//check if point is in the cutoff 
			float cutoff = sl->cutoff;
			float theta = acos(dot(D, -L));
			if (theta <= cutoff) {
				point3 pPrime;
				shadow = intersection(temp, point3(0, 0, 0), pPrime, p, L, pick);
				shadow = (shadow && (distance(p, L) > distance(p, pPrime)));
			}
		}
		//sum the diffuse and specular components of our scene if its not a shadow 
		if (!shadow) {
			R = normalize(2 * dot(N, L) * N - L);

			id = light->color; 
			is = light->color;

			diffuse = id * kd * max(dot(N, L), 0.0f);
			specular = is * ks * pow( max(dot(R, V), 0.0f), alpha);

			lightEffect += clamp(diffuse, 0.0f, 1.0f) + clamp(specular, 0.0f, 1.0f);
			lightEffect = clamp(lightEffect, 0.0f, 1.0f);

			if (pick) {
				printf_material(material);
				cout << endl;
				cout << "light type: " << light->type << endl;
				cout << "diffuse: " << glm::to_string(diffuse) << endl;
				cout << "specular: " << glm::to_string(specular) << endl;
				cout << "not a shadow" << endl;
				cout << endl;
				cout << endl;
			}
		}	
	} 

	point3 I = ambient + lightEffect;

	return colour3(glm::clamp(I.x, 0.0f, 1.0f), glm::clamp(I.y, 0.0f, 1.0f), glm::clamp(I.z, 0.0f, 1.0f)); 
} 


// calculate reflection and the corresponding colour
colour3 reflect(point3 N, point3 p, point3 e, colour3 reflective, int reflectionDepth, bool pick) {
	point3 V = normalize(e - p);
	point3 R = normalize(2 * dot(N, V) * N - V);
	colour3 reflectColour;

	if (lightEffects(p, R, reflectColour, reflectionDepth, false, pick)) 
		reflectColour = reflectColour * reflective;
	
	return reflectColour;
}


// returns a boolean indicating if a point should be refracting 
bool refract(const point3& Vi, point3 N, const float& etaI, const float& etaR, point3& Vr, bool pick, bool passing) {
	bool refract = false;
	float ViN = dot(normalize(Vi), N);

	if (passing)
		N *= -1.0f;
	else
		ViN *= -1.0f;

	float disc = 1.0 - (pow(etaI, 2) * (1.0 - pow(dot(Vi,N), 2))) / pow(etaR, 2);

	if (disc < 0) {
		Vr = normalize(((etaI * normalize(Vi - (N * ViN))) / etaR) - N * sqrt(disc));
	}
	
	return refract;
}