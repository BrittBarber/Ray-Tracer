// Ray Tracer Assignment JSON Schema
// This one works for all the basic features of the assignment
// Winter 2021

#ifndef SCHEMA_H
#define SCHEMA_H

#include <string>
#include <vector>
#include <glm/glm.hpp>

typedef glm::vec3 RGB;
typedef glm::vec3 Vertex;
typedef glm::vec3 Vector;

struct Camera {
  float field;
  RGB background;
  
  Camera() : field(0), background(RGB(0,0,0)) {}
  Camera(float _field, RGB _background) : field(_field), background(_background) {}
};

struct Material {
  RGB ambient;
  RGB diffuse;
  RGB specular;
  float shininess;
  
  RGB reflective;
  RGB transmissive;
  float refraction;

  Material() :
    ambient(RGB(0,0,0)), diffuse(RGB(0,0,0)), specular(RGB(0,0,0)), shininess(1),
    reflective(RGB(0,0,0)), transmissive(RGB(0,0,0)), refraction(0) {}
  Material(RGB _ambient, RGB _diffuse) :
    ambient(_ambient), diffuse(_diffuse), specular(RGB(0,0,0)), shininess(1),
    reflective(RGB(0,0,0)), transmissive(RGB(0,0,0)), refraction(0) {}
  Material(RGB _ambient, RGB _diffuse, RGB _specular, float _shininess) :
    ambient(_ambient), diffuse(_diffuse), specular(_specular), shininess(_shininess),
    reflective(RGB(0,0,0)), transmissive(RGB(0,0,0)), refraction(0) {}
  Material(RGB _ambient, RGB _diffuse, RGB _specular, float _shininess,
    RGB _reflective, RGB _transmissive, float _refraction) :
      ambient(_ambient), diffuse(_diffuse), specular(_specular), shininess(_shininess),
      reflective(_reflective), transmissive(_transmissive), refraction(_refraction) {}
};

struct Object {
  std::string type;
  Material material;
  
  Object(std::string _type, Material _material) :
    type(_type), material(_material) {}
};

struct Sphere : public Object {
  float radius;
  Vertex position;

  Sphere(Material _material, float _radius, Vertex _position) :
    Object("sphere", _material), radius(_radius), position(_position) {}
};

struct Plane : public Object {
  Vertex position;
  Vector normal;
  
  Plane(Material _material, Vertex _position, Vector _normal) :
    Object("plane", _material), position(_position), normal(_normal) {}
};

struct Triangle {
  Vertex vertices[3];
};

struct Mesh : public Object {
  std::vector<Triangle> triangles;

  Mesh(Material _material, std::vector<Triangle> _triangles) :
    Object("mesh", _material), triangles(_triangles) {}
};

struct Light {
  std::string type;
  // for ambient lights, color is ia
  // for all other kinds of lights, color is both id and is
  // you could separate out those two if necessary
  RGB color;
  
  Light(std::string _type, RGB _color) : type(_type), color(_color) {}
};

struct AmbientLight : public Light {
  AmbientLight(RGB _color) : Light( "ambient", _color) {}
};

struct DirectionalLight : public Light {
  Vector direction;
  
  DirectionalLight(RGB _color, Vector _direction) :
    Light( "directional", _color), direction(_direction) {}
};

struct PointLight : public Light {
  Vertex position;

  PointLight(RGB _color, Vertex _position) :
    Light( "point", _color), position(_position) {}
};

struct SpotLight : public Light {
  Vertex position;
  Vector direction;
  // cutoff is half angle of cone, in degrees
  float cutoff;

  SpotLight(RGB _color, Vertex _position, Vector _direction, float _cutoff) :
    Light( "spot", _color), position(_position), direction(_direction), cutoff(_cutoff) {}
};

struct Scene {
  Camera camera;
  std::vector<Object*> objects;
  std::vector<Light*> lights;
};

#endif

/*** The following is c.json hard-coded using this schema ***/

/*
Scene scene_c = {
  // camera
  Camera( 60, RGB( 0, 0, 0.1 ) ),
  // objects
  {
    new Sphere( Material( RGB( 0.2, 0.2, 0.2 ), RGB( 0.56, 0.24, 0.12 ), RGB( 0.4, 0.4, 0.4 ), 120.888832 ),
      2.5, Vertex( 1.5, 1.5, -8.5 )),
    new Sphere( Material( RGB( 0.3, 0.3, 0.3 ), RGB( 1, 0, 0 ), RGB( 1, 1, 1 ), 25.6 ),
      1.0, Vertex( 0, 0, -4 )),
    new Plane( Material( RGB( 0.2, 0.2, 0.2 ), RGB( 0.56, 0.24, 0.12 ) ),
      Vertex( 0, -2, 0 ), Vector( 0, 1, 0 )),
    new Plane( Material( RGB( 0.5, 0.5, 0.7 ), RGB( 0.72, 0.75, 0.84 ) ),
      Vertex( -2.5, -1, 1.5 ), Vector( 1, 0, -0.3 )),
    new Mesh( Material( RGB( 0.15, 0.15, 0.35 ), RGB( 0.7, 0, 0.8 ) ),
      {
        { Vertex( -0.5, 1.0, -2.1  ), Vertex( -1.0, -0.7, -2.1 ), Vertex( -0.5, -0.7, -1.6 ) },
        { Vertex( -0.5, 1.0, -2.1  ), Vertex( -0.5, -0.7, -1.6 ), Vertex( 0.0, -0.7, -2.1  ) },
        { Vertex( -0.5, 1.0, -2.1  ), Vertex( -0.5, -0.7, -2.6 ), Vertex( -1.0, -0.7, -2.1 ) },
        { Vertex( -0.5, 1.0, -2.1  ), Vertex( 0.0, -0.7, -2.1  ), Vertex( -0.5, -0.7, -2.6 ) },
        { Vertex( -0.5, -0.7, -1.6 ), Vertex( -1.0, -0.7, -2.1 ), Vertex( -0.5, -0.7, -2.6 ) },
        { Vertex( -0.5, -0.7, -1.6 ), Vertex( -0.5, -0.7, -2.6 ), Vertex( 0.0, -0.7, -2.1  ) }
      } )
  },
  // lights
  {
    new AmbientLight( RGB( 0.7, 0.7, 0.7 ) ),
    new DirectionalLight( RGB( 0.9, 0.9, 0.9 ), Vector( -0.1, -1, 0 ) ),
    new PointLight( RGB( 0.6, 0.6, 0.0 ), Vertex( -1, 8, 1 ) ),
    new SpotLight( RGB( 0.2, 0.32, 0.85 ), Vertex( -2, 3.9, 1.5 ),
      Vector( 0.4, -0.8, -1 ), 10 )
  }
};
*/