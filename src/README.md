# Ray Tracer Assignment files README

Do not hand in this README! Make your own.

## Included files:

* `common.h` and `main.cpp` are as before, except modified to handle command-line arguments
* `q1.cpp` does the OpenGL rendering. You probably don't want to modify this, unless you are improving the ray tracer to the project. It handles casting rays into the view plane, and draws the output from your ray tracer.
* `raytracer.h` has declarations for the two functions that do the actual work, as described in the assignment: `choose_scene` and `trace`.
* `raytracer.cpp` provides a default implementation for the functions from `raytracer.h`. They are not complete, but they give you some output.
* `json.hpp` is a third-party JSON parser for C++.

## Implementation notes

The `choose_scene()` function given in `raytracer.cpp` opens up a JSON file that corresponds to the selected scene; e.g. if you run it as `q1 c` then it will load and parse a scene from the JSON file `scenes/c.json`. It will set the camera parameters (the FOV angle and background colour), which are used by the given `q1.cpp` program. If you replace this function, you should still set those values (found as extern declarations in `raytracer.h`).

The given `trace()` function is (of course) incomplete. It demonstrates how to access the parsed JSON data to get data describing the scene.

You can continue to use the parsed JSON data directly from the `json scene` object. However, it will be slow. Eventually it will be a bottleneck, especially if you are going to develop your work into a project. Here are some options:

* You can copy the data from the JSON file into your own data structures for the scene. These data structures are entirely your choice. Especially if you build the raytracer project, you will need to be very familiar with these data structures because you will be re-organizing them to build acceleration data structures.
* You can use the schema that I have provided. The details are described elsewhere in this README. You can modify, rearrange, alter, etc. these data structures as much as you like.

To get the data into your data structures, whichever you choose, you can either read the JSON files and store them into your data structures at run time, or you can pre-process the JSON data (automatically or manually) and hard-code the initialization of your data structures. There is code to do both of these for the schema I have provided. You are free to use, modify, etc. this code however you like, for your chosen data structures.

## Utils directory

This directory contains some code you may find useful. All of it is provided without support, and it is not thoroughly tested.

* `schema.h` is a C++ struct-based implementation of the schema used by the JSON scene files. You may find it helps you understand the JSON files better. Or you may use it directly. Included at the bottom of the file is a commented-out hard-coded version of `scenes/c.json` represented using the schema. You can take this code and move it into a `.cpp` file and use it directly for scene data.
* `json2cpp.cpp` is a command-line utility that will convert a JSON file to the C++ schema format, like the hard-coded scene in `schema.h`. It contains `json_to_scene()`, which will read in a JSON file and store it using the schema. It also contains `scene_to_cpp()`, which will take an in-memory schema structure and write it out as C++ source code. You may find the first function especially useful to read in JSON scenes at run-time. Again, you may re-use or modify these anywhere you need in your assignment.
* `obj2json.py` is a Python script that will convert a [Wavefront OBJ file](https://en.wikipedia.org/wiki/Wavefront_.obj_file) to a triangle mesh in our JSON format. You may find this a convenient way to get large models into your ray tracer, especially if you are extending it into a project. It can transform the mesh as it reads it.

Note that you should be citing your sources for OBJ files or other geometries, and any tools you use to generate them.
