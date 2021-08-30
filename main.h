#pragma once

#include <vector>
#include <iostream>
using namespace std;

#include "marching_cubes.h"
using namespace marching_cubes;



bool write_triangles_to_binary_stereo_lithography_file(const vector<triangle> &triangles, const char *const file_name)
{
	cout << "Triangle count: " << triangles.size() << endl;

	if (0 == triangles.size())
		return false;

	// Write to file.
	ofstream out(file_name, ios_base::binary);

	if (out.fail())
		return false;

	const size_t header_size = 80;
	vector<char> buffer(header_size, 0);
	const unsigned int num_triangles = static_cast<unsigned int>(triangles.size()); // Must be 4-byte unsigned int.
	vertex_3 normal;

	// Write blank header.
	out.write(reinterpret_cast<const char *>(&(buffer[0])), header_size);

	// Write number of triangles.
	out.write(reinterpret_cast<const char *>(&num_triangles), sizeof(unsigned int));

	// Copy everything to a single buffer.
	// We do this here because calling ofstream::write() only once PER MESH is going to 
	// send the data to disk faster than if we were to instead call ofstream::write()
	// thirteen times PER TRIANGLE.
	// Of course, the trade-off is that we are using 2x the RAM than what's absolutely required,
	// but the trade-off is often very much worth it (especially so for meshes with millions of triangles).
	cout << "Generating normal/vertex/attribute buffer" << endl;

	// Enough bytes for twelve 4-byte floats plus one 2-byte integer, per triangle.
	const size_t data_size = (12 * sizeof(float) + sizeof(short unsigned int)) * num_triangles;
	buffer.resize(data_size, 0);

	// Use a pointer to assist with the copying.
	// Should probably use std::copy() instead, but memcpy() does the trick, so whatever...
	char *cp = &buffer[0];

	for (vector<triangle>::const_iterator i = triangles.begin(); i != triangles.end(); i++)
	{
		// Get face normal.
		vertex_3 v0 = i->vertex[1] - i->vertex[0];
		vertex_3 v1 = i->vertex[2] - i->vertex[0];
		normal = v0.cross(v1);
		normal.normalize();

		memcpy(cp, &normal.x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.z, sizeof(float)); cp += sizeof(float);

		memcpy(cp, &i->vertex[0].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[0].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[0].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[1].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[2].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[2].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &i->vertex[2].z, sizeof(float)); cp += sizeof(float);

		cp += sizeof(short unsigned int);
	}

	cout << "Writing " << data_size / 1048576 << " MB of data to binary Stereo Lithography file: " << file_name << endl;

	out.write(reinterpret_cast<const char *>(&buffer[0]), data_size);
	out.close();

	return true;
}






quintonion conj_number_type(quintonion& in)
{
	quintonion out;

	out.vertex_data[0] = in.vertex_data[0];

	for (size_t i = 1; i < in.vertex_length; i++)
		out.vertex_data[i] = -in.vertex_data[i];

	return out;
}

quintonion pow_number_type(quintonion& in, float exponent)
{
	const float beta = exponent;

	float all_self_dot = 0;
	float imag_self_dot = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float all_len = sqrtf(all_self_dot);
	const float imag_len = sqrtf(imag_self_dot);
	const float self_dot_beta = powf(all_self_dot, beta / 2.0f);

	out.vertex_data[0] = self_dot_beta * std::cos(beta * std::acos(in.vertex_data[0] / all_len));

	if (imag_len != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] * self_dot_beta * sin(beta * acos(in.vertex_data[0] / all_len)) / imag_len;
	}

	return out;
}




inline float iterate(
	quintonion Z,
	quintonion C,
	const short unsigned int max_iterations,
	const float threshold)
{
	// Uncomment these lines for the Mandelbrot set
	//C = Z;

	//Z.vertex_data[0] = 0.0f;
	//Z.vertex_data[1] = 0.0f;
	//Z.vertex_data[2] = 0.0f;
	//Z.vertex_data[3] = 0.0f;
	//Z.vertex_data[4] = 0.0f;


	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = pow_number_type(Z, 2.0) + C;

		if (Z.magnitude() >= threshold)
			break;
	}

	return Z.magnitude();
}


