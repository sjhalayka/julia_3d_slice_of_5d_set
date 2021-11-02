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


quintonion sin(const quintonion& in)
{
	//	float d = in.x * in.x + in.y * in.y + in.z * in.z + in.w * in.w;


	float e =	in.vertex_data[1] * in.vertex_data[1] +
				in.vertex_data[2] * in.vertex_data[2] +
				in.vertex_data[3] * in.vertex_data[3] +
				in.vertex_data[4] * in.vertex_data[4];

	//	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;

	out.vertex_data[0] = sin(in.vertex_data[0]) * cosh(l_e);
	out.vertex_data[1] = in.vertex_data[1] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[2] = in.vertex_data[2] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[3] = in.vertex_data[3] / l_e * cos(in.vertex_data[0]) * sinh(l_e);
	out.vertex_data[4] = in.vertex_data[4] / l_e * cos(in.vertex_data[0]) * sinh(l_e);

	return out;
}

quintonion exp(const quintonion& in)
{
	//	float d = in.x * in.x + in.y * in.y + in.z * in.z + in.w * in.w;
	float d = in.vertex_data[0] * in.vertex_data[0] +
		in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;
	
	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = exp(in.vertex_data[0]) * cos(l_e);
	}

	if (l_e != 0)
	{
		out.vertex_data[1] = in.vertex_data[1] / l_e * exp(in.vertex_data[0]) * sin(l_e);
		out.vertex_data[2] = in.vertex_data[2] / l_e * exp(in.vertex_data[0]) * sin(l_e);
		out.vertex_data[3] = in.vertex_data[3] / l_e * exp(in.vertex_data[0]) * sin(l_e);
		out.vertex_data[4] = in.vertex_data[4] / l_e * exp(in.vertex_data[0]) * sin(l_e);
	}

	return out;
}

quintonion ln(const quintonion& in)
{
	float d = in.vertex_data[0] * in.vertex_data[0] +
		in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float e = in.vertex_data[1] * in.vertex_data[1] +
		in.vertex_data[2] * in.vertex_data[2] +
		in.vertex_data[3] * in.vertex_data[3] +
		in.vertex_data[4] * in.vertex_data[4];

	float l_d = sqrtf(d);
	float l_e = sqrtf(e);

	quintonion out;

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		out.vertex_data[1] = in.vertex_data[1] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[2] = in.vertex_data[2] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[3] = in.vertex_data[3] / l_e * acos(in.vertex_data[0] / l_d);
		out.vertex_data[4] = in.vertex_data[4] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}

quaternion traditional_mul(const quaternion& in_a, const quaternion& in_b)
{
	quaternion out;

	// perform traditional multiply
	out.x = in_a.x * in_b.x - in_a.y * in_b.y - in_a.z * in_b.z - in_a.w * in_b.w;
	out.y = in_a.x * in_b.y + in_a.y * in_b.x + in_a.z * in_b.w - in_a.w * in_b.z;
	out.z = in_a.x * in_b.z - in_a.y * in_b.w + in_a.z * in_b.x + in_a.w * in_b.y;
	out.w = in_a.x * in_b.w + in_a.y * in_b.z - in_a.z * in_b.y + in_a.w * in_b.x;

	return out;
}

quintonion mul(const quintonion& in_a, const quintonion& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	quintonion ln_a = ln(in_a);
	quintonion ln_b = ln(in_b);

	quintonion out;

	out.vertex_data[0] = ln_a.vertex_data[0] + ln_b.vertex_data[0];
	out.vertex_data[1] = ln_a.vertex_data[1] + ln_b.vertex_data[1];
	out.vertex_data[2] = ln_a.vertex_data[2] + ln_b.vertex_data[2];
	out.vertex_data[3] = ln_a.vertex_data[3] + ln_b.vertex_data[3];
	out.vertex_data[4] = ln_a.vertex_data[4] + ln_b.vertex_data[4];

	return exp(out);
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

	float d = 0;
	float e = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		d += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		e += (in.vertex_data[i] * in.vertex_data[i]);

	if (d == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(d);
	const float l_e = sqrtf(e);
	const float d_b2 = powf(d, beta / 2.0f);

	const float theta = beta * acos(in.vertex_data[0] / l_d);

	out.vertex_data[0] = d_b2 * cos(theta);

	if (e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = (in.vertex_data[i]/l_e) * d_b2 * sin(theta);
	}

	return out;
}




inline float iterate(
	quintonion Z,
	quintonion C,
	const short unsigned int max_iterations,
	const float threshold)
{
	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		quintonion Z_orig = Z;

		quintonion Z_base = Z;
		Z = mul(Z, Z_base);
		Z = mul(Z, Z_base);
		Z = mul(Z, Z_base);
		Z = Z + C;

		// Z = pow_number_type(Z_orig, 4.0) + C;

	//	Z = sin(Z) + mul(sin(Z), C);


		//quaternion qc;
		//qc.x = C.vertex_data[0];
		//qc.y = C.vertex_data[1];
		//qc.z = C.vertex_data[2];
		//qc.w = C.vertex_data[3];

		//quintonion s = sin(Z);

		//quaternion qs;
		//qs.x = s.vertex_data[0];
		//qs.y = s.vertex_data[1];
		//qs.z = s.vertex_data[2];
		//qs.w = s.vertex_data[3];

		//quaternion f = traditional_mul(qs, qc);

		//f.x += qs.x;
		//f.y += qs.y;
		//f.z += qs.z;
		//f.w += qs.w;

		//Z.vertex_data[0] = f.x;
		//Z.vertex_data[1] = f.y;
		//Z.vertex_data[2] = f.z;
		//Z.vertex_data[3] = f.w;
		//Z.vertex_data[4] = 0;


		if (Z.magnitude() >= threshold)
			break;
	}

	return Z.magnitude();
}


