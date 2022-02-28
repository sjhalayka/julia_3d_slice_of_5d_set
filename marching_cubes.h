// Modified from Paul Bourke, Polygonising a Scalar Field
// Source code by Shawn Halayka
// Source code is in the public domain


#ifndef MARCHING_CUBES_H
#define MARCHING_CUBES_H


#include "primitives.h"


#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;

#include <iomanip>
using std::ios_base;

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <map>
using std::map;
using std::multimap;

#include <cmath>
using std::sin;
using std::exp;

#include <utility> 
using std::pair;



namespace marching_cubes
{
	class grid_cube
	{
	public:
		vertex_3 vertex[8];
		float value[8];
	};


	short unsigned int tesselate_grid_cube(vector<trajectory>& trajectory_data, quintonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, const grid_cube& grid, triangle* const triangles);
	void tesselate_adjacent_xy_plane_pair(vector<trajectory>& trajectory_data, quintonion C, float z_w, const float isovalue, const float upper_threshold, const float lower_threshold, short unsigned int max_iterations, vector<float> xyplane0, vector<float> xyplane1, const size_t z, vector<triangle> &triangles, const float x_grid_min, const float x_grid_max, const size_t x_res, const float y_grid_min, const float y_grid_max, const size_t y_res, const float z_grid_min, const float z_grid_max, const size_t z_res);


	quintonion sin(const quintonion& in);

	quintonion exp(const quintonion& in);

	octonion exp(const octonion& in);

	quintonion ln(const quintonion& in);

	
	octonion ln(const octonion& in);
	
	octonion trad_mul_5D(octonion qA, octonion qB);


	quaternion traditional_mul(const quaternion& in_a, const quaternion& in_b);

	quintonion mul(const quintonion& in_a, const quintonion& in_b);
	quintonion div(const quintonion& in_a, const quintonion& in_b);
	octonion div(const octonion& in_a, const octonion& in_b);

	quintonion conj_number_type(quintonion& in);

	quintonion pow_number_type(quintonion& in, float exponent);

	quintonion pow_number_type(quintonion& in, quintonion& exponent);

	octonion get_commutator(octonion in_a, octonion in_b);

	inline float iterate(
		vector<trajectory> &trajectory_data,
		quintonion Z,
		quintonion C,
		float z_w,
		const short unsigned int max_iterations,
		const float threshold)
	{
		Z.vertex_data[3] = z_w;
		Z.vertex_data[4] = z_w;

		trajectory t;
//		const quintonion Z_0 = Z;
		
		//vertex_3 v;
		//v.x = Z_0.vertex_data[0];
		//v.y = Z_0.vertex_data[1];
		//v.z = Z_0.vertex_data[2];

		//t.traj_data.push_back(v);

		for (short unsigned int i = 0; i < max_iterations; i++)
		{
			quintonion Z_base = Z;







			//quintonion mul_left_q = C; // Z;
			//quintonion mul_right_q = sin(Z_base); // Z_base;

			//float s = 0;
			//
			//if (mul_right_q.magnitude() < mul_left_q.magnitude())
			//	s = mul_right_q.magnitude();
			//else
			//	s = mul_left_q.magnitude();

			//float s2 = 2.0f / (s);

			//mul_left_q = mul_left_q * s2;
			//mul_right_q = mul_right_q * s2;

			//float scale = log(mul_left_q.magnitude()) / log(mul_right_q.magnitude());

			//mul_left_q = mul_left_q * (1/s2);
			//mul_right_q = mul_right_q * (1/s2);

			////float rounded = floor(scale + 0.5f);

			////if (rounded == 0)
			////	rounded = 1;

			////float f = fabs(rounded - scale);

			//quintonion mul_quint = mul(mul_left_q, mul_right_q);

			//octonion mul_oct;
			//mul_oct.vertex_data[0] = mul_quint.vertex_data[0];
			//mul_oct.vertex_data[1] = mul_quint.vertex_data[1];
			//mul_oct.vertex_data[2] = mul_quint.vertex_data[2];
			//mul_oct.vertex_data[3] = mul_quint.vertex_data[3];
			//mul_oct.vertex_data[4] = mul_quint.vertex_data[4];

			//octonion mul_left;
			//mul_left.vertex_data[0] = mul_left_q.vertex_data[0];
			//mul_left.vertex_data[1] = mul_left_q.vertex_data[1];
			//mul_left.vertex_data[2] = mul_left_q.vertex_data[2];
			//mul_left.vertex_data[3] = mul_left_q.vertex_data[3];
			//mul_left.vertex_data[4] = mul_left_q.vertex_data[4];

			//octonion mul_right;
			//mul_right.vertex_data[0] = mul_right_q.vertex_data[0];
			//mul_right.vertex_data[1] = mul_right_q.vertex_data[1];
			//mul_right.vertex_data[2] = mul_right_q.vertex_data[2];
			//mul_right.vertex_data[3] = mul_right_q.vertex_data[3];
			//mul_right.vertex_data[4] = mul_right_q.vertex_data[4];

			//octonion mul_oct_trad = trad_mul_5D(mul_left, mul_right);

			//octonion commutator = get_commutator(mul_left, mul_right);
			//octonion commutator2 = get_commutator(mul_oct, mul_oct_trad);


			//Z = sin(Z_base) + mul(sin(Z_base), C);
			






			quintonion sinz = sin(Z_base);

			octonion sinz_oct;
			sinz_oct.vertex_data[0] = sinz.vertex_data[0];
			sinz_oct.vertex_data[1] = sinz.vertex_data[1];
			sinz_oct.vertex_data[2] = sinz.vertex_data[2];
			sinz_oct.vertex_data[3] = sinz.vertex_data[3];
			sinz_oct.vertex_data[4] = sinz.vertex_data[4];

			octonion C_oct;
			C_oct.vertex_data[0] = C.vertex_data[0];
			C_oct.vertex_data[1] = C.vertex_data[1];
			C_oct.vertex_data[2] = C.vertex_data[2];
			C_oct.vertex_data[3] = C.vertex_data[3];
			C_oct.vertex_data[4] = C.vertex_data[4];

			quintonion new_mul = mul(C, sinz);
			//quintonion new_mul = mul(Z_base, Z_base);


			octonion new_mul_oct;
			new_mul_oct.vertex_data[0] = new_mul.vertex_data[0];
			new_mul_oct.vertex_data[1] = new_mul.vertex_data[1];
			new_mul_oct.vertex_data[2] = new_mul.vertex_data[2];
			new_mul_oct.vertex_data[3] = new_mul.vertex_data[3];
			new_mul_oct.vertex_data[4] = new_mul.vertex_data[4];

			//octonion Z_base_oct;
			//Z_base_oct.vertex_data[0] = Z_base.vertex_data[0];
			//Z_base_oct.vertex_data[1] = Z_base.vertex_data[1];
			//Z_base_oct.vertex_data[2] = Z_base.vertex_data[2];
			//Z_base_oct.vertex_data[3] = Z_base.vertex_data[3];
			//Z_base_oct.vertex_data[4] = Z_base.vertex_data[4];

			octonion trad_mul = trad_mul_5D(C_oct, sinz_oct);
			//octonion trad_mul = trad_mul_5D(Z_base_oct, Z_base_oct);

			octonion commutator = get_commutator(C_oct, sinz_oct);
			octonion commutator2 = get_commutator(trad_mul, new_mul_oct);

			//octonion new_mul_oct_normalized = new_mul_oct;
			//float m = new_mul_oct_normalized.magnitude();
			//new_mul_oct_normalized = new_mul_oct_normalized * (1.0f / m);

			//octonion trad_mul_normalized = trad_mul;
			//m = trad_mul_normalized.magnitude();
			//trad_mul_normalized = trad_mul_normalized * (1.0f / m);

			//float dot = 0;

			//for (size_t i = 0; i < 8; i++)
			//	dot += trad_mul_normalized.vertex_data[i] * new_mul_oct_normalized.vertex_data[i];

	/*		octonion diff;
			diff.vertex_data[0] = fabs(trad_mul.vertex_data[0] - new_mul_oct.vertex_data[0]);
			diff.vertex_data[1] = fabs(trad_mul.vertex_data[1] - new_mul_oct.vertex_data[1]);
			diff.vertex_data[2] = fabs(trad_mul.vertex_data[2] - new_mul_oct.vertex_data[2]);
			diff.vertex_data[3] = fabs(trad_mul.vertex_data[3] - new_mul_oct.vertex_data[3]);
			diff.vertex_data[4] = fabs(trad_mul.vertex_data[4] - new_mul_oct.vertex_data[4]);

			cout << diff.vertex_data[1] / commutator.vertex_data[1] << endl;*/

			quintonion z2 = mul(mul(Z_base, Z_base), Z_base);

			octonion z2_oct;
			z2_oct.vertex_data[0] = z2.vertex_data[0];
			z2_oct.vertex_data[1] = z2.vertex_data[1];
			z2_oct.vertex_data[2] = z2.vertex_data[2];
			z2_oct.vertex_data[3] = z2.vertex_data[3];
			z2_oct.vertex_data[4] = z2.vertex_data[4];

			octonion Z_base_oct;
			Z_base_oct.vertex_data[0] = Z_base.vertex_data[0];
			Z_base_oct.vertex_data[1] = Z_base.vertex_data[1];
			Z_base_oct.vertex_data[2] = Z_base.vertex_data[2];
			Z_base_oct.vertex_data[3] = Z_base.vertex_data[3];
			Z_base_oct.vertex_data[4] = Z_base.vertex_data[4];

			//cout << get_commutator(z2_oct, Z_base_oct).magnitude() << endl;


			

			Z = sinz + mul(C, sinz);



			//Z = mul(Z_base, Z_base) + C;




			octonion q = div(trad_mul, new_mul_oct);

			//quintonion q_;
			//q_.vertex_data[0] = q.vertex_data[0];
			//q_.vertex_data[1] = q.vertex_data[1];
			//q_.vertex_data[2] = q.vertex_data[2];
			//q_.vertex_data[3] = q.vertex_data[3];
			//q_.vertex_data[4] = q.vertex_data[4];





			//cout << q.vertex_data[0] << " " << q.vertex_data[1] << " " << q.vertex_data[2] << " " << q.vertex_data[3] << " " << q.vertex_data[4] << endl;

			//float b_norm = q.vertex_data[0] * q.vertex_data[0] +
			//	q.vertex_data[1] * q.vertex_data[1] +
			//	q.vertex_data[2] * q.vertex_data[2] +
			//	q.vertex_data[3] * q.vertex_data[3] +
			//	q.vertex_data[4] * q.vertex_data[4];

			//quintonion q_inv;
			//q_inv.vertex_data[0] = q.vertex_data[0];// / b_norm;
			//q_inv.vertex_data[1] = -q.vertex_data[1];// / b_norm;
			//q_inv.vertex_data[2] = -q.vertex_data[2];// / b_norm;
			//q_inv.vertex_data[3] = -q.vertex_data[3];// / b_norm;
			//q_inv.vertex_data[4] = -q.vertex_data[4];// / b_norm;


			//for (size_t j = 0; j < 5; j++)
			//	cout << q.vertex_data[j] << ' ';

			//cout << endl;

			//for (size_t j = 0; j < 5; j++)
			//	cout << q_inv.vertex_data[j] << ' ';

			//cout << endl;

//			quintonion qp = mul(q_, new_mul);

			//for (size_t j = 0; j < 5; j++)
			//	cout << qp.vertex_data[j] << ' ';

//			cout << endl;

			//for (size_t j = 0; j < 5; j++)
			//	cout << new_mul.vertex_data[j] << ' ';

			//cout << endl;	

			//for (size_t j = 0; j < 5; j++)
			//	cout << trad_mul.vertex_data[j] << ' ';

			//cout << endl;

			//cout << q.magnitude() << " " << trad_mul.magnitude() << endl;

			//cout << commutator.magnitude() << endl;
			//cout << commutator2.magnitude() << endl;


			//cout << endl;


			//octonion m_oct;
			//m_oct.vertex_data[0] = m.vertex_data[0];
			//m_oct.vertex_data[1] = m.vertex_data[1];
			//m_oct.vertex_data[2] = m.vertex_data[2];
			//m_oct.vertex_data[3] = m.vertex_data[3];
			//m_oct.vertex_data[4] = m.vertex_data[4];

			//octonion x = div(trad_mul, m_oct);

			//cout << x.vertex_data[0] << " " << trad_mul.vertex_data[0] << endl;
			////cout << endl << endl;


			//octonion new_mul_oct_normalized = new_mul_oct;
			//float m = new_mul_oct_normalized.magnitude();
			//new_mul_oct_normalized = new_mul_oct_normalized * (1.0f/m);

			//octonion trad_mul_normalized = trad_mul;
			//m = trad_mul_normalized.magnitude();
			//trad_mul_normalized = trad_mul_normalized * (1.0f / m);


			//float b_norm = trad_mul_normalized.vertex_data[0] * trad_mul_normalized.vertex_data[0] +
			//	trad_mul_normalized.vertex_data[1] * trad_mul_normalized.vertex_data[1] +
			//	trad_mul_normalized.vertex_data[2] * trad_mul_normalized.vertex_data[2] +
			//	trad_mul_normalized.vertex_data[3] * trad_mul_normalized.vertex_data[3] +
			//	trad_mul_normalized.vertex_data[4] * trad_mul_normalized.vertex_data[4];

			//octonion temp_b;
			//temp_b.vertex_data[0] = trad_mul_normalized.vertex_data[0] / b_norm;
			//temp_b.vertex_data[1] = -trad_mul_normalized.vertex_data[1] / b_norm;
			//temp_b.vertex_data[2] = -trad_mul_normalized.vertex_data[2] / b_norm;
			//temp_b.vertex_data[3] = -trad_mul_normalized.vertex_data[3] / b_norm;
			//temp_b.vertex_data[4] = -trad_mul_normalized.vertex_data[4] / b_norm;

			//octonion prod = trad_mul_5D(new_mul_oct_normalized, temp_b);

			//cout << prod.vertex_data[0] << endl;
			//cout << endl;





			//float dot = 0;
			//float vec_dot = 0;

			//for (size_t i = 0; i < 8; i++)
			//	dot += trad_mul_normalized.vertex_data[i] * new_mul_oct_normalized.vertex_data[i];
			//
			//float u = 1.0f;

			//float theta = acosf(dot);
			//float sinTheta = sinf(theta);


			//cout << theta << endl;
			//cout << sinTheta << endl;
			//cout << sinTheta/theta << endl;

			//cout << endl;



			if (Z.magnitude() >= threshold)
				break;
			else
			{
				vertex_3 v;
				v.y = new_mul_oct.vertex_data[0];// magnitude();// vertex_data[4]; //magnitude();//static_cast<float>(i + 1);
				v.x = trad_mul.vertex_data[0];// magnitude();// vertex_data[4];// magnitude();// vertex_data[0];
				v.z = 0.0;

				t.traj_data.push_back(v);

				//vertex_3 v;
				//v.x = commutator2.vertex_data[1] + Z_0.vertex_data[0];
				//v.y = commutator2.vertex_data[2] + Z_0.vertex_data[1];
				//v.z = commutator2.vertex_data[3] + Z_0.vertex_data[2];

				//t.traj_data.push_back(v);
			}
		}
	
		trajectory_data.push_back(t);

		return Z.magnitude();
	}

	vertex_3 vertex_interp_refine(
		vector<trajectory>& trajectory_data,
		quintonion C,
		float z_w,
		float isovalue,
		float upper_threshold,
		float lower_threshold,
		short unsigned int max_iterations,
		vertex_3 v0, vertex_3 v1,
		float val_v0, float val_v1);

};

#endif
