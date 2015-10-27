/*
 *  Copyright (C) 2015 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 *  hk_encoder.h
 *
 *  Created on: 2015/04/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef SRC_HK_ENCODER_H_
#define SRC_HK_ENCODER_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include "pq_utilities.h"
#include "encoder.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace PQLearn {

class TreeEncoder : public Encoder {
protected:
	float * tree_cb;
	int * hlb, * llb;
	unsigned char * tree_codes;
	vector<TreeBucket> tree_ivf;
public:
	TreeEncoder();
	virtual ~TreeEncoder();

	void load_codebooks(const char *, const char *, bool);
	template<typename DataType>
	inline void tree_encode(const char *, int, bool);
	template<typename DataType>
	inline void tree_rencode(const char *, int, bool);

	void output(const char *, bool);
	void distribution(bool);
};

template<typename DataType>
inline void TreeEncoder::tree_encode(
		const char * filename,
		int offset,
		bool verbose) {
	// Load data from file
	DataType * data;
	config.N = load_data<DataType>(filename,data,offset,config.dim,verbose);
	size_t N = static_cast<size_t>(config.N);
	if(config.mc <= 0 || config.mp <= 0
			|| config.L <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}

	SimpleCluster::init_array(hlb, N * config.mc);
	SimpleCluster::init_array(llb, N * config.mc);
	SimpleCluster::init_array(tree_codes, N * config.mp);

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data

	double d = 0.0, d_tmp = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	int base = 0;
	size_t p = N / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(base,d,d_tmp)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp1,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp2,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * config.mc;
			size_t base_u2 = start * config.mc;
			size_t base_c = start * config.mp;
			size_t base_p = start;
			size_t i, j, k;
			ushort min_p, min_p2;

			if(end > N || i0 ==  max_threads - 1)
				end = N;

			v_tmp = data + start * static_cast<size_t>(config.dim);
			for(i = start; i < end; i++) {
				for(j = 0; j < config.dim; j++) {
					v_tmp0[j] = static_cast<float>(v_tmp[j]);
				}
				v_tmp5 = v_tmp0;

				// Find the coarse quantizer code
				base = 0;
				v_tmp3 = tree_cb; // fixed this pointer to point into the begin of tree_cb
				d = DBL_MAX;
				// Linear search
				v_tmp4 = v_tmp3;
				for(k = 0; k < config.kc; k++) {
					d_tmp = SimpleCluster::distance_l2<float>(
							v_tmp5, v_tmp4,bsc);
					if(d > d_tmp) {
						d = d_tmp;
						min_p = k;
					}
					v_tmp4 += bsc;
				}
				v_tmp1 = v_tmp3 + min_p * bsc;
				hlb[base_u] = min_p;
				// Residual vector
				for(k = 0; k < bsc; k++) {
					v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
				}
				v_tmp3 = v_tmp4;

				// Find the low level code
				v_tmp5 = v_tmp4 + min_p * config.L * bsc;
				d = DBL_MAX;
				for(k = 0; k < config.L; k++) {
					d_tmp = SimpleCluster::distance_l2<float>(
							v_tmp5, v_tmp2,bsc);
					if(d > d_tmp) {
						d = d_tmp;
						min_p2 = k;
					}
					v_tmp5 += bsc;
				}

				llb[base_u2] = min_p2;
				v_tmp1 = v_tmp4 + min_p * config.L * bsc + min_p2 * bsc;
				base = 0;
				// Residual vector
				for(k = 0; k < bsc; k++) {
					v_tmp2[base++] -= v_tmp1[k];
				}

				base_u++;
				base_u2++;

				// Encoding
				// Point this to the begin of pq
				v_tmp3 = pq;
				v_tmp4 = v_tmp2;
				for(j = 0; j < config.mp; j++) {
					d = DBL_MAX;
					// Linear search
					for(k = 0; k < config.kp; k++) {
						d_tmp = SimpleCluster::distance_l2<float>(
								v_tmp4, v_tmp3,bsp);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp3 += bsp;
					}
					tree_codes[base_c++] = min_p;
					v_tmp4 += bsp;
				}

				v_tmp += config.dim;
			}
		}
#ifdef _OPENMP
	}
#endif
}

template<typename DataType>
inline void TreeEncoder::tree_rencode(
		const char * filename,
		int offset,
		bool verbose) {
	// Load data from file
	DataType * data;
	config.N = load_data<DataType>(filename,data,offset,config.dim,verbose);
	size_t N = static_cast<size_t>(config.N);
	if(config.mc <= 0 || config.mp <= 0
			|| config.L <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}

	if(verbose) {
		cout << "This one has " << config.kc << " high-level cells" << endl;
	}

	SimpleCluster::init_array(hlb, N * config.mc);
	SimpleCluster::init_array(llb, N * config.mc);
	SimpleCluster::init_array(tree_codes, N * config.mp);

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif

	double d = 0.0, d_tmp = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp;
	int base = 0;
	size_t p = N / max_threads;
	size_t count = 0, id;
	int range;

	for(int i = 0; i < config.kc; i++) {
		range = L[i];
		for(int j = 0; j < range; j++) {
			id = pid[count++];
			hlb[id] = i;
		}
	}

	if(count != N) {
		cout << "Ng" << endl;
	}
#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(base,d,d_tmp)
#endif
		for(int i0 = 0; i0 < max_threads; i0++) {
			DataType * v_tmp;
			float * v_tmp0, * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
			if(!SimpleCluster::init_array<float>(v_tmp0,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp1,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}
			if(!SimpleCluster::init_array<float>(v_tmp2,config.dim)) {
				cerr << "Cannot allocate memory" << endl;
				exit(1);
			}

			// Range definition
			size_t start = p * static_cast<size_t>(i0);
			size_t  end = start + p;
			size_t base_u = start * config.mc;
			size_t base_c = start * config.mp;
			size_t base_p = start;
			size_t i, j, k;
			ushort min_p, min_p2;

			if(end > N || i0 ==  max_threads - 1)
				end = N;

			v_tmp = data + start * static_cast<size_t>(config.dim);
			v_tmp3 = tree_cb + config.kc * config.dim;
			for(i = start; i < end; i++) {
				for(j = 0; j < config.dim; j++) {
					v_tmp0[j] = static_cast<float>(v_tmp[j]);
				}
				v_tmp5 = v_tmp0;

				// Find the coarse quantizer code
				base = 0;
				min_p = hlb[i];
				v_tmp1 = tree_cb + min_p * bsc;
				// Residual vector
				for(k = 0; k < bsc; k++) {
					v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
				}

				// Find the low level code
				v_tmp5 = v_tmp3 + min_p * config.L * bsc;
				d = DBL_MAX;
				for(k = 0; k < config.L; k++) {
					d_tmp = SimpleCluster::distance_l2<float>(
							v_tmp5, v_tmp2,bsc);
					if(d > d_tmp) {
						d = d_tmp;
						min_p2 = k;
					}
					v_tmp5 += bsc;
				}

				llb[base_u] = min_p2;
				v_tmp1 = v_tmp3 + min_p * config.L * bsc + min_p2 * bsc;
				base = 0;
				// Residual vector
				for(k = 0; k < bsc; k++) {
					v_tmp2[base++] -= v_tmp1[k];
				}

				base_u++;

				// Encoding
				// Point this to the begin of pq
				v_tmp4 = pq;
				v_tmp5 = v_tmp2;
				for(j = 0; j < config.mp; j++) {
					d = DBL_MAX;
					// Linear search
					for(k = 0; k < config.kp; k++) {
						d_tmp = SimpleCluster::distance_l2<float>(
								v_tmp5, v_tmp4,bsp);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp4 += bsp;
					}
					tree_codes[base_c++] = min_p;
					v_tmp5 += bsp;
				}

				v_tmp += config.dim;
			}
		}
#ifdef _OPENMP
	}
#endif
}

} /* namespace PQLearn */

#endif /* SRC_HK_ENCODER_H_ */
