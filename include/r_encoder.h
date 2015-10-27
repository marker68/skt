/*
 * r_encoder.h
 *
 *  Created on: 2014/12/31
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_R_ENCODER_H_
#define INCLUDE_R_ENCODER_H_

#include <iostream>
#include "encoder.h"

namespace PQLearn {

class REncoder : public Encoder {
protected:
	float * rq;
	PQConfig config2;
public:
	REncoder() : Encoder::Encoder() {
		rq = nullptr;
	}
	virtual ~REncoder() {
		::delete rq;
		rq = nullptr;
	}

	void load_codebooks(
			const char * cq_path,
			const char * pq_path,
			const char * rerank,
			bool verbose) {
		load_codebook<float>(cq_path,config,cq,0,verbose);
		load_codebook<float>(pq_path,config,pq,1,verbose);
		load_codebook<float>(rerank,config2,rq,1,verbose);
		size = static_cast<int>(pow(config.kc,config.mc));
		if(verbose)
			cout << "The size of ivf:" << size << endl;
		cout << "--> Settings: (kc,mc,kp,mp)=" << config.kc << " "
				<< config.mc << " " << config.kp << " " << config.mp << endl;
	}

	void load_encoded_data(const char *, bool);
	void output(const char*, const char *, bool);

	template<typename DataType>
	inline void encode(const char *, int, bool);
};

template<typename DataType>
inline void REncoder::encode(
		const char * filename,
		int offset,
		bool verbose) {
	// Load data from file
	DataType * data;
	config.N = load_data<DataType>(filename,data,offset,config.dim,verbose);
	if(config.mc <= 0 || config.mp <= 0) {
		cerr << "Nothing to do" << endl;
		return;
	}

	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	// Encode data

	double d = 0.0, d_tmp = 0.0;
	int bsc = config.dim/config.mc,
			bsp = config.dim/config.mp,
			bsr = config.dim/config2.mp;
	int base = 0;
	size_t p = static_cast<size_t>(config.N) / max_threads;

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
			size_t base_c = start * (config.mp + config2.mp);
			size_t base_p = start;
			size_t i, j, k;
			ushort min_p;

			if(end > config.N || i0 ==  max_threads - 1)
				end = config.N;

//			v_tmp = data + start * static_cast<size_t>(config.dim);
			for(i = start; i < end; i++) {
				v_tmp = data + static_cast<size_t>(config.dim) * pid[base_p++];
				for(j = 0; j < config.dim; j++) {
					v_tmp0[j] = static_cast<float>(v_tmp[j]);
				}
				v_tmp5 = v_tmp0;

				// Find the coarse quantizer code
				base = 0;
				v_tmp3 = cq; // fixed this pointer to point into the begin of cq
				for(j = 0; j < config.mc; j++) {
					v_tmp1 = v_tmp3 + cid[base_u] * bsc;
					// Residual vector
					for(k = 0; k < bsc; k++) {
						v_tmp2[base++] = v_tmp5[k] - v_tmp1[k];
					}
					v_tmp5 += bsc;
					v_tmp3 += config.kc * bsc;
					base_u++;
				}

				// Encoding
				// Point this to the begin of pq
				base = 0;
				v_tmp3 = pq;
				for(j = 0; j < config.mp; j++) {
					v_tmp1 = v_tmp3 + codes[base_c] * bsp;
					// Residual vector
					for(k = 0; k < bsp; k++) {
						v_tmp2[base++] -= v_tmp1[k];
					}
					v_tmp3 += config.kp * bsp;
					base_c++;
				}

				// Encoding
				// Point this to the begin of rq
				v_tmp3 = rq;
				v_tmp4 = v_tmp2;
				for(j = 0; j < config2.mp; j++) {
					d = DBL_MAX;
					// Linear search
					for(k = 0; k < config2.kp; k++) {
						d_tmp = SimpleCluster::distance_l2<float>(
								v_tmp4, v_tmp3,bsr);
						if(d > d_tmp) {
							d = d_tmp;
							min_p = k;
						}
						v_tmp3 += bsr;
					}
					codes[base_c++] = min_p;
					v_tmp4 += bsr;
				}

//				v_tmp += config.dim;
			}
		}
#ifdef _OPENMP
	}
#endif

	// Update the mp
	config.mp += config2.mp;
	if(verbose)
		cout << "Now config.mp=" << config.mp << endl;

	// A minor distribution procedure
}


} /* namespace PQLearn */

#endif /* INCLUDE_R_ENCODER_H_ */
