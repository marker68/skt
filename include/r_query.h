/*
 * r_query.h
 *
 *  Created on: 2015/01/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef R_QUERY_H_
#define R_QUERY_H_

#include <iostream>
#include "pq_query.h"
#include "pq_utilities.h"
#include "pq_algorithm.h"

using namespace std;

namespace PQLearn
{

class PQRQuery : public PQQuery
{
protected:
	float * rq;
	float * norm_r1, * norm_r2;
	float * diff_qr1, * diff_qr2, * dot_cr1, * dot_cr2, * dot_rr;
	PQConfig config2;
public:
	PQRQuery();
	virtual ~PQRQuery();
	void load_codebooks(const char *, const char *,
			const char *, bool);
	int load_encoded_data(const char *, bool);
	inline void pre_compute1();
	inline void pre_compute2(float *);
	inline void search_ivfadc(
			float *,
			float *&, float *&,
			int *&, int *&, int *&,
			int&, int, int,int, bool, bool);
};

/**
 * Precompute all things that are query-independent
 */
inline void PQRQuery::pre_compute1() {
	size_t c = static_cast<size_t>(config.mc) * config.kc;
	size_t p = static_cast<size_t>(config.mp) * config.kp;
	SimpleCluster::init_array(norm_c, c);
	SimpleCluster::init_array(norm_r1, p);
	SimpleCluster::init_array(norm_r2, p);
	SimpleCluster::init_array(dot_cr1, config.kc * p);
	SimpleCluster::init_array(dot_cr2, config.kc * p);
	SimpleCluster::init_array(dot_rr, p * config.kp);
	SimpleCluster::init_array(diff_qc, c);
	SimpleCluster::init_array(diff_qr1, p);
	SimpleCluster::init_array(diff_qr2, p);

	float * v_tmp1, * v_tmp2, * v_tmp3, * v_tmp4, * v_tmp5;
	float d_tmp, d, d2;
	size_t i, j, i1, j1, k;
	size_t base, base_c, base_p;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	// Calculate all coarse center norms
	v_tmp1 = cq;
	for(i = 0; i < c; i++) {
		d = 0.0;
		for(j = 0; j < bsc; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_c[i] = d;
	}

	// Calculate all product center norms
	v_tmp1 = pq;
	for(i = 0; i < p; i++) {
		d = 0.0;
		for(j = 0; j < bsp; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_r1[i] = d;
	}
	v_tmp1 = rq;
	for(i = 0; i < p; i++) {
		d = 0.0;
		for(j = 0; j < bsp; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_r2[i] = d;
	}

	// Calculate all product center norms
	v_tmp1 = pq;
	v_tmp2 = rq;
	base = 0;
	for(i1 = 0; i1 < config.mp; i1++) {
		for(i = 0; i < config.kp; i++) {
			v_tmp3 = v_tmp2;
			for(j = 0; j < config.kp ; j++) {
				d = 0.0;
				for(k = 0; k < bsp; k++) {
					d_tmp = v_tmp1[k] * *(v_tmp3++);
					d += d_tmp;
				}
				dot_rr[base++] = 2.0f * d;
			}
			v_tmp1 += bsp;
		}
		v_tmp2 = v_tmp3;
	}
	assert(base == p * config.kp);

	// Calculate all dot-products
	base = 0;
	v_tmp1 = cq;
	v_tmp2 = pq;
	v_tmp4 = rq;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < config.kc; j++) {
			v_tmp3 = v_tmp2;
			v_tmp5 = v_tmp4;
			for(i1 = 0; i1 < config.mp/config.mc; i1++) {
				for(j1 = 0; j1 < config.kp; j1++) {
					d = d2 = 0.0;
					for(k = 0; k < bsp; k++) {
						d += v_tmp1[k] * v_tmp3[k];
						d2 += v_tmp1[k] * v_tmp5[k];
					}
					dot_cr1[base] = 2.0f * d;
					dot_cr2[base++] = 2.0f * d2;
					v_tmp3 += bsp;
					v_tmp5 += bsp;
				}
				v_tmp1 += bsp;
			}
		}
		v_tmp2 += config.kp * bsc;
		v_tmp4 += config.kp * bsc;
	}
	assert(base == p * config.kc);
}

/**
 * Precompute all things that are query dependent
 */
inline void PQRQuery::pre_compute2(float * query) {
	size_t i, j, k;
	float *  v_tmp1, * v_tmp2, * v_tmp3;
	float d, d_tmp = 0.0;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	v_tmp1 = query;
	v_tmp2 = cq;
	size_t base = 0;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < config.kc; j++) {
			d = cblas_sdot(bsc,v_tmp1,1,v_tmp2,1);
			diff_qc[base] = norm_c[base++] - d - d;
			v_tmp2 += bsc;
		}
		v_tmp1 += bsc;
	}

	v_tmp1 = query;
	v_tmp2 = pq;
	v_tmp3 = rq;
	base = 0;
	for(i = 0; i < config.mp; i++) {
		for(j = 0; j < config.kp; j++) {
			d = cblas_sdot(bsp,v_tmp1,1,v_tmp2,1);
			diff_qr1[base] = norm_r1[base] - d - d;
			d = cblas_sdot(bsp,v_tmp1,1,v_tmp3,1);
			diff_qr2[base] = norm_r2[base++] - d - d;
			v_tmp2 += bsp;
			v_tmp3 += bsp;
		}
		v_tmp1 += bsp;
	}
}

/**
 * Search method: A demo on single thread mode
 * @param query the query vector
 * @param result the result by identifiers
 * @param R the number of top retrieved results
 * @param verbose to enable verbose mode
 */
inline void PQRQuery::search_ivfadc(float * query,
		float *& v_tmp, float *& dist,
		int *& result, int *& buckets, int *& prebuck,
		int& sum, int R, int w, int T, bool real_dist, bool verbose) {
	if(config.mc != 1) {
		cerr << "This search method is for IVFADC only" << endl;
		return;
	}

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	float * v_tmp2;
	int * i_tmp;
	unsigned char * c_tmp1 = codes;
	unsigned char * c_tmp2 = codes + config.mp;

	// Step 1: assign the query to coarse quantizer
	int i, j, k, l, count = 0, count2,
			base = 0, base1, base_c1, base_c2, c1, c2, bid, sum2 = 0; // 8 * 4 = 32 bytes
	float d_tmp, d_tmp1; // 4 bytes

	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}

	v_tmp1 = diff_qc;
	for(i = 0; i < config.kc; i++) {
		d_tmp = q_sum + *(v_tmp1++);
		pq_insert(v_tmp,buckets,d_tmp,i,sum2,config.kc,verbose);
	}

	if(verbose) {
		cout << "Finished STEP 1" << endl;
		for(i = 0; i < config.kc; i++) {
			cout << buckets[i] << " ";
		}
		cout << endl;
	}


	// Step 2: Local search
	sum2 = config.kc;
	sum = 0;

	for(i = 0; i < w; i++) {
		bid = pq_pop(v_tmp,buckets,sum2,verbose);
		if(bid > 0) {
			l = L[bid] - L[bid-1];
		} else {
			l = L[bid];
		}
		prebuck[count++] = bid;
		sum += l;
		if(sum >= T) break;
	}

	// Allocate the memory to store search results
	// Remember to free them after used
	SimpleCluster::init_array(result,sum << 1);
	SimpleCluster::init_array(dist,sum << 1);
	count = 0;
	int bs = config.dim / config.mp, bsz = bs * config.kp
			,bs1 = config.kp * config.mp / config.mc; // 8 bytes

	for(i = 0; i < w; i++) {
		bid = prebuck[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp1 = codes + static_cast<size_t>(config.mp << 1) *
					static_cast<size_t>(L[bid-1]);
		} else {
			l = L[0];
			i_tmp = pid;
			c_tmp1 = codes;
		}
		c_tmp2 = c_tmp1 + config.mp;

		if(verbose) {
			cout << "Searching in bucket " << bid
					<< " that has " << l << " elements" << endl;
		}

		d_tmp1 = q_sum + diff_qc[bid];
		base1 = bid * bs1;

		// Calculate all l distances
		memcpy(&result[count],i_tmp,l * sizeof(int));
		memcpy(&result[count+sum],i_tmp,l * sizeof(int));
		count2 = count;
		if(!real_dist) {
			for(j = 0; j < l; j++) {
				base = 0;
				d_tmp = d_tmp1;
				for(k = 0; k < config.mp; k++) {
					c2 = *(c_tmp2++);
					base_c1 = base + *(c_tmp1++);
					base_c2 = base + c2;
					d_tmp += (diff_qr1[base_c1] + dot_cr1[base1 + base_c1]);
					d_tmp += (diff_qr2[base_c2] + dot_cr2[base1 + base_c2]);
					d_tmp += dot_rr[base_c1 * config.kp + c2];
					base += config.kp;
				}
				dist[count++] = d_tmp;
				c_tmp1 = c_tmp2;
				c_tmp2 = c_tmp1 + config.mp;
			}
		} else {
			for(j = 0; j < l; j++) {
				dist[count++] = SimpleCluster::distance_l2(
						query,raw_data + i_tmp[j] * config.dim,config.dim);
			}
		}
		memcpy(&dist[count2+sum],&dist[count2],l*sizeof(float));
		if(verbose)
			cout << "Searched all " << l << " elements" << endl;
		if(count >= T) break;
	}

	// Step 3: Extract the top R
	if(sum >= R) {
		nth_element(dist+sum,dist+sum+R-1,dist + (sum << 1));
		count = 0;
		for(i = 0; i < sum; i++) {
			if(dist[i] <= dist[sum+R-1]) {
				result[count++] = result[i+sum];
			}
		}
	}
}
} /* namespace PQLearn */

#endif /* R_QUERY_H_ */
