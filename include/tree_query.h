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
 *  hk_query.h
 *
 *  Created on: 2015/04/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_TREE_QUERY_H_
#define INCLUDE_TREE_QUERY_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cstring>
#include <string>
#include <queue>
#include <climits>
#include <cerrno>
#include <cassert>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#endif
// Using OpenBLAS for better performance
#include <cblas.h>
#include "pq_utilities.h"
#include "pq_algorithm.h"
#include "pq_query.h"

using namespace std;
using namespace SimpleCluster;

namespace PQLearn {

class TreeQuery : public PQQuery {
private:
	float * low_centers;
	float * high_centers;
	float * diff_qch;
	float * norm_ch;
	int high, low;
public:
	TreeQuery();
	virtual ~TreeQuery();

	void load_codebooks(const char *, const char *, bool);
	inline void search(
			float *, int, int, int, int,
			float *&, int *&, float *&, int *&,
			float *&, int *&, float *&, int *&,
			int &, double &, double &, double &,
			bool,
			bool);
	inline void search_mem(
			float *, int, int, int, int,
			float *&, int *&, float *&, int *&,
			float *&, int *&, float *&, int *&,
			int &, double &, double &, double &,
			bool,
			bool);
	inline void pre_compute1();
	inline void pre_compute2(float *);
	inline void pre_compute3(float *, int);
	inline void pre_compute4(float *);
	inline void pre_compute5();
	//void search_multi(float *, bool);
};

inline void TreeQuery::pre_compute1() {
	size_t i, j, k;
	float * tmp1 = low_centers;
	float * tmp2 = high_centers;
	for(i = 0; i < high; i++) {
		for(j = 0; j < low; j++) {
			for(k = 0; k < this->config.dim; k++) {
				tmp1[k] += tmp2[k];
			}
			tmp1 +=  this->config.dim;
		}
		tmp2 += this->config.dim;
	}

	init_array(norm_c, config.mc * high * low);
	init_array(norm_ch, config.mc * high);
	init_array(norm_r, config.mp * config.kp);
	init_array(diff_qc, config.mc * high * low);
	init_array(diff_qch, config.mc * high);
	init_array(diff_qr, config.mp * config.kp);

	float * v_tmp1, * v_tmp2, * v_tmp3;
	float d_tmp, d;
	size_t i1, j1;
	size_t base, base_c, base_p;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	// Calculate all coarse center norms
	v_tmp1 = low_centers;
	for(i = 0; i < config.mc * high * low; i++) {
		d = 0.0;
		for(j = 0; j < bsc; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_c[i] = d;
	}

	v_tmp1 = high_centers;
	for(i = 0; i < config.mc * high; i++) {
		d = 0.0;
		for(j = 0; j < bsc; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_ch[i] = d;
	}

	// Calculate all product center norms
	v_tmp1 = pq;
	for(i = 0; i < config.mp * config.kp; i++) {
		d = 0.0;
		for(j = 0; j < bsp; j++) {
			d_tmp = *(v_tmp1++);
			d += d_tmp * d_tmp;
		}
		norm_r[i] = d;
	}
}

inline void TreeQuery::pre_compute2(float * q) {
	size_t i, j, k;
	float *  v_tmp1, * v_tmp2;
	float d, d_tmp = 0.0;
	size_t base;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	//	v_tmp1 = q;
	//	v_tmp2 = low_centers;
	//	base = 0;
	//	for(i = 0; i < config.mc; i++) {
	//		for(j = 0; j < high * low; j++) {
	//			d = cblas_sdot(bsc,v_tmp1,1,v_tmp2,1);
	//			diff_qc[base] = norm_c[base++] - d - d;
	//			v_tmp2 += bsc;
	//		}
	//		v_tmp1 += bsc;
	//	}

	v_tmp1 = q;
	v_tmp2 = high_centers;
	base = 0;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < high; j++) {
			d = cblas_sdot(bsc,v_tmp1,1,v_tmp2,1);
			diff_qch[base] = norm_ch[base++] - d - d;
			v_tmp2 += bsc;
		}
		v_tmp1 += bsc;
	}

	v_tmp1 = q;
	v_tmp2 = pq;
	base = 0;
	for(i = 0; i < config.mp; i++) {
		for(j = 0; j < config.kp; j++) {
			d = cblas_sdot(bsp,v_tmp1,1,v_tmp2,1);
			diff_qr[base] = norm_r[base++] - d - d;
			v_tmp2 += bsp;
		}
		v_tmp1 += bsp;
	}
}

inline void TreeQuery::pre_compute3(float * q, int h) {
	size_t i, j, k;
	float *  v_tmp1, * v_tmp2;
	float d, d_tmp = 0.0;
	size_t base;
	int bsc = config.dim / config.mc;

	v_tmp1 = q;
	v_tmp2 = low_centers + h * low * config.dim;
	base = h * low;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < low; j++) {
			d = cblas_sdot(bsc,v_tmp1,1,v_tmp2,1);
			diff_qc[base] = norm_c[base++] - d - d;
			v_tmp2 += bsc;
		}
		v_tmp1 += bsc;
	}
}

inline void TreeQuery::pre_compute4(float * q) {
	size_t i, j, k;
	float *  v_tmp1, * v_tmp2;
	float d, d_tmp = 0.0;
	size_t base;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	v_tmp1 = q;
	v_tmp2 = pq;
	base = 0;
	for(i = 0; i < config.mp; i++) {
		for(j = 0; j < config.kp; j++) {
			d = cblas_sdot(bsp,v_tmp1,1,v_tmp2,1);
			diff_qr[base] = norm_r[base++] - d - d;
			v_tmp2 += bsp;
		}
		v_tmp1 += bsp;
	}
}

/**
 * Memory consumed pre-computation
 */
inline void TreeQuery::pre_compute5() {
	size_t i, j, k;
	init_array(dot_cr, static_cast<size_t>(high * low) * config.kp * config.mp);

	float * v_tmp1, * v_tmp2, * v_tmp3;
	float d_tmp, d;
	size_t i1, j1;
	size_t base, base_c, base_p;
	int bsc = config.dim / config.mc;
	int bsp = config.dim / config.mp;

	// Calculate all dot-products
	base = 0;
	v_tmp1 = low_centers;
	v_tmp2 = pq;
	for(i = 0; i < config.mc; i++) {
		for(j = 0; j < high * low; j++) {
			v_tmp3 = v_tmp2;
			for(i1 = 0; i1 < config.mp/config.mc; i1++) {
				for(j1 = 0; j1 < config.kp; j1++) {
					d = 0.0;
					for(k = 0; k < bsp; k++) {
						d += v_tmp1[k] * v_tmp3[k];
					}
					dot_cr[base++] = 2.0f * d;
					v_tmp3 += bsp;
				}
				v_tmp1 += bsp;
			}
		}
		v_tmp2 += config.kp * bsc;
	}
	assert(base == static_cast<size_t>(high) * low * config.kp * config.mp);
}

inline void TreeQuery::search(
		float * query,
		int T,
		int m, // extract top m high-level cells
		int n, // extract top n low-level cells in each high-level cells
		int R,
		float *& dist_h, int *& id_h,
		float *& dist_l, int *& id_l,
		float *& dists, int *& ids,
		float *& distf, int *& idf,
		int & sum, double& t0, double& t1, double& t2,
		bool real_dist,
		bool verbose) {
	if(config.mc != 1) {
		cerr << "This search method is for IVFADC only" << endl;
		return;
	}

	if(n >= low) n = low;
	if(m >= high) m = high;

	if(T <= 0) T = config.N;

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	float * v_tmp2;
	int * i_tmp;
	unsigned char * c_tmp = codes;
	float * q = (float *)::operator new(config.dim * sizeof(float));

	// First, extract top m high-level cells
	size_t i, j, k, l, count = 0, count2,
			base = 0, base1, base2, base_c, c, bid; // 8 * 4 = 32 bytes
	float d_tmp, d_tmp1; // 4 bytes
	clock_t st, ed;

	st = clock();
	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}
	ed = clock();
	t0 += (ed - st);

	st = clock();
	v_tmp1 = diff_qch;
	for(i = 0; i < high; i++) {
		dist_h[i] = q_sum + v_tmp1[i];
		id_h[i] = i;
	}

	nth_element_id(dist_h,dist_h + high,id_h,m-1);
	sort_id(dist_h,dist_h + m, id_h);
	ed = clock();

	t1 += (ed - st);

	if(verbose) {
		for(i = 0; i < m; i++) {
			cout << "pos " << i << ":(" << id_h[i] << "," << dist_h[i] << ")" << endl;
		}
	}

	// Second, extract top n low-level cells in each high-level cells in top m
	base = 0;
	int base_h = 0;
	sum = 0;
	for(i = 0 ; i < m; i++) {
		st = clock();
		pre_compute3(query,id_h[i]);
		ed = clock();
		t0 += (ed -st);
		st = clock();
		base_h = id_h[i] * low;
		v_tmp1 = diff_qc + base_h;
		for(j = 0; j < low; j++) {
			dist_l[j] = q_sum + v_tmp1[j];
			id_l[j] = j;
		}

		nth_element_id(dist_l,dist_l + low,id_l,n-1);
		sort_id(dist_l,dist_l + n, id_l);
		for(j = 0; j < n; j++) {
			count = ids[base] = base_h + id_l[j];
			dists[base] = dist_l[j];
			if(count == 0) sum += L[0];
			else sum += (L[count] - L[count-1]);
			base++;
			if(sum >= T) break;
		}
		ed = clock();
		t1 += (ed -st);
		if(sum >= T) break;
	}

	if(verbose) {
		for(i = 0; i < base; i++) {
			cout << "Bucket " << i << ":(" << ids[i] << "," << dists[i] << ")" << endl;
		}
		cout << "We have " << sum << " candidate(s) in " << base << " buckets" << endl;
	}

	// Perform the search in specified low-level cells
	// Allocate the memory to store search results
	// Remember to free them after used
	init_array(distf,sum);
	init_array(idf,sum);
	count = 0;
	int bs = config.dim / config.mp, bsz = bs * config.kp
			,bs1 = config.kp * config.mp / config.mc; // 8 bytes
	int limit = config.kp * config.dim / (config.dim - config.mp);

	st = clock();
	for(i = 0; i < base; i++) {
		bid = ids[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp = codes + static_cast<size_t>(config.mp) *
					static_cast<size_t>(L[bid-1]);
		} else {
			l = L[0];
			i_tmp = pid;
			c_tmp = codes;
		}

		if(verbose) {
			cout << i << ":Searching in bucket " << bid
					<< " that has " << l << " elements" << endl;
		}

		// Subtracting q = query - low_centers[bid]
		v_tmp1 = low_centers + bid * config.dim;
		d_tmp1 = 0.0;
		for(j = 0; j < config.dim; j++) {
			q[j] = query[j] - *(v_tmp1++);
			d_tmp1 += q[j] * q[j];
		}

		// Compute the cache values
		if(l >= limit) {
			pre_compute4(q);
		}

		// Calculate all l distances
		memcpy(&idf[count],i_tmp,l * sizeof(int));

		if(!real_dist) {
			if(l >= limit) {
				for(j = 0; j < l; j++) {
					base2 = 0;
					d_tmp = d_tmp1;
					for(k = 0; k < config.mp; k++) {
						c = static_cast<int>(*(c_tmp++));
						base_c = base2 + c;
						d_tmp += diff_qr[base_c];
						base2 += config.kp;
					}
					distf[count] = d_tmp;
					count++;
				}
			} else {
				for(j = 0; j < l; j++) {
					base2 = 0;
					v_tmp1 = q;
					d_tmp = 0.0;
					for(k = 0; k < config.mp; k++) {
						c = static_cast<int>(*(c_tmp++));
						base_c = base2 + c;
						d_tmp += distance_l2(v_tmp1, pq + base_c * bs, bs);
						base2 += config.kp;
						v_tmp1 += bs;
					}
					distf[count] = d_tmp;
					count++;
				}
			}
		} else {
			for(j = 0; j < l; j++) {
				distf[count++] = distance_l2(
						query,raw_data + i_tmp[j] * config.dim,config.dim);
			}
		}
		if(verbose)
			cout << "Searched all " << count << " elements" << endl;
		if(count >= sum) break;
	}
	ed = clock();
	t2 += (ed -st);

	if(verbose) {
		for(i = 0; i < count; i++) {
			cout << "Candidate " << i << " : (" << idf[i] << "," << distf[i] << ")" << endl;
		}
	}

	// Step 3: Extract the top R
	st = clock();
	if(sum >= R) {
		nth_element_id(distf,distf + sum,idf,R-1);
		sort_id(distf,distf + R, idf);
	} else sort_id(distf,distf+sum,idf);
	ed = clock();
	t1 += (ed - st);
}

inline void TreeQuery::search_mem(
		float * query,
		int T,
		int m, // extract top m high-level cells
		int n, // extract top n low-level cells in each high-level cells
		int R,
		float *& dist_h, int *& id_h,
		float *& dist_l, int *& id_l,
		float *& dists, int *& ids,
		float *& distf, int *& idf,
		int & sum, double& t0, double& t1, double& t2,
		bool real_dist,
		bool verbose) {
	if(config.mc != 1) {
		cerr << "This search method is for IVFADC only" << endl;
		return;
	}

	if(n >= low) n = low;
	if(m >= high) m = high;

	if(T <= 0) T = config.N;

	// Temporary pointers: 8 * 8 = 64 bytes
	float * v_tmp1;
	float * v_tmp2;
	int * i_tmp;
	unsigned char * c_tmp = codes;
	float * q = (float *)::operator new(config.dim * sizeof(float));

	// First, extract top m high-level cells
	size_t i, j, k, l, count = 0, count2,
			base = 0, base1, base2, base_c, c, bid; // 8 * 4 = 32 bytes
	float d_tmp, d_tmp1; // 4 bytes
	clock_t st, ed;

	st = clock();
	pre_compute2(query);
	float q_sum = 0.0;
	v_tmp1 = query;
	for(i = 0; i < config.dim; i++) {
		d_tmp = *(v_tmp1++);
		q_sum += d_tmp * d_tmp;
	}
	ed = clock();
	t0 += (ed - st);

	st = clock();
	v_tmp1 = diff_qch;
	for(i = 0; i < high; i++) {
		dist_h[i] = q_sum + v_tmp1[i];
		id_h[i] = i;
	}

	nth_element_id(dist_h,dist_h + high,id_h,m-1);
	sort_id(dist_h,dist_h + m, id_h);
	ed = clock();

	t1 += (ed - st);

	if(verbose) {
		for(i = 0; i < m; i++) {
			cout << "pos " << i << ":(" << id_h[i] << "," << dist_h[i] << ")" << endl;
		}
	}

	// Second, extract top n low-level cells in each high-level cells in top m
	base = 0;
	int base_h = 0;
	sum = 0;
	for(i = 0 ; i < m; i++) {
		st = clock();
		pre_compute3(query,id_h[i]);
		ed = clock();
		t0 += (ed -st);
		st = clock();
		base_h = id_h[i] * low;
		v_tmp1 = diff_qc + base_h;
		for(j = 0; j < low; j++) {
			dist_l[j] = q_sum + v_tmp1[j];
			id_l[j] = j;
		}

		nth_element_id(dist_l,dist_l + low,id_l,n-1);
		sort_id(dist_l,dist_l + n, id_l);
		for(j = 0; j < n; j++) {
			count = ids[base] = base_h + id_l[j];
			dists[base] = dist_l[j];
			if(count == 0) sum += L[0];
			else sum += (L[count] - L[count-1]);
			base++;
			if(sum >= T) break;
		}
		ed = clock();
		t1 += (ed -st);
		if(sum >= T) break;
	}

	if(verbose) {
		for(i = 0; i < base; i++) {
			cout << "Bucket " << i << ":(" << ids[i] << "," << dists[i] << ")" << endl;
		}
		cout << "We have " << sum << " candidate(s) in " << base << " buckets" << endl;
	}

	// Perform the search in specified low-level cells
	// Allocate the memory to store search results
	// Remember to free them after used
	init_array(distf,sum);
	init_array(idf,sum);
	count = 0;
	int bs = config.dim / config.mp, bsz = bs * config.kp
			,bs1 = config.kp * config.mp / config.mc; // 8 bytes
	int limit = config.kp * config.dim / (config.dim - config.mp);

	st = clock();
	for(i = 0; i < base; i++) {
		bid = ids[i];
		if(bid > 0) {
			l = L[bid] - L[bid-1];
			i_tmp = pid + L[bid-1];
			c_tmp = codes + static_cast<size_t>(config.mp) *
					static_cast<size_t>(L[bid-1]);
		} else {
			l = L[0];
			i_tmp = pid;
			c_tmp = codes;
		}

		if(verbose) {
			cout << i << ":Searching in bucket " << bid
					<< " that has " << l << " elements" << endl;
		}

		d_tmp1 = dists[i];
		base1 = bid * bs1;

		// Calculate all l distances
		memcpy(&idf[count],i_tmp,l * sizeof(int));

		for(j = 0; j < l; j++) {
			base2 = 0;
			d_tmp = d_tmp1;
			for(k = 0; k < config.mp; k++) {
				c = static_cast<int>(*(c_tmp++));
				base_c = base2 + c;
				d_tmp += (diff_qr[base_c] + dot_cr[base1 + base_c]);
				base2 += config.kp;
			}
			distf[count] = d_tmp;
			count++;
		}

		if(verbose)
			cout << "Searched all " << count << " elements" << endl;
		if(count >= T) break;
	}
	ed = clock();
	t2 += (ed -st);

	if(verbose) {
		for(i = 0; i < count; i++) {
			cout << "Candidate " << i << " : (" << idf[i] << "," << distf[i] << ")" << endl;
		}
	}

	// Step 3: Extract the top R
	st = clock();
	if(count >= R) {
		nth_element_id(distf,distf + count,idf,R-1);
		sort_id(distf,distf + R, idf);
	} else sort_id(distf,distf+sum,idf);
	ed = clock();
	t1 += (ed - st);
}

} /* namespace PQLearn */

#endif /* INCLUDE_TREE_QUERY_H_ */
