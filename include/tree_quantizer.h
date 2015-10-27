/*
 *  SIMPLE CLUSTERS: A simple library for clustering works.
 *  Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
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
 *  hk_quantizer.h
 *
 *  Created on: 2015/04/01
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_TREE_QUANTIZER_H_
#define INCLUDE_TREE_QUANTIZER_H_

#include <iostream>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "quantizer.h"

using namespace std;

namespace PQLearn {

template<typename DataType>
class TreeQuantizer : public PQQuantizer<DataType> {
private:
	int H, L;
	float * low_centers;
	int * low_labels;
	unsigned int * low_ulabels;
	float * low_data;
public:
	TreeQuantizer(int,int,int,int,int,bool);
	virtual ~TreeQuantizer();

	inline void create_low_level_cells(bool);
	inline void reload_high_level_cells(const char *,bool);
	inline void tree_output(const char *, bool);
};

template<typename DataType>
TreeQuantizer<DataType>::TreeQuantizer (
		int _dim,
		int _m,
		int _h,
		int _l,
		int _os,
		bool verbose) : PQQuantizer<DataType>::PQQuantizer(_dim,_m,_h,_os,verbose) {
	H = _h;
	L = _l;
	if(H <= 0 || L <= 0) return;
	low_labels = nullptr;
	low_ulabels = nullptr;
	low_data = nullptr;

	if(!SimpleCluster::init_array<float>(low_centers,static_cast<size_t>(H) * L * this->dim)) {
		if(verbose)
			cerr << "There are some errors occurred while initializing low-level centers" << endl;
		exit(EXIT_FAILURE);
	}
}

template<typename DataType>
TreeQuantizer<DataType>::~TreeQuantizer () {
	::delete low_centers;
	low_centers = nullptr;
	::delete low_labels;
	low_labels = nullptr;
	::delete low_ulabels;
	low_ulabels = nullptr;
	::delete low_data;
	low_data = nullptr;
}

template<typename DataType>
inline void TreeQuantizer<DataType>::create_low_level_cells(bool verbose) {
	if(verbose) {
		cout << "Create low level quantizers for:"<< endl;
		cout << "H=" << this->nsc << endl;
		cout << "L=" << L << endl;
		cout << "N=" << this->N << endl;
	}

	if(this->type == 1) {
		if(!SimpleCluster::init_array<unsigned int>(low_ulabels,this->N)) {
			if(verbose)
				cerr << "There are some errors occurred while initializing low-level labels" << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		if(!SimpleCluster::init_array<int>(low_labels,this->N)) {
			if(verbose)
				cerr << "There are some errors occurred while initializing low-level labels" << endl;
			exit(EXIT_FAILURE);
		}
	}
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif

	int i, j, i0;
	size_t p = static_cast<size_t>(this->N) / static_cast<size_t>(max_threads);

	low_data = (float *)::operator new(static_cast<size_t>(this->N) * this->dim * sizeof(float));
	size_t * tmp_count = (size_t *)::operator new(this->nsc * sizeof(size_t));
	int * tmp_cur = (int *)::operator new(this->nsc * sizeof(int));

	// Count the number of elements in each clusters
	for(i = 0; i < this->nsc; i++) {
		tmp_count[i] = 0;
		tmp_cur[i] = 0;
	}
	if(this->type == 1) {
		for(i = 0; i < this->N; i++) {
			tmp_count[this->ulabels[i]]++;
		}
	} else {
		for(i = 0; i < this->N; i++) {
			tmp_count[this->labels[i]]++;
		}
	}
	for(i = 1; i < this->nsc; i++) {
		tmp_count[i] += tmp_count[i-1];
	}
	for(i = this->nsc - 1; i > 0; i--) {
		tmp_count[i] = tmp_count[i-1];
	}
	tmp_count[0] = 0;

	if(verbose) {
		cout << "Counted all elements" << endl;
		for(i = 0; i < this->nsc; i++) {
			cout << "Group " << i << " starts from " << tmp_count[i] << endl;
		}
	}

	// Calculate the residual vectors
#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(i, j, i0)
#endif
		for(i0 = 0; i0 < max_threads; i0++) {
			size_t start = p * static_cast<size_t>(i0);
			size_t end = start + p;
			if(end > this->N || i0 == max_threads - 1)
				end = this->N;
			float * tmp = this->data + start * this->dim;
			float * tmp1;
			int id;
			for(i = start; i < end; i++) {
				if(this->type == 1)
					id = this->ulabels[i];
				else id = this->labels[i];
				tmp1 = this->centers + id * this->dim;
				for(j = 0; j < this->dim; j++) {
					// Calculate the residual vector
					*tmp -= *tmp1;
					tmp++;
					tmp1++;
				}
			}
		}
#ifdef _OPENMP
	}
#endif

	if(verbose) {
		cout << "Computed the residual vectors" << endl;
	}

	// Regrouping
	float * tmp = this->data;
	float * tmp1;
	int id;
	size_t pos;
	for(i = 0; i < this->N; i++) {
		if(this->type == 1)
			id = this->ulabels[i];
		else id = this->labels[i];
		pos = tmp_count[id] + tmp_cur[id];
		tmp1 = low_data + pos * this->dim;
		memcpy(tmp1,tmp,this->dim * sizeof(float));
		tmp_cur[id]++;
		tmp += this->dim;
	}

	if(verbose) {
		cout << "Regrouped" << endl;
		for(i = 0; i < this->nsc; i++) {
			cout << "Group " << i << " has " << tmp_cur[i] << " elements" << endl;
		}
	}

	// Do the k-means in each groups
	float * _centers, * _seeds = nullptr, * _data;
	int * _labels;
	unsigned int * _ulabels;
	KmeansCriteria criteria = {2.0,1.0,1000};
	_centers = low_centers;
	if(this->type == 1)
		_ulabels = low_ulabels;
	else
		_labels = low_labels;
	_data = low_data;
	float * _distances;
	double energy;
	size_t cl = L * this->dim;
	for(i = 0; i < this->nsc; i++) {
		if(verbose) {
			cout << "Creating codebook " << i  << "/" << this->nsc << endl;
			cout << "Group " << i << " has " << tmp_cur[i] << " elements" << endl;
		}
		if(this->type == 1) {
			vl_kmeans_exec(
					_data, _centers, _ulabels, _distances,
					tmp_cur[i], L, this->dim, 1, VlDistanceL2, energy, verbose);
		} else {
			greg_kmeans<float>(
					_data,_centers,_labels,_seeds,
					KmeansType::KMEANS_PLUS_SEEDS,
					criteria,
					DistanceType::NORM_L2,
					EmptyActs::SINGLETON,
					tmp_cur[i],L,this->dim,max_threads,
					verbose);
		}
		if(verbose) {
			cout << "Finished subcodebook " << i << endl;
		}


		// Calculate the residuals
		if(this->type == 1) {
			for(j = 0; j < tmp_cur[i]; j++) {
				_data[j] -= _centers[_ulabels[j]];
			}
		} else {
			for(j = 0; j < tmp_cur[i]; j++) {
				_data[j] -= _centers[_labels[j]];
			}
		}
		_centers += cl;
		if(this->type == 1)
			_ulabels += tmp_cur[i];
		else _labels += tmp_cur[i];
		_data += static_cast<size_t>(tmp_cur[i]) * this->dim;
	}

	// Remember to update the residual data
	memcpy(this->data,low_data,static_cast<size_t>(this->N) * this->dim * sizeof(float));
	::delete low_data;
	low_data = nullptr;
}

template<typename DataType>
inline void TreeQuantizer<DataType>::reload_high_level_cells(
		const char * filename,
		bool verbose) {
	PQConfig config;
	load_codebook<float>(filename,config,this->centers,0,verbose);

	this->nsc = config.kc;
	if(verbose) {
		cout << "nsc=" << this->nsc << endl;
		cout << "dim=" << this->dim << endl;
		float * tmp = this->centers;
		for(int i = 0; i < this->dim; i++) {
			cout << *(tmp++) << " ";
		}
		cout << endl;
	}

	// Recalculate the labels
	int i, j, i0;
	int max_threads = 1;
#ifdef _OPENMP
	max_threads = omp_get_max_threads();
#endif
	if(max_threads <= 0) max_threads = 1;
	size_t p = static_cast<size_t>(this->N) / max_threads;

#ifdef _OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel
	{
#pragma omp for private(i, j, i0)
#endif
		for(i0 = 0; i0 < max_threads; i0++) {
			size_t start = p * static_cast<size_t>(i0);
			size_t end = start + p;
			if(end > this->N || i0 == max_threads - 1)
				end = this->N;
			float * tmp = this->data + start * this->dim;

			float * tmp1;
			int id = -1;
			float d_tmp, d;
			for(i = start; i < end; i++) {
				d = FLT_MAX;
				tmp1 = this->centers;
				id = -1;
				for(j = 0; j < this->nsc; j++) {
					d_tmp = SimpleCluster::distance_l2(tmp,tmp1,this->dim);
					if(d >= d_tmp) {
						d = d_tmp;
						id = j;
					}
					tmp1 += this->dim;
				}
				if(this->type == 1)
					this->ulabels[i] = id;
				else this->labels[i] = id;
				tmp += this->dim;
			}
		}
#ifdef _OPENMP
	}
#endif
}

template<typename DataType>
inline void TreeQuantizer<DataType>::tree_output(
		const char * filename,
		bool verbose) {
	char fn[256];
	// First, write the centers to file named filename.ctr_
	sprintf(fn,"%s.ctr_",filename);
#ifdef _WIN32
	ofstream output;
	output.open(fn, ios::out | ios::binary);
	output << this->nsc << this->part << this->dim;
	int i;
	int base = 3;
	float * tmp = this->centers;
	for(int i = 0; i < this->dim * this->nsc; i++) {
		output << *tmp;
		tmp++;
	}
	output.close();
#else
	int fd = open(fn, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		cerr << "Cannot open the file " << fn << endl;
		exit(EXIT_FAILURE);
	}

	size_t size = (static_cast<size_t>(this->nsc) *
			(this->dim + L * this->dim) + 4) * sizeof(float);
	int result = lseek(fd, size, SEEK_SET);
	if (result == -1) {
		close(fd);
		cerr << "Error calling lseek() to 'stretch' the file" <<  endl;
		exit(EXIT_FAILURE);
	}

	int status = write(fd, "", 1);
	if(status != 1) {
		cerr << "Cannot write to file" << endl;
		exit(EXIT_FAILURE);
	}

	float * fd_map = (float *)mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (fd_map == MAP_FAILED) {
		close(fd);
		cerr << "Error mmapping the file" << endl;
		exit(EXIT_FAILURE);
	}
	fd_map[0] = this->nsc;
	fd_map[1] = this->part;
	fd_map[2] = this->dim;
	fd_map[3] = L;
	size_t base = 4;
	float * tmp = this->centers;

	for(int i = 0; i < this->dim * this->nsc; i++) {
		fd_map[base++] = *tmp;
		tmp++;
	}

	tmp = low_centers;
	for(int i = 0; i < H * L * this->dim; i++) {
		fd_map[base++] = *tmp;
		tmp++;
	}

	if (munmap(fd_map, size) == -1) {
		cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);
#endif
}

} /* namespace PQLearn */

#endif /* INCLUDE_TREE_QUANTIZER_H_ */
