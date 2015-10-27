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
 *  hk_encoder.cpp
 *
 *  Created on: 2015/04/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cmath>

#include "tree_encoder.h"
#include "encoder.h"

using namespace std;

namespace PQLearn {

TreeEncoder::TreeEncoder() : Encoder::Encoder() {
	tree_cb = nullptr;
	tree_codes = nullptr;
	hlb = nullptr;
	llb = nullptr;
}

TreeEncoder::~TreeEncoder() {
	::delete tree_cb;
	tree_cb = nullptr;
	::delete tree_codes;
	tree_codes = nullptr;
	::delete hlb;
	hlb = nullptr;
	::delete llb;
	llb = nullptr;
}

void TreeEncoder::load_codebooks(
		const char * cb_path,
		const char * pq_path,
		bool verbose) {
	load_codebook<float>(cb_path,config,tree_cb,3,verbose);
	load_codebook<float>(pq_path,config,pq,1,verbose);
	size = static_cast<int>(config.kc * config.L);
	if(verbose)
		cout << "The size of ivf:" << size << endl;
	cout << "--> Settings: (kc,mc,L,kp,mp)=" << config.kc << " "
			<< config.mc << " " << config.L << " " << config.kp << " "
			<< config.mp << endl;
}

void TreeEncoder::distribution(bool verbose) {
	size_t i, j;
	size_t base_u = 0, base_c = 0, base_p = 0;
	non_empty_bucket = 0;
	int id;
	TreeBucket bk_tmp;
	tree_ivf.clear();
	for(i = 0; i < size; i++) {
		tree_ivf.push_back(bk_tmp);
	}

	for(i = 0; i < config.N; i++) {
		// Check whether if the bucket was created or not?
		id = hlb[i] * config.L + llb[i];
		// The bucket
		if(id < size) {
			// Set the bucket id
			tree_ivf[id].pid.insert(tree_ivf[id].pid.end(),base_p++);
			tree_ivf[id].codes.insert(tree_ivf[id].codes.end(),
					&tree_codes[base_c], &tree_codes[base_c + config.mp]);
			tree_ivf[id].L++;
			if(tree_ivf[id].L == 1) non_empty_bucket++;
		}
		base_c += config.mp;
	}
}

void TreeEncoder::output(
		const char * filename,
		bool verbose) {
#ifdef _WIN32
#else
	int fd = open(filename, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file " << filename << endl;
		exit(1);
	}

	// The size of codes file
	size_t f_size = static_cast<size_t>(config.N)
					* static_cast<size_t>(config.mp)
					+ static_cast<size_t>(config.N + size + 2)
					* static_cast<size_t>(sizeof(int));
	int result = lseek(fd, f_size, SEEK_SET);
	if (result == -1) {
		close(fd);
		if(verbose)
			cerr << "Error calling lseek() to 'stretch' the file" <<  endl;
		exit(1);
	}

	int status = write(fd, "", 1);
	if(status != 1) {
		if(verbose)
			cerr << "Cannot write to file" << endl;
		exit(1);
	}

	unsigned char * fd_map = (unsigned char *)mmap(0, f_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (fd_map == MAP_FAILED) {
		close(fd);
		if(verbose)
			cerr << "Error mmapping the file" << endl;
		exit(1);
	}
	// Output the encoded data
	size_t bytes = 0;
	// The size of codes
	int N = config.N;
	// The first 8 bytes will be
	// the number of buckets
	memcpy(&fd_map[bytes],&(non_empty_bucket),sizeof(int));
	bytes += sizeof(int);
	// and the size of the database
	memcpy(&fd_map[bytes],&N,sizeof(int));
	bytes += sizeof(int);

	int L, tmp, i, j;
	unsigned char c;
	vector<int> pid;
	vector<unsigned char> code2;
	for(i = 0; i < size; i++) {
		L = tree_ivf[i].L;
		if(verbose) {
			cout << "This bucket has "<< L << " element(s)" << endl;
		}
		// Write out the length of the bucket
		memcpy(&fd_map[bytes],&L,sizeof(int));
		bytes += sizeof(int);

		// Write the pid
		pid = tree_ivf[i].pid;
		for(j = 0; j < L; j++) {
			tmp = pid[j];
			memcpy(&fd_map[bytes],&tmp,sizeof(int));
			bytes += sizeof(int);
		}

		// Write the code
		code2 = tree_ivf[i].codes;
		if(code2.size() != L * config.mp) {
			cerr << "Wrong data" << endl;
			exit(EXIT_FAILURE);
		}
		for(j = 0; j < code2.size(); j++) {
			c = code2[j];
			memcpy(&fd_map[bytes],&c,sizeof(unsigned char));
			bytes += sizeof(unsigned char);
		}
		pid.clear();
		code2.clear();
	}

	if (munmap(fd_map, f_size) == -1) {
		if(verbose)
			cerr << "Error un-mmapping the file" << endl;
	}
	close(fd);

	cout << "Read " << bytes << " byte(s) and wrote out " << f_size << " byte(s)" << endl;
#endif
}

} /* namespace PQLearn */
