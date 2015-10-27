/*
 * r_encoder.cpp
 *
 *  Created on: 2014/12/31
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include "r_encoder.h"
#include <cassert>
#include <iostream>

using namespace std;

namespace PQLearn {
void REncoder::load_encoded_data(
		const char * filename,
		bool verbose) {
#ifdef _WIN32
#else
	int fd = open(filename, O_RDONLY);
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file" << endl;
		exit(1);
	}

	struct stat s;
	int status = fstat(fd, &s);
	if(status < 0) {
		if(verbose)
			cerr << "Cannot get statistics of file" << endl;
		exit(1);
	}

	size_t f_size = s.st_size; // The size of file

	unsigned char * mapped, * temp;
	/* Mapping the file */
	mapped = (unsigned char *)mmap(0, f_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if(mapped == MAP_FAILED) {
		if(verbose)
			cerr << "Cannot map the file " << filename << endl;
		exit(1);
	}

	int i, j, k;
	size_t base_pid = 0, base_c = 0, base_code = 0;
	temp = mapped;
	int not_empty, p, p1, l;

	// Read the number of buckets and the size of database
	memcpy(&non_empty_bucket,temp,sizeof(int));
	cout << "The number of non empty buckets: " << non_empty_bucket << endl;
	temp += sizeof(int);
	memcpy(&(config.N),temp,sizeof(int));
	temp += sizeof(int);

	// Memory allocation
	SimpleCluster::init_array(cid,static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mc)); // Only restore the length of non-empty buckets
	SimpleCluster::init_array(codes,static_cast<size_t>(config.N)
			* static_cast<size_t>(config.mp + config2.mp)); // Only restore the length of non-empty buckets
	SimpleCluster::init_array(L,size);
	SimpleCluster::init_array(pid,static_cast<size_t>(config.N) * sizeof(int));
	ushort c[config.mc];

	for(i = 0; i < size; i++) {
		// Read the length of the bucket
		memcpy(&l,temp,sizeof(int));
		L[i] = l;
		temp += sizeof(int);
		if(verbose)
			cout << "This bucket contains " << l << " vector(s)" << endl;

		p1 = i;

		for(j = config.mc - 1; j >= 0; --j) {
			c[j] = p1 % config.kc;
			p1 = (p1 - c[j]) / config.kc;
		}

		memcpy(&pid[base_pid],temp,l * sizeof(int));
		base_pid += l;
		temp += l * sizeof(int);

		for(j = 0; j < l; j++) {
			for(k = 0; k < config.mc; k++) {
				cid[base_c++] = c[k];
			}
		}
		for(j = 0; j < l; j++) {
			memcpy(&codes[base_code],temp,config.mp * sizeof(unsigned char));
			temp += config.mp * sizeof(unsigned char);
			base_code += (config.mp + config2.mp);
		}
	}

	if(munmap(mapped,f_size) != 0) {
		cerr << "Cannot munmap file data" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	cout << "Read " << config.N << " data  from " << filename << endl;
#endif
}

void REncoder::output(
		const char * db_path,
		const char * db_prefix,
		bool verbose) {
	// Now we output the ivf structure to file
	char fname[256];
	sprintf(fname,"%s/%s_ivf.edat_",db_path,db_prefix);
	// Output the encoded data
#ifdef _WIN32
#else
	int fd = open(fname, O_RDWR | O_CREAT | O_TRUNC,(mode_t)0600); // file description
	if(fd < 0) {
		if(verbose)
			cerr << "Cannot open the file " << fname << endl;
		exit(1);
	}

	// The size of codes file
	size_t f_size = static_cast<size_t>(config.N)
											* static_cast<size_t>(config.mp)
											+ static_cast<size_t>(config.N + size + 2) * static_cast<size_t>(sizeof(int));
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
	memcpy(&fd_map[bytes],&non_empty_bucket,sizeof(int));
	bytes += sizeof(int);
	// and the size of the database
	memcpy(&fd_map[bytes],&N,sizeof(int));
	bytes += sizeof(int);

	// Output the buckets data
	int l, tmp, i, j;
	unsigned char c;
	int * _pid = pid;
	unsigned char * _code2 = codes;
	size_t ss = config.mp;
	for(i = 0; i < size; i++) {
		l = L[i];
		// Write out the length of the bucket
		memcpy(&fd_map[bytes],&l,sizeof(int));
		bytes += sizeof(int);

		// Write the pid
		memcpy(&fd_map[bytes],_pid, l * sizeof(int));
		bytes += l * sizeof(int);
		_pid += l;

		// Write the code
		memcpy(&fd_map[bytes],_code2,ss *  l * sizeof(unsigned char));
		bytes += ss *  l * sizeof(unsigned char);
		_code2 += ss * l;
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
