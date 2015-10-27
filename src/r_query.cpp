/*
 * r_query.cpp
 *
 *  Created on: 2015/01/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include "r_query.h"
#include <iostream>
#include <cassert>

using namespace std;

namespace PQLearn
{

PQRQuery::PQRQuery() : PQQuery::PQQuery() {
	rq = nullptr;
	diff_qr1 = nullptr;
	diff_qr2 = nullptr;
	dot_cr1 = nullptr;
	dot_cr2 = nullptr;
	dot_rr = nullptr;
	norm_r1 = nullptr;
	norm_r2 = nullptr;
}

PQRQuery::~PQRQuery(){
	::delete rq;
	::delete diff_qr1;
	::delete diff_qr2;
	::delete dot_cr1;
	::delete dot_cr2;
	::delete dot_rr;
	::delete norm_r1;
	::delete norm_r2;
	rq = nullptr;
	diff_qr1 = nullptr;
	diff_qr2 = nullptr;
	dot_cr1 = nullptr;
	dot_cr2 = nullptr;
	dot_rr = nullptr;
	norm_r1 = nullptr;
	norm_r2 = nullptr;
}

void PQRQuery::load_codebooks(
		const char * cq_path,
		const char * pq_path,
		const char * rq_path,
		bool verbose) {
	load_codebook<float>(cq_path,config,cq,0,verbose);
	load_codebook<float>(pq_path,config,pq,1,verbose);
	load_codebook<float>(rq_path,config2,rq,1,verbose);
	size = static_cast<int>(pow(config.kc,config.mc));
	cout << "--> Settings: (kc,mc,kp,mp)=" << config.kc << " "
			<< config.mc << " " << config.kp << " " << config.mp << endl;
}

/**
 * Load the encoded data from binary file
 * @param filename path to the encoded data file
 * @param verbose enable verbose mode
 */
int PQRQuery::load_encoded_data(const char * filename, bool verbose) {
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

	size_t l, i, j, k, count = 0, tmp;
	size_t base_pid = 0, base_code = 0;
	temp = mapped;

	// Read the number of buckets and the size of database
	memcpy(&not_empty,temp,sizeof(int));
	cout << "The number of non empty buckets: " << not_empty << endl;
	temp += sizeof(int);
	memcpy(&(config.N),temp,sizeof(int));
	temp += sizeof(int);

	size_t p = config.mp << 1;

	// Memory allocation
	// SimpleCluster::init_array<int>(bits,size >> 5 + 1); // The list of non-empty buckets as a bit set
	SimpleCluster::init_array(L,size); // Only restore the length of non-empty buckets
	SimpleCluster::init_array(pid,config.N);
	codes = (unsigned char *)::operator new(static_cast<size_t>(config.N) *
			p * sizeof(unsigned char));

	for(i = 0; i < size; i++) {
		// Read the length of the bucket
		memcpy(&l,temp,sizeof(int));
		temp += sizeof(int);
		//		if(l > 0) {
		//			tmp = i & 32;
		//			bits[i>>5] |= (1 << (31 - tmp));
		//			L[count++] = l;
		//		}
		if(i > 0) L[i] = l + L[i-1];
		else L[i] = l;
		if(verbose)
			cout << "This bucket contains " << l << " vector(s)" << endl;

		if(l > 0) {
			// Read the pid
			memcpy(pid + base_pid,temp,l * sizeof(int));
			temp += l * sizeof(int);
			base_pid += l;

			// Read the codes
			memcpy(codes + base_code,temp,l * p * sizeof(unsigned char));
			temp += l * p * sizeof(unsigned char);
			base_code += l * p;
		}
	}
	assert(base_code == p * config.N);


	if(munmap(mapped,f_size) != 0) {
		cerr << "Cannot munmap file data" << endl;
		exit(EXIT_FAILURE);
	}
	close(fd);

	cout << "Read " << config.N << " data  from " << filename << endl;

	return config.N;
#endif
}

} /* namespace PQLearn */
