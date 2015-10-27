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
 *  test_hk_query.cpp
 *
 *  Created on: 2015/04/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstring>
#include <cstdio>
#include <climits>
#include <cerrno>
#include <gtest/gtest.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

#include "tree_query.h"

using namespace std;
using namespace PQLearn;

char setting[256];

inline void create_log_dir(const char * log_path, time_t timer) {
	struct stat sb;
	char dir[256];
	if(stat(log_path, &sb) != 0) {
#ifdef _WIN32
		if(CreateDirectory(log_path, nullptr) == 0) {
#else
		if(mkdir(log_path,00777) == -1) {
#endif
			cerr << "Failed to create data folder search" << log_path << endl;
#ifdef _WIN32
			cerr << "Error code: " << GetLastError() << endl;
#endif
			exit(1);
		}
	}

	sprintf(dir,"%s/%d",log_path,static_cast<int>(timer));
#ifdef _WIN32
	if(CreateDirectory(dir, nullptr) == 0) {
#else
	if(mkdir(dir,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << dir << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}

class TreeQueryTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {
		strcpy(setting_path, setting);
	}

	static void TearDownTestCase() {
		::delete hk1;
		hk1 = nullptr;
		::delete query;
		query = nullptr;
	}

	virtual void SetUp() {}
	virtual void TearDown() {}

public:
	static char
	setting_path[256],
	base_path[256],
	query_path[256],
	cq_path[256],
	pq_path[256],
	code_path[256],
	log_path[256];
	static TreeQuery * hk1;
	static float * query;
	static int N, h, l, H, L, T, M;
};

char TreeQueryTest::base_path[256];
char TreeQueryTest::setting_path[256];
char TreeQueryTest::query_path[256];
char TreeQueryTest::cq_path[256];
char TreeQueryTest::pq_path[256];
char TreeQueryTest::code_path[256];
char TreeQueryTest::log_path[256];
TreeQuery * TreeQueryTest::hk1;
float * TreeQueryTest::query;
int TreeQueryTest::N;
int TreeQueryTest::h;
int TreeQueryTest::l;
int TreeQueryTest::H;
int TreeQueryTest::L;
int TreeQueryTest::T;
int TreeQueryTest::M;

/**
 * Load configuration
 */
TEST_F(TreeQueryTest, test0) {
	ifstream input;
	input.open(setting_path, ios::in);
	input >> base_path >> query_path
	>> cq_path >> pq_path >> code_path
	>> log_path
	>> h >> l >> H >> L >> T >> M;
	cout << "h=" << h << ";l=" << l << ";T=" << T << endl;
	input.close();
}

TEST_F(TreeQueryTest, test1) {
	hk1 = new TreeQuery();
}

TEST_F(TreeQueryTest, test2) {
	hk1->load_codebooks(
			cq_path,
			pq_path,
			true);
}

TEST_F(TreeQueryTest, test3) {
	hk1->load_encoded_data(
			code_path,
			false);
}

TEST_F(TreeQueryTest, test4) {
	hk1->pre_compute1();
}

TEST_F(TreeQueryTest, test5) {
	N = load_and_convert_data<unsigned char,float>(query_path,query,4,128,true);
	M = std::min(M,N);
	//hk1->load_data<unsigned char>(base_path,4,true);
}

TEST_F(TreeQueryTest, DISABLED_test6) {
	clock_t st, ed;
	time_t timer = time(NULL);
	create_log_dir(log_path,timer);
	double t = 0.0;
	double t0, t1, t2;
	int R, r[] = {1,10,100,1000,10000};
	ofstream output;
	char filename[256];
	int * id_h, * id_l, * ids, * idf;
	float * q, * dist_h, * dist_l, * dists, * distf;

	SimpleCluster::init_array(dist_h,H);
	SimpleCluster::init_array(dist_l,L);
	SimpleCluster::init_array(id_h,H);
	SimpleCluster::init_array(id_l,L);
	SimpleCluster::init_array(ids,H * L);
	SimpleCluster::init_array(dists,H * L);

	int sum = 0;
	for(int i = 0; i < 5; i++) {
		t = 0.0;
		t0 = t1 = t2 = 0.0;
		R = r[i];
		q = query;
		sprintf(filename,"%s/%d/search_result_%d.txt",
				log_path,
				static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < M; j++) {
			st = clock();
			hk1->search(
					q,
					T, h, l, R,
					dist_h, id_h,
					dist_l, id_l,
					dists, ids,
					distf, idf,
					sum, t0, t1, t2,
					false,
					false);
			ed = clock();
			t += static_cast<double>(ed - st);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << idf[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
			output << endl;
			q += 128;
			::delete distf;
			distf = nullptr;
			::delete idf;
			idf =  nullptr;
		}

		output.close();
		cout << "Finished search@" << R << " in " <<
				1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Precomputing time: " << 1000.0f * t0 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Sorting time: " << 1000.0f * t1 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Linear Searching time: " << 1000.0f * t2 / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

TEST_F(TreeQueryTest, test7) {
	hk1->pre_compute5();
}

TEST_F(TreeQueryTest, test8) {
	clock_t st, ed;
	time_t timer = time(NULL);
	create_log_dir(log_path,timer);
	double t = 0.0;
	double t0, t1, t2;
	int R, r[] = {1,10,100,1000,10000};
	ofstream output;
	char filename[256];
	int * id_h, * id_l, * ids, * idf;
	float * q, * dist_h, * dist_l, * dists, * distf;

	SimpleCluster::init_array(dist_h,H);
	SimpleCluster::init_array(dist_l,L);
	SimpleCluster::init_array(id_h,H);
	SimpleCluster::init_array(id_l,L);
	SimpleCluster::init_array(ids,H * L);
	SimpleCluster::init_array(dists,H * L);

	int sum = 0;
	for(int i = 0; i < 5; i++) {
		t = 0.0;
		t0 = t1 = t2 = 0.0;
		R = r[i];
		q = query;
		sprintf(filename,"%s/%d/search_result_%d.txt",
				log_path,
				static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < M; j++) {
			st = clock();
			hk1->search_mem(
					q,
					T, h, l, R,
					dist_h, id_h,
					dist_l, id_l,
					dists, ids,
					distf, idf,
					sum, t0, t1, t2,
					false,
					false);
			ed = clock();
			t += static_cast<double>(ed - st);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << idf[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
			output << endl;
			q += 128;
			::delete distf;
			distf = nullptr;
			::delete idf;
			idf =  nullptr;
		}

		output.close();
		cout << "Finished search@" << R << " in " <<
				1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Precomputing time: " << 1000.0f * t0 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Sorting time: " << 1000.0f * t1 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Linear Searching time: " << 1000.0f * t2 / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

TEST_F(TreeQueryTest, DISABLED_test9) {
	clock_t st, ed;
	time_t timer = time(NULL);
	create_log_dir(log_path,timer);
	double t = 0.0;
	double t0, t1, t2;
	int R, r[] = {1,10,100,1000,10000};
	ofstream output;
	char filename[256];
	int * id_h, * id_l, * ids, * idf;
	float * q, * dist_h, * dist_l, * dists, * distf;

	SimpleCluster::init_array(dist_h,H);
	SimpleCluster::init_array(dist_l,L);
	SimpleCluster::init_array(id_h,H);
	SimpleCluster::init_array(id_l,L);
	SimpleCluster::init_array(ids,H * L);
	SimpleCluster::init_array(dists,H * L);

	int sum = 0;
	for(int i = 0; i < 5; i++) {
		t = 0.0;
		t0 = t1 = t2 = 0.0;
		R = r[i];
		q = query;
		sprintf(filename,"%s/%d/search_result_%d.txt",
				log_path,
				static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < M; j++) {
			st = clock();
			hk1->search(
					q,
					T, h, l, R,
					dist_h, id_h,
					dist_l, id_l,
					dists, ids,
					distf, idf,
					sum, t0, t1, t2,
					true,
					false);
			ed = clock();
			t += static_cast<double>(ed - st);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << idf[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
			output << endl;
			q += 128;
			::delete distf;
			distf = nullptr;
			::delete idf;
			idf =  nullptr;
		}

		output.close();
		cout << "Finished search@" << R << " in " <<
				1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Precomputing time: " << 1000.0f * t0 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Sorting time: " << 1000.0f * t1 / CLOCKS_PER_SEC << "[ms]" << endl;
		cout << "Linear Searching time: " << 1000.0f * t2 / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cerr << "ERROR: Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " /path/to/the/config/file" << endl;
		exit(EXIT_FAILURE);
	}
	// Parse the arguments
	snprintf(setting,sizeof(setting),"%s",argv[1]);

	/**
	 * RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	 * It's must be called only once in the code because multiple calls lead to conflicts and,
	 * therefore, are not supported.
	 **/
	return RUN_ALL_TESTS();
}

