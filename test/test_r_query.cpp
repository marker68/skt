/*
 * test_r_query.cpp
 *
 *  Created on: 2015/01/03
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include "r_query.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <gtest/gtest.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

using namespace std;
using namespace PQLearn;

/**
 * Customized test case for testing
 */
class QueryTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		strcpy(setting_path, "./data/config/test_r_query.txt");
		sprintf(db_prefix,"test");
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		::delete worker;
		worker = nullptr;
		::delete data;
		data = nullptr;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static char
	setting_path[256],
	base_dir[256],
	query_path[256],
	db_prefix[256];
	static PQRQuery * worker;
	static float * data;
	static int d, w, T, N, kc;
};

// Global variables
char QueryTest::base_dir[256];
char QueryTest::setting_path[256];
char QueryTest::query_path[256];
char QueryTest::db_prefix[256];
PQRQuery * QueryTest::worker;
float * QueryTest::data;
int QueryTest::d;
int QueryTest::w;
int QueryTest::T;
int QueryTest::N;
int QueryTest::kc;

/**
 * Load configuration
 */
TEST_F(QueryTest, test0) {
	ifstream input;
	input.open(setting_path, ios::in);
	input >> base_dir >> query_path
	>> d >> w >> T >> kc;
	input.close();
}

TEST_F(QueryTest, test1) {
	worker = new PQRQuery();
}

TEST_F(QueryTest, test2) {
	worker->load_codebooks(
			"./data/codebooks/sift_test1_cq.ctr_",
			"./data/codebooks/sift_test1_4_pq.ctr_",
			"./data/codebooks/sift_test1_4_rq.ctr_",
			true);
	worker->load_encoded_data("./data/codebooks/sift_test1_8_r_ivf.edat_",false);
	cout << "Entropy: " << worker->entropy(worker->get_full_size()) << endl;
	cout << "Max entropy:" << log2(worker->get_size()) << endl;
}

TEST_F(QueryTest, test3) {
	N = load_data<float>(query_path,data,4,d,true);
	worker->load_data<unsigned char>("./data/sift/sift_base.bvecs",4,true);
	worker->pre_compute1();
}

TEST_F(QueryTest, test4) {
	clock_t st, ed;
	time_t timer = time(NULL);
	char data_folder[256];
	sprintf(data_folder,"%s/search",base_dir);
	struct stat sb;
if(stat(data_folder, &sb) != 0) {
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder search" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}
sprintf(data_folder,"%s/search/%d",base_dir,static_cast<int>(timer));
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
	double t = 0.0;
	int R, r[] = {1,10,100,1000,10000,100000};
	ofstream output;
	char filename[256];
	int * result, * buckets, * i_tmp, * bk;
	float * v_tmp, * dist, * tmp;
	SimpleCluster::init_array(v_tmp,kc);
	SimpleCluster::init_array(buckets,kc);
	SimpleCluster::init_array(bk,kc);

	int sum = 0;
	for(int i = 0; i < 5; i++) {
		t = 0.0;
		R = r[i];
		tmp = data;
		sprintf(filename,"%s/search/%d/search_result_%d.txt",base_dir,static_cast<int>(timer),R);
		output.open(filename,ios::out);
		for(int j = 0; j < N; j++) {
			st = clock();
			worker->search_ivfadc(
					tmp,
					v_tmp,dist,
					result,buckets,bk,
					sum,R,w,T,false,false);
			ed = clock();
			t += static_cast<double>(ed - st);
			for(int i = 0; i < (R>sum?sum:R); i++)
				output << result[i] << " ";
			if(sum < R)
				for(int i = 0; i < R - sum; i++)
					output << "-1 ";
			output << endl;
			tmp += d;
			::delete result;
			result = nullptr;
			::delete dist;
			dist =  nullptr;
		}

		output.close();
		cout << "Finished search@" << R << " in " <<
				1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
	}
}

TEST_F(QueryTest, DISABLED_test5) {
	clock_t st, ed;
	time_t timer = time(NULL);
	char data_folder[256];
	sprintf(data_folder,"%s/search",base_dir);
	struct stat sb;
if(stat(data_folder, &sb) != 0) {
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder search" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
}
sprintf(data_folder,"%s/search/%d",base_dir,static_cast<int>(timer));
#ifdef _WIN32
if(CreateDirectory(data_folder, nullptr) == 0) {
#else
	if(mkdir(data_folder,00777) == -1) {
#endif
		cerr << "Failed to create data folder" << endl;
#ifdef _WIN32
		cerr << "Error code: " << GetLastError() << endl;
#endif
		exit(1);
	}
	double t = 0.0;
	ofstream output;
	char filename[256];
	int * result, * buckets, * i_tmp, * bk;
	float * v_tmp, * dist, * tmp;
	SimpleCluster::init_array(v_tmp,kc);
	SimpleCluster::init_array(buckets,kc);
	SimpleCluster::init_array(bk,kc);

	int sum = 0;
	t = 0.0;
	tmp = data;
	sprintf(filename,"%s/search/%d/search_result_%d.txt",base_dir,static_cast<int>(timer),1);
	output.open(filename,ios::out);
	for(int j = 0; j < N; j++) {
		st = clock();
		worker->search_ivfadc(
				tmp,
				v_tmp,dist,
				result,buckets,bk,
				sum,1,w,T,false,false);
		ed = clock();
		t += static_cast<double>(ed - st);
		for(int i = 0; i < 1; i++)
			output << result[i] << " ";
		output << endl;
		tmp += d;
		::delete result;
		result = nullptr;
		::delete dist;
		dist =  nullptr;
	}
	output.close();
	cout << "Finished search@1 in " <<
			1000.0f * t / CLOCKS_PER_SEC << "[ms]" << endl;
}

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
It's must be called only once in the code because multiple calls lead to conflicts and,
therefore, are not supported.
	 */
	return RUN_ALL_TESTS();
}

