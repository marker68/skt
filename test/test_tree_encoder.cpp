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
 *  test_hk_encoder.cpp
 *
 *  Created on: 2015/04/05
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

#include "tree_encoder.h"

using namespace std;
using namespace PQLearn;

char config_file[256];

class TreeEncoderTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {
		ifstream input;
		input.open(config_file, ios::in);
		input >> data_in >> tree_cq_in
		>> pq_in >> code_in >> code_out >> code_out_r
		>> m >> h >> l >> offset;
		input.close();
	}
	static void TearDownTestCase() {
		::delete hk1;
		hk1 = nullptr;
		::delete hk2;
		hk2 = nullptr;
	}

	virtual void SetUp() {}
	virtual void TearDown() {}

public:
	static TreeEncoder * hk1, * hk2;
	static char tree_cq_in[256], pq_in[256],
	data_in[256], code_out[256], code_out_r[256], code_in[256];
	static int m, h, l, offset;
};


TreeEncoder * TreeEncoderTest::hk1;
TreeEncoder * TreeEncoderTest::hk2;
char TreeEncoderTest::data_in[256];
char TreeEncoderTest::tree_cq_in[256];
char TreeEncoderTest::pq_in[256];
char TreeEncoderTest::code_in[256];
char TreeEncoderTest::code_out[256];
char TreeEncoderTest::code_out_r[256];
int TreeEncoderTest::m;
int TreeEncoderTest::h;
int TreeEncoderTest::l;
int TreeEncoderTest::offset;


TEST_F(TreeEncoderTest, test1) {
	hk1 = new TreeEncoder();
	hk2 = new TreeEncoder();
}

TEST_F(TreeEncoderTest, test2) {
	hk1->load_codebooks(
			tree_cq_in,
			pq_in,
			true);
	hk2->load_codebooks(
			tree_cq_in,
			pq_in,
			true);
}

TEST_F(TreeEncoderTest, test3) {
	hk1->tree_encode<float>(
			data_in,
			offset,
			true);
}

TEST_F(TreeEncoderTest, test4) {
	hk1->distribution(false);
}

TEST_F(TreeEncoderTest, test5) {
	hk1->output(code_out,false);
}

//TEST_F(TreeEncoderTest, test6) {
//	hk2->set_mp(m);
//	hk2->set_size(h);
//	hk2->load_encoded_data(code_in,false);
//}
//
//TEST_F(TreeEncoderTest, test7) {
//	hk2->tree_rencode<unsigned char>(
//			data_in,
//			4,
//			true);
//}
//
//TEST_F(TreeEncoderTest, test8) {
//	hk2->set_size(h * l);
//	hk2->distribution(true);
//}
//
//TEST_F(TreeEncoderTest, test9) {
//	hk2->output(code_out_r,false);
//}

int main(int argc, char * argv[]) {
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cerr << "Usage: " << argv[0] << " /path/to/the/config/file" << endl;
		exit(EXIT_FAILURE);
	}
	snprintf(config_file,sizeof(config_file), "%s",argv[1]);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
