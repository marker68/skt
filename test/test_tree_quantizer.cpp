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
 *  test_hk_quantizer.cpp
 *
 *  Created on: 2015/04/01
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

#include "tree_quantizer.h"

using namespace std;
using namespace PQLearn;

char config_file[256];

class TreeQuantizerTest : public ::testing::Test {
protected:
	static void SetUpTestCase() {
		ifstream input;
		input.open(config_file, ios::in);
		input >> data_in >> cq_in >> pq_out >> tree_cq_out
		>> d >> m >> h >> l >> mp >> kp >> offset;
		input.close();
	}

	static void TearDownTestCase() {

	}

	virtual void SetUp() {}
	virtual void TearDown() {}

public:
	static char data_in[256], cq_in[256], pq_out[256], tree_cq_out[256];
	static TreeQuantizer<float> * hk1;
	static int d, m, h, l, mp, kp, offset;
};

TreeQuantizer<float> * TreeQuantizerTest::hk1;
int TreeQuantizerTest::d;
int TreeQuantizerTest::m;
int TreeQuantizerTest::h;
int TreeQuantizerTest::l;
int TreeQuantizerTest::mp;
int TreeQuantizerTest::kp;
int TreeQuantizerTest::offset;
char TreeQuantizerTest::data_in[256];
char TreeQuantizerTest::cq_in[256];
char TreeQuantizerTest::pq_out[256];
char TreeQuantizerTest::tree_cq_out[256];

TEST_F(TreeQuantizerTest, test1) {
	hk1 = new TreeQuantizer<float>(d,m,h,l,offset,true);
	hk1->set_type(1);
}

TEST_F(TreeQuantizerTest, test2) {
	hk1->load_data(data_in,true);
}

TEST_F(TreeQuantizerTest, DISABLED_test3) {
	hk1->create_sub_quantizers(true);
}

TEST_F(TreeQuantizerTest, DISABLED_test4) {
	cout << "Distortion is " << hk1->distortion(false) << endl;
}

TEST_F(TreeQuantizerTest, DISABLED_test5) {
	hk1->create_low_level_cells(false);
}

TEST_F(TreeQuantizerTest, DISABLED_test6) {
	hk1->tree_output(tree_cq_out,true);
}

TEST_F(TreeQuantizerTest, DISABLED_test7) {
	hk1->set_params(0,8,256);
	hk1->create_sub_quantizers(false);
}

TEST_F(TreeQuantizerTest, DISABLED_test8) {
	hk1->output(pq_out,true);
}

TEST_F(TreeQuantizerTest, test9) {
	hk1->reload_high_level_cells(cq_in,true);
}

TEST_F(TreeQuantizerTest, test10) {
	hk1->create_low_level_cells(true);
}

TEST_F(TreeQuantizerTest, test11) {
	hk1->tree_output(tree_cq_out,true);
}

TEST_F(TreeQuantizerTest, test12) {
	hk1->set_params(0,mp,kp);
	hk1->create_sub_quantizers(true);
}

TEST_F(TreeQuantizerTest, test13) {
	hk1->output(pq_out,true);
}

int main(int argc, char * argv[]) {
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	if(argc != 2) {
		cerr << "ERROR: Not enough arguments" << endl;
		cerr << "Usage: " << argv[0] << " /path/to/the/config/file" << endl;
		exit(EXIT_FAILURE);
	}
	// Parse the arguments
	snprintf(config_file,sizeof(config_file), "%s",argv[1]);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
