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
 *  hk_query.cpp
 *
 *  Created on: 2015/04/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */


#include <iostream>

#include "tree_query.h"
#include "pq_utilities.h"

using namespace std;

namespace PQLearn {

TreeQuery::TreeQuery() : PQQuery::PQQuery() {
	high_centers =  nullptr;
	low_centers =  nullptr;
	norm_ch = nullptr;
	diff_qch = nullptr;
	high = low = 0;
}

TreeQuery::~TreeQuery() {
//	::delete high_centers;
//	high_centers = nullptr;
//	::delete low_centers;
//	low_centers = nullptr;
}


void TreeQuery::load_codebooks(
		const char * cq_path,
		const char * pq_path,
		bool verbose) {
	load_codebook<float>(cq_path,this->config,this->cq,3,verbose);
	high_centers = this->cq;
	low_centers = this->cq + this->config.kc * this->config.dim;
	high = this->config.kc;
	low = this->config.L;
	load_codebook<float>(pq_path,this->config,this->pq,1,verbose);
	this->size = static_cast<int>(high * low);
	if(verbose)
		cout << "The size of ivf:" << this->size << endl;
	cout << "--> Settings: (kc,mc,L,kp,mp)=" << high << " "
			<< this->config.mc << " " << low << " " << this->config.kp << " "
			<< this->config.mp << endl;
	cout << high << endl;
}


} /* namespace PQLearn */
