Selective K-means Tree
================================

## Introduction

This project contains the source code of the paper "**Selective K-means Tree Search**".
In this project, we achieved to create a library for nearest neighbors search with the following features:

* An implementation of [Product Quantization for Approximated Nearest Neighbor Search](https://hal.inria.fr/inria-00514462/document)
* An implementation of [The Inverted Multi-index](http://sites.skoltech.ru/app/data/uploads/sites/25/2014/12/TPAMI14.pdf)
* An implementation of [Re-ranking](https://hal.inria.fr/inria-00566883/PS/paper.ps)
* The [Selective K-means Tree Search](http://dl.acm.org/citation.cfm?id=2806353) which apply hierarchical k-means to nearest neighbor search.

## Installation

#### Prerequisites

* An Unix Operating System. Tested on Mac OS X 10.9 and Ubuntu 14.04.
* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An GNU C++ Compiler or LLVM compiler that support C++ 11 and OpenMP.
* [BOOST](http://www.boost.org) library installed on the machine.

#### Build

We use CMake as the build system. On terminal,
```bash
$ git clone git@github.com:marker68/skt.git
$ cd ./skt
$ git submodule update --init
$ cmake -H. -Bbuild && cmake --build build -- -j4
# OR IF YOU WANT TO BUILD TESTS
$ cmake -DBUILD_TEST=ON -H. -Bbuild && cmake --build build -- -j4
```
This script will create binaries in your `bin/` and `lib/` directories. 

## Documentation

We are using Doxygen to generate the documentation. You can find that the style of comments for every methods that are implemented in these source codes are matching with Doxygen style.
To generate the documentation:

```bash
$ cd ./doc
$ doxygen config.doxygen
```

Doxygen will create the documentation files in HTML and LaTeX format. You can also find samples and documentation in the Wiki of this project.

## To developers

* This project was built under Eclipse CDT so we leave the configuration files `.cproject` and `.project` of the project here. If you are using Eclipse, just import the folder that contains the source code.

## Related Publications
If you use this source code, please also cite the following publication:
```
@inproceedings{Nguyen2015,
 author = {Nguyen, Tuan Anh and Matsui, Yusuke and Yamasaki, Toshihiko and Aizawa, Kiyoharu},
 title = {Selective K-means Tree Search},
 booktitle = {Proceedings of the 23rd Annual ACM Conference on Multimedia Conference},
 year = {2015},
 pages = {875--878},
 publisher = {ACM},
} 
```

