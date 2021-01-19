cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(rhf_mp2 CXX)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(rhf_mp2
		tensors.cc
        rhf_mp2.cc
		main.cc
		)
