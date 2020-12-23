cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(scf_plugin CXX)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(scf_plugin
		tensors.cc
        scf_plugin.cc
		main.cc
		)
