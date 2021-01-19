cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(mp_plugin CXX)

find_package(psi4 1.1 REQUIRED)

add_psi4_plugin(mp_plugin
		tensors.cc
        mp_plugin.cc
		main.cc
		)
