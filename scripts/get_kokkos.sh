#!/bin/bash

# when changing dependency versions,
# make sure that the related paths in the Makefile are still valid.
kokkos_ver="4.2.01"

get_kokkos() {
	local ver="$1"
	# clone KOKKOS, if no kokkos dir exists
	if [ -d "kokkos" ]; then
		echo "Existing directory \"kokkos\" found. Skipping download..."
	else
		echo "Cloning kokkos..."
		git clone --depth 1 --branch $ver https://github.com/kokkos/kokkos.git
	fi
}

get_kokkoskernels() {
	local ver="$1"
	# clone KOKKOS Kernels, if no kokkos dir exists
	if [ -d "kokkos-kernels" ]; then
		echo "Existing directory \"kokkos-kernels\" found. Skipping download..."
	else
		echo "Cloning kokkos-kernels..."
		git clone --depth 1 --branch $ver https://github.com/kokkos/kokkos-kernels.git
	fi
}

get_kokkos "$kokkos_ver"
get_kokkoskernels "$kokkos_ver"
