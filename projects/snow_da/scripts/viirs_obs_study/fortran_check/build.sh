#!/bin/bash
# Build test_driver from the SCF reader routines of a GEOSldas_GridComp checkout.
#
# usage: build.sh <GEOSlandassim_GridComp dir> <build dir>
#   e.g. build.sh ~/Desktop/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp /tmp/scf_check
#
# Needs the Mac Baselibs (HDF4, HDF5) and Homebrew gfortran/mpifort/clang; see README.md.
set -euo pipefail
G=$(cd "$1" && pwd)
B=$2
HERE=$(cd "$(dirname "$0")" && pwd)
P=${BASELIBS:-/Users/amfox/Desktop/GEOSldas_develop/aarch64-apple-darwin24.5.0/gfortran/Darwin}

F="gfortran -O0 -g -fcheck=all -ffree-line-length-none -I$P/include/hdf5"
L="-L$P/lib -lmfhdf -ldf -ljpeg -lhdf5_fortran -lhdf5 -lz -lsz -laec -lcurl -lm"

mkdir -p "$B"; cd "$B"

# test module = header with stubs + the two reader routines copied verbatim from the checkout
{ cat "$HERE/test_mod_header.F90"
  awk '/^  subroutine read_JPSS_SCF_hdf\(/{p=1}  p{print} /^  end subroutine read_JPSS_SCF_hdf/{p=0}'  "$G/clsm_ensupd_read_obs.F90"
  awk '/^  subroutine read_MODIS_SCF_hdf\(/{p=1} p{print} /^  end subroutine read_MODIS_SCF_hdf/{p=0}' "$G/clsm_ensupd_read_obs.F90"
  echo "end module test_mod"; } > test_mod.F90

grep -q "subroutine read_JPSS_SCF_hdf"  test_mod.F90 || { echo "read_JPSS_SCF_hdf not found in $G"; exit 1; }
grep -q "subroutine read_MODIS_SCF_hdf" test_mod.F90 || { echo "read_MODIS_SCF_hdf not found in $G"; exit 1; }

$F -c "$G/io_hdf5.F90"          -o io_hdf5.o
$F -c "$G/hdf4_fortran_api.F90" -o hdf4_fortran_api.o
clang -c "$G/hdf4_c_api.c" -I"$P/include/hdf" -o hdf4_c_api.o
$F -c test_mod.F90              -o test_mod.o
$F -c "$HERE/test_driver.F90"   -o test_driver.o
/opt/homebrew/bin/mpifort test_driver.o test_mod.o io_hdf5.o hdf4_fortran_api.o hdf4_c_api.o $L -o test_driver 2>&1 | grep -v "ld: warning" || true

# record which code was built
git -C "$G" log -1 --format='GridComp %h %s' > built_from.txt
cat built_from.txt
ls -l test_driver
