#!/bin/bash
#SBATCH --job-name=build_jpss_ca6b377
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:45:00
#SBATCH --ntasks=24
#SBATCH --constraint=mil
#SBATCH --output=/gpfsm/dnb06/projects/p284/M21C_testing/build_jpss_ca6b377.%j.log
# Builds feature/amfox/jpss_hdf5_reader at ca6b377 (cee4997 + VIIRS SCF QC: keep Basic_QA 0,1,3; Clear_Index<=100; comments)
# from the shared land_da checkout into new build/install dirs; earlier installs stay untouched.
set -e
G=/discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas
restore() {
  git -C $G/src/Components/@GEOSldas_GridComp checkout -q feature/amfox/matlab-obsfcstana-nc4
  echo RESTORED
}
trap restore EXIT
source $G/@env/g5_modules.sh
git -C $G/src/Components/@GEOSgcm_GridComp log --oneline -1
name=jpss-ca6b377
git -C $G/src/Components/@GEOSldas_GridComp checkout -q ca6b377
git -C $G/src/Components/@GEOSldas_GridComp log --oneline -1
mkdir -p $G/build-$name && cd $G/build-$name
cmake .. -DBASEDIR=$BASEDIR/Linux -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_BUILD_TYPE=Aggressive \
  -DCMAKE_INSTALL_PREFIX=$G/install-$name
make -j 24 install
ls -l $G/install-$name/bin/GEOSldas.x $G/install-$name/bin/ldas_setup $G/install-$name/bin/preprocess_ldas.x
