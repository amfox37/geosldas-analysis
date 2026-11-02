#!/bin/bash
#SBATCH --job-name=build_jpss_5d4d636
#SBATCH --account=g0610
#SBATCH --qos=debug
#SBATCH --time=0:45:00
#SBATCH --ntasks=24
#SBATCH --constraint=mil
#SBATCH --output=/gpfsm/dnb06/projects/p284/M21C_testing/build_jpss_5d4d636.%j.log
# Incremental build of feature/amfox/jpss_hdf5_reader at 5d4d636 (ca6b377 + cap VIIRS SCF obs at 1)
# reusing build-jpss-ca6b377 (only clsm_ensupd_read_obs.F90 changed), installed into a new install-jpss-5d4d636.
set -e
G=/discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas
restore() {
  git -C $G/src/Components/@GEOSldas_GridComp checkout -q feature/amfox/matlab-obsfcstana-nc4
  echo RESTORED
}
trap restore EXIT
source $G/@env/g5_modules.sh
name=jpss-5d4d636
git -C $G/src/Components/@GEOSldas_GridComp checkout -q 5d4d636
git -C $G/src/Components/@GEOSldas_GridComp log --oneline -1
cd $G/build-jpss-ca6b377
cmake . -DCMAKE_INSTALL_PREFIX=$G/install-$name
make -j 24 install
ls -l $G/install-$name/bin/GEOSldas.x $G/install-$name/bin/ldas_setup $G/install-$name/bin/preprocess_ldas.x
