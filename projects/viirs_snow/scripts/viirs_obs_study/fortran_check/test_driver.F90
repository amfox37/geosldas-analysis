! Run the GEOSldas SCF hdf readers (read_JPSS_SCF_hdf, read_MODIS_SCF_hdf) on one
! daily CMG file for the whole globe and write the kept obs as raw real*4:
!   N_good (int32), then lon(N_good), lat(N_good), scf(N_good) (real32), in reader order
!
! usage: test_driver VIIRS|MODIS <cmg_file> <out_file>

program test_driver
  use test_mod
  use, intrinsic :: iso_fortran_env, only: int32, real32
  implicit none
  character(400) :: fname, sensor, fout
  real, parameter :: lon_min = -180., lon_max = 179.999, lat_min = -89.99, lat_max = 90.
  integer :: N_lon, N_lat, N_good, start(2), last(2), iu
  real, allocatable :: lon(:), lat(:), scf(:)

  call get_command_argument(1, sensor)
  call get_command_argument(2, fname)
  call get_command_argument(3, fout)

  ! same sub-array logic as read_obs_*_SCF (whole globe here)
  start(1) = (lon_min + 180.)/0.05 ; last(1) = (lon_max + 180.)/0.05
  start(2) = (90. - lat_max)/0.05  ; last(2) = (90. - lat_min)/0.05
  N_lon = last(1)-start(1)+1 ; N_lat = last(2)-start(2)+1
  allocate(lon(N_lon*N_lat), lat(N_lon*N_lat), scf(N_lon*N_lat))

  if     (trim(sensor) == 'VIIRS') then
     call read_JPSS_SCF_hdf( trim(fname), lon_min, lon_max, lat_min, lat_max, N_good, lon, lat, scf)
  elseif (trim(sensor) == 'MODIS') then
     call read_MODIS_SCF_hdf(trim(fname), lon_min, lon_max, lat_min, lat_max, N_good, lon, lat, scf)
  else
     stop 'first argument must be VIIRS or MODIS'
  end if

  open(newunit=iu, file=trim(fout), access='stream', form='unformatted', status='replace')
  write(iu) int(N_good, int32)
  write(iu) real(lon(1:N_good), real32), real(lat(1:N_good), real32), real(scf(1:N_good), real32)
  close(iu)

end program test_driver
