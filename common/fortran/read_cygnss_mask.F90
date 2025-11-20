program read_cygnss_mask
    
    use netcdf

    implicit none

! Declare the dimids and dim_len arrays if not already declared
integer :: ndims, ierr
integer :: ncid, small_SM_range_varid, poor_SMAP_varid, high_ubrmsd_varid, few_obs_varid, low_signal_varid
integer :: lat_dimid, lon_dimid
integer :: N_lat_m, N_lon_m
real, dimension(:,:), allocatable :: latitudes_m, longitudes_m
integer, dimension(:,:), allocatable :: small_SM_range, poor_SMAP, high_ubrmsd, few_obs, low_signal

! Open the NetCDF file
ierr = nf90_open('/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/cygnss_test/ucar_cu_cygnss_sm_v1_static_flags.nc', nf90_nowrite, ncid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error opening file:', ierr
    stop
end if

! Inquire variable IDs
ierr = nf90_inq_varid(ncid, 'latitude', lat_dimid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring latitude variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'longitude', lon_dimid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring longitude variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'flag_small_SM_range', small_SM_range_varid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring flag_small_SM_range variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'flag_poor_SMAP', poor_SMAP_varid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring flag_poor_SMAP variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'flag_high_ubrmsd', high_ubrmsd_varid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring flag_high_ubrmsd variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'flag_few_obs', few_obs_varid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring flag_few_obs variable ID:', ierr
    stop
end if

ierr = nf90_inq_varid(ncid, 'flag_low_signal', low_signal_varid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error inquiring flag_low_signal variable ID:', ierr
    stop
end if

    ! dimensions sizes
ierr = nf90_inquire_dimension(ncid, lon_dimid, len=N_lon_m)
ierr = nf90_inquire_dimension(ncid, lat_dimid, len=N_lat_m)

! Ensure the arrays are allocated with the correct dimensions
allocate(latitudes_m(N_lon_m, N_lat_m))
allocate(longitudes_m(N_lon_m, N_lat_m))
allocate(small_SM_range(N_lon_m, N_lat_m))
allocate(poor_SMAP(N_lon_m, N_lat_m))
allocate(high_ubrmsd(N_lon_m, N_lat_m))
allocate(few_obs(N_lon_m, N_lat_m))
allocate(low_signal(N_lon_m, N_lat_m))

! Read the variables
ierr = nf90_get_var(ncid, lat_dimid, latitudes_m)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading latitude:', ierr
    stop
end if

ierr = nf90_get_var(ncid, lon_dimid, longitudes_m)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading longitude:', ierr
    stop
end if

ierr = nf90_get_var(ncid, small_SM_range_varid, small_SM_range)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading flag_small_SM_range:', ierr
    stop
end if

ierr = nf90_get_var(ncid, poor_SMAP_varid, poor_SMAP)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading flag_poor_SMAP:', ierr
    stop
end if

ierr = nf90_get_var(ncid, high_ubrmsd_varid, high_ubrmsd)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading flag_high_ubrmsd:', ierr
    stop
end if

ierr = nf90_get_var(ncid, few_obs_varid, few_obs)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading flag_few_obs:', ierr
    stop
end if

ierr = nf90_get_var(ncid, low_signal_varid, low_signal)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error reading flag_low_signal:', ierr
    stop
end if

! Close the NetCDF file
ierr = nf90_close(ncid)
if (ierr /= nf90_noerr) then
    write(*,*) 'Error closing file:', ierr
    stop
end if

! Print the first 10 values of each variable
write(*,*) 'latitudes_m:', latitudes_m(1:10, 1)
write(*,*) 'longitudes_m:', longitudes_m(1:10, 1)
write(*,*) 'small_SM_range:', small_SM_range(1:10, 1)
write(*,*) 'poor_SMAP:', poor_SMAP(1:10, 1)
write(*,*) 'high_ubrmsd:', high_ubrmsd(1:10, 1)
write(*,*) 'few_obs:', few_obs(1:10, 1)
write(*,*) 'low_signal:', low_signal(1:10, 1)

end program read_cygnss_mask