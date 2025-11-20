program read_netcdf
    use netcdf
    implicit none
    integer :: ncid, varid, ierr
    integer, dimension(:), allocatable :: asc_lon, asc_lat, cold_mask, wet_mask, veg_mask, subsurface_mask

    ! Open the NetCDF file
    ierr = nf90_open('../test_data/clsm/subsurface_scattering_ASCAT_ERA5_Land.nc', nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) stop 'Error opening file'

    ! Get the variable IDs and read the variables
    ierr = nf90_inq_varid(ncid, 'asc_lon', varid)
    ierr = nf90_get_var(ncid, varid, asc_lon)
    ierr = nf90_inq_varid(ncid, 'asc_lat', varid)
    ierr = nf90_get_var(ncid, varid, asc_lat)
    ierr = nf90_inq_varid(ncid, 'cold_mask', varid)
    ierr = nf90_get_var(ncid, varid, cold_mask)
    ierr = nf90_inq_varid(ncid, 'wet_mask', varid)
    ierr = nf90_get_var(ncid, varid, wet_mask)
    ierr = nf90_inq_varid(ncid, 'veg_mask', varid)
    ierr = nf90_get_var(ncid, varid, veg_mask)
    ierr = nf90_inq_varid(ncid, 'subsurface_mask', varid)
    ierr = nf90_get_var(ncid, varid, subsurface_mask)

    ! Close the NetCDF file
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) stop 'Error closing file'
end program read_netcdf