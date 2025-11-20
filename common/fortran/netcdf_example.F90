! example.f90
program main
    use :: netcdf
    implicit none
    character(len=*), parameter :: FILE_NAME = 'data.nc' ! Export/import file.
    integer,          parameter :: NX        = 12        ! Number of columns.
    integer,          parameter :: NY        = 6         ! Number of rows.

    integer :: array(NX, NY)
    integer :: x, y

    ! Create and output sample data.
    print '("Data to be written to NetCDF file:")'
    print '(a)', repeat('-', 64)

    do y = 1, NY
        do x = 1, NX
            array(x, y) = (y - 1) * NY + (x - 1)
            write (*, '(i4)', advance='no') array(x, y)
        end do
        write (*, *)
    end do

    print '(a)', repeat('-', 64)

    ! Write data to NetCDF file (into variable `data`).
    call export_int_array(FILE_NAME, array, 'data')

    array = 0

    ! Read data back in from the same NetCDF file (from variable `data`).
    call import_int_array(FILE_NAME, array, 'data')

    ! Print read data to stdout.
    print '(/, "Data read from NetCDF file:")'
    print '(a)', repeat('-', 64)

    do y = 1, NY
        do x = 1, NX
            write (*, '(i4)', advance='no') array(x, y)
        end do
        write (*, *)
    end do

    print '(a)', repeat('-', 64)
contains
    subroutine check(stat)
        integer, intent(in) :: stat

        if (stat /= NF90_NOERR) then
            print '(a)', trim(nf90_strerror(stat))
            stop
        end if
    end subroutine check

    subroutine export_int_array(file_name, array, var_name)
        character(len=*), intent(in)    :: file_name
        integer,          intent(inout) :: array(:, :)
        character(len=*), intent(in)    :: var_name

        integer :: ncid, varid
        integer :: x_dimid, y_dimid

        ! Create the NetCDF file. Override file, if it already exists.
        call check(nf90_create(file_name, NF90_CLOBBER, ncid))

        ! Define the dimensions. NetCDF returns the IDs `x_dimid` and `y_dimid`.
        call check(nf90_def_dim(ncid, 'x', size(array, 1), x_dimid))
        call check(nf90_def_dim(ncid, 'y', size(array, 2), y_dimid))

        ! Define the variable type (NF90_INT: 4-byte integer).
        call check(nf90_def_var(ncid, var_name, NF90_INT, [ x_dimid, y_dimid ], varid))

        ! End define mode.
        call check(nf90_enddef(ncid))

        ! Write the data to the file.
        call check(nf90_put_var(ncid, varid, array))

        ! Close the file.
        call check(nf90_close(ncid))
    end subroutine export_int_array

    subroutine import_int_array(file_name, array, var_name)
        character(len=*), intent(in)    :: file_name
        integer,          intent(inout) :: array(:, :)
        character(len=*), intent(in)    :: var_name
        integer                         :: ncid, varid

        ! Open the NetCDF file read-only.
        call check(nf90_open(file_name, NF90_NOWRITE, ncid))

        ! Get the `varid` of the data variable, based on its name.
        call check(nf90_inq_varid(ncid, var_name, varid))

        ! Read the data.
        call check(nf90_get_var(ncid, varid, array))

        ! Close the file.
        call check(nf90_close(ncid))
    end subroutine import_int_array
end program main