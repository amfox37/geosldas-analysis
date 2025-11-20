module shared_variables
    implicit none

    integer, parameter :: N_obs_ang_max = 7

    type :: date_time_type              
        integer :: year               ! 4-digit year
        integer :: month              ! month in year
        integer :: day                ! day in month
        integer :: hour               ! hour of day
        integer :: min                ! minute of hour
        integer :: sec                ! seconds of minute
        integer :: pentad             ! pentad of year
        integer :: dofyr              ! day of year
    end type date_time_type

    type :: obs_param_type

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! WARNING: When modifying this derived type make sure that the corresponding
     !          MPI STRUCTURE in module CLSM_ENSDRV_MPI is also updated, as are
     !          any subroutines or operators defined herein
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

     character(40)                 :: descr     ! description
     integer                       :: species   ! identifier for type of measurement

     integer                       :: orbit     ! type of (half-)orbit
                                                !                     0 = n/a  [eg., in situ obs]
                                                !                     1 = ascending
                                                !                     2 = descending
                                                !                     3 = ascending or descending
                                                !                     4 = geostationary

     integer                       :: pol       ! polarization
                                                ! 0 = n/a  [eg., multi-pol. retrieval]
                                                ! 1 = horizontal
                                                ! 2 = vertical 
                                                ! 3 = ...
                                                ! [add 3rd/4th Stokes, HH, HV, VH, VV]

     integer                       :: N_ang     ! # satellite viewing angles in species (radiance obs only)

     real, &
          dimension(N_obs_ang_max) :: ang       ! vector of satellite viewing angles
     
     real                          :: freq      ! frequency [Hz]

     real                          :: FOV       ! field-of-view *radius* 
                                                ! if FOV==0. equate obs footprint w/ tile
                                                ! for details see LDASsa_DEFAULT_inputs ensupd.nml
     character(40)                 :: FOV_units ! FOV units ('km' or 'deg') 

     logical                       :: assim     ! assimilate yes/no? (see also "obs_type")
     logical                       :: scale     ! scale yes/no?
     logical                       :: getinnov  ! compute innovs? (.T. if assim==.T.)

     integer                       :: RTM_ID    ! ID of radiative transfer model 

     integer                       :: bias_Npar ! number of bias states tracked per day
     integer                       :: bias_trel ! e-folding time scale of obs bias memory [s]
     integer                       :: bias_tcut ! cutoff time for confident obs bias est [s]

     real                          :: nodata    ! no-data-value

     character(40)                 :: varname   ! equivalent model variable name (Obs_pred)
     character(40)                 :: units     ! units (eg., 'K' or 'm3/m3')

     character(200)                :: path      ! path to measurements file 
     character(80)                 :: name      ! name identifier for measurements 
     character(200)                :: maskpath  ! path to obs mask file
     character(80)                 :: maskname  ! filename for obs mask
     character(200)                :: scalepath ! path to file with scaling parameters
     character(80)                 :: scalename ! filename for scaling parameters
     character(200)                :: flistpath ! path to file with list of obs file names
     character(80)                 :: flistname ! name of file with list of obs file names

     real                          :: errstd    ! default obs error std
                                  
     real                          :: std_normal_max  ! see pert_param_type
     logical                       :: zeromean        ! see pert_param_type
     logical                       :: coarsen_pert    ! see pert_param_type ("%coarsen")
     real                          :: xcorr           ! see pert_param_type
     real                          :: ycorr           ! see pert_param_type
     
     integer                       :: adapt     ! identifier for adaptive filtering
     
  end type obs_param_type

    type(date_time_type)            :: date_time
    integer                         :: dtstep_assim, N_catd
    type(obs_param_type)            :: this_obs_param
    logical                         :: found_obs
    real, dimension(:), allocatable :: CYG_sm, CYG_sm_std, CYG_lon, CYG_lat, CYG_time
end module shared_variables

program main
    use shared_variables
    implicit none

    ! Initialize or allocate variables as needed
    date_time = date_time_type(2018, 08, 01, 03, 00, 00, -9999, -9999)
    dtstep_assim = 3
    N_catd = 100

    this_obs_param%path = '/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/cygnss_test/data'
    this_obs_param%name = 'cyg.ddmi.s20180801-030000-e20180801-210000.l3.grid-soil-moisture-36km.a32.d33'
    this_obs_param%nodata = -9999.0
    this_obs_param%maskpath = '/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/cygnss_test'
    this_obs_param%maskname = 'ucar_cu_cygnss_sm_v1_static_flags'

    call read_obs_sm_CYGNSS()

contains
  subroutine read_obs_sm_CYGNSS()
    use shared_variables
    use netcdf   
    implicit none

    integer, parameter :: max_obs = 20000       ! max number of obs read by subroutine (expecting < 6 hr assim window)

    character(300)       :: tmpfname, tmpname, err_msg, tmpmaskname
    character(  2)       :: MM, DD, HH, MI
    character(  4)       :: YYYY
    integer              :: pos, i, idx, N_obs, j
    integer              :: ierr, ncid
    integer              :: lon_dimid, lat_dimid, time_dimid, timeslices_dimid, startstop_dimid
    integer              :: sm_d_varid, sm_subd_varid, sigma_d_varid, sigma_subd_varid, timeintervals_varid, lat_varid, lon_varid
    integer              :: N_lon, N_lat, N_time, N_timeslices, N_startstop

    ! for mask read
    integer              :: lon_dimid_m, lat_dimid_m
    integer              :: longitudes_m_varid, latitudes_m_varid, small_SM_range_varid, poor_SMAP_varid, high_ubrmsd_varid, few_obs_varid, low_signal_varid
    integer              :: N_lon_m, N_lat_m
    real, allocatable    :: latitudes_m(:,:), longitudes_m(:,:)
    integer, allocatable :: small_SM_range(:,:), poor_SMAP(:,:), high_ubrmsd(:,:), few_obs(:,:), low_signal(:,:)

    logical              :: file_exists
    real, allocatable    :: sm_d(:,:,:), sm_subd(:,:,:), sigma_d(:,:,:), sigma_subd(:,:,:)
    real, allocatable    :: timeintervals(:,:), latitudes(:,:), longitudes(:,:)



    real,      dimension(:),     pointer     :: tmp_lon(:), tmp_lat(:), tmp_obs(:), tmp_err(:), tmp_time(:)
    real                 :: tmp1_lon(max_obs), tmp1_lat(max_obs), tmp1_obs(max_obs), tmp1_err(max_obs), tmp1_time(max_obs)

    nullify( tmp_lon, tmp_lat, tmp_obs, tmp_err, tmp_time )

    write (YYYY,'(i4.4)') date_time%year
    write (MM,  '(i2.2)') date_time%month
    write (DD,  '(i2.2)') date_time%day 
    write (HH,  '(i2.2)') date_time%hour 
    write (MI,  '(i2.2)') date_time%min

    tmpname = trim(this_obs_param%name)

    ! Replace s20180801 with sYYYYMMDD
    pos = index(tmpname, "s20180801")
    if (pos > 0) then
       tmpname = tmpname(1:pos-1) // "s" // YYYY // MM // DD // tmpname(pos+9:)
    endif

    ! Replace e20180801 with eYYYYMMDD
    pos = index(tmpname, "e20180801")
    if (pos > 0) then
       tmpname = tmpname(1:pos-1) // "e" // YYYY // MM // DD // tmpname(pos+9:)
    endif
        
    tmpfname = trim(this_obs_param%path) // '/Y' // YYYY // '/M' // MM // '/' // trim(tmpname) // '.nc'

    write(*,*) 'Reading CYGNSS soil moisture data from file: ', trim(tmpfname)

    ! Check if file exists

    inquire(file=tmpfname, exist=file_exists)

    if (.not. file_exists) then
       err_msg = 'CYGNSS SM obs file not found!'
!           call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! Open the NetCDF observation file
    ierr = nf90_open(trim(tmpfname), nf90_nowrite, ncid)

    ! get variable dimension IDs
    ierr = nf90_inq_dimid(ncid, 'lon',        lon_dimid)
    ierr = nf90_inq_dimid(ncid, 'lat',        lat_dimid)
    ierr = nf90_inq_dimid(ncid, 'time',       time_dimid)    
    ierr = nf90_inq_dimid(ncid, 'timeslices', timeslices_dimid)
    ierr = nf90_inq_dimid(ncid, 'startstop',  startstop_dimid)

    ! dimensions sizes
    ierr = nf90_inquire_dimension(ncid, lon_dimid,        len=N_lon)
    ierr = nf90_inquire_dimension(ncid, lat_dimid,        len=N_lat)
    ierr = nf90_inquire_dimension(ncid, time_dimid,       len=N_time)
    ierr = nf90_inquire_dimension(ncid, timeslices_dimid, len=N_timeslices)
    ierr = nf90_inquire_dimension(ncid, startstop_dimid,  len=N_startstop)

    ! Print dimensions for debugging
    write(*,*) 'Dimensions:'
    write(*,*) 'N_lon = ', N_lon
    write(*,*) 'N_lat = ', N_lat
    write(*,*) 'N_time = ', N_time
    write(*,*) 'N_timeslices = ', N_timeslices
    write(*,*) 'N_startstop = ', N_startstop

    ! get variable IDs
    ierr = nf90_inq_varid(ncid, 'SM_daily',       sm_d_varid)
    ierr = nf90_inq_varid(ncid, 'SM_subdaily',    sm_subd_varid)
    ierr = nf90_inq_varid(ncid, 'SIGMA_daily',    sigma_d_varid)
    ierr = nf90_inq_varid(ncid, 'SIGMA_subdaily', sigma_subd_varid)
    ierr = nf90_inq_varid(ncid, 'timeintervals',  timeintervals_varid)
    ierr = nf90_inq_varid(ncid, 'latitude',       lat_varid)
    ierr = nf90_inq_varid(ncid, 'longitude',      lon_varid)

    ! allocate memory for the variables
    allocate(sm_d(N_lon, N_lat, N_time))
    allocate(sm_subd(N_lon, N_lat, N_timeslices))
    allocate(sigma_d(N_lon, N_lat, N_time))
    allocate(sigma_subd(N_lon, N_lat, N_timeslices))
    allocate(timeintervals(N_startstop, N_timeslices))
    allocate(latitudes(N_lon, N_lat))
    allocate(longitudes(N_lon, N_lat))

    ! read the variables
    ierr = nf90_get_var(ncid, sm_d_varid,       sm_d)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading sm_d:', ierr
    end if
    ierr = nf90_get_var(ncid, sm_subd_varid,    sm_subd)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading sm_subd:', ierr
    end if
    ierr = nf90_get_var(ncid, sigma_d_varid,    sigma_d)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading sigma_d:', ierr
    end if
    ierr = nf90_get_var(ncid, sigma_subd_varid, sigma_subd)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading sigma_subd:', ierr
    end if
    ierr = nf90_get_var(ncid, timeintervals_varid, timeintervals)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading timeintervals:', ierr
    end if
    ierr = nf90_get_var(ncid, lat_varid, latitudes)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading latitudes:', ierr
    end if
    ierr = nf90_get_var(ncid, lon_varid, longitudes)
    if (ierr /= nf90_noerr) then
        write(*,*) 'Error reading longitudes:', ierr
    end if
    
    ! write the first 10 values of sm_d
    write(*,*) 'SM_daily:'
    do i = 1, 10
       write(*,*) sm_d(1,i,1)
    end do

    ! write the first 10 values of sm_subd
    write(*,*) 'SM_subdaily:'
    do i = 1, 10
       write(*,*) sm_subd(1,i,1)
    end do

    ! write the first 10 values of sigma_d
    write(*,*) 'SIGMA_daily:'
    do i = 1, 10
       write(*,*) sigma_d(1,i,1)
    end do

    ! write the first 10 values of sigma_subd
    write(*,*) 'SIGMA_subdaily:'
    do i = 1, 10
       write(*,*) sigma_subd(1,i,1)
    end do

    ! Write out time intervals
    write(*,*) 'Time intervals:' 
    do i = 1, N_timeslices
       write(*,*) timeintervals(1,i), timeintervals(2,i)
    end do

    ! get name for CYGNSS mask file

    tmpmaskname = trim(this_obs_param%maskpath) // '/' // trim(this_obs_param%maskname) // '.nc'

    inquire(file=tmpfname, exist=file_exists)

    if (.not. file_exists) then
       err_msg = 'CYGNSS mask file not found!'
!           call ldas_abort(LDAS_GENERIC_ERROR, Iam, err_msg)
    end if

    ! open the CYGNSS mask file

    ierr = nf90_open(trim(tmpmaskname), nf90_nowrite, ncid)

    ! get variable dimension IDs
    ierr = nf90_inq_dimid(ncid, 'lon', lon_dimid)
    ierr = nf90_inq_dimid(ncid, 'lat', lat_dimid)

    ! dimensions sizes
    ierr = nf90_inquire_dimension(ncid, lon_dimid, len=N_lon_m)
    ierr = nf90_inquire_dimension(ncid, lat_dimid, len=N_lat_m)

    ! get variable IDs
    ierr = nf90_inq_varid(ncid, 'longitude',           longitudes_m_varid)
    ierr = nf90_inq_varid(ncid, 'latitude',            latitudes_m_varid)
    ierr = nf90_inq_varid(ncid, 'flag_small_SM_range', small_SM_range_varid)
    ierr = nf90_inq_varid(ncid, 'flag_poor_SMAP',      poor_SMAP_varid)
    ierr = nf90_inq_varid(ncid, 'flag_high_ubrmsd',    high_ubrmsd_varid)
    ierr = nf90_inq_varid(ncid, 'flag_few_obs',        few_obs_varid)
    ierr = nf90_inq_varid(ncid, 'flag_low_signal',     low_signal_varid)


    ! allocate memory for the variables
    allocate(latitudes_m(N_lon_m, N_lat_m))
    allocate(longitudes_m(N_lon_m, N_lat_m))
    allocate(small_SM_range(N_lon_m, N_lat_m))
    allocate(poor_SMAP(N_lon_m, N_lat_m))
    allocate(high_ubrmsd(N_lon_m, N_lat_m))
    allocate(few_obs(N_lon_m, N_lat_m))
    allocate(low_signal(N_lon_m, N_lat_m))

    ! read the variables
    ierr = nf90_get_var(ncid, latitudes_m_varid,    latitudes_m)
    ierr = nf90_get_var(ncid, longitudes_m_varid,   longitudes_m)
    ierr = nf90_get_var(ncid, small_SM_range_varid, small_SM_range)
    ierr = nf90_get_var(ncid, poor_SMAP_varid,      poor_SMAP)
    ierr = nf90_get_var(ncid, high_ubrmsd_varid,    high_ubrmsd)
    ierr = nf90_get_var(ncid, few_obs_varid,        few_obs)
    ierr = nf90_get_var(ncid, low_signal_varid,     low_signal)

    ! Print the first 10 values of each variable
    write(*,*) 'latitudes_m:', latitudes_m(1:10, 1)
    write(*,*) 'longitudes_m:', longitudes_m(1:10, 1)
    write(*,*) 'small_SM_range:', small_SM_range(1:10, 1)
    write(*,*) 'poor_SMAP:', poor_SMAP(1:10, 1)
    write(*,*) 'high_ubrmsd:', high_ubrmsd(1:10, 1)
    write(*,*) 'few_obs:', few_obs(1:10, 1)
    write(*,*) 'low_signal:', low_signal(1:10, 1) 

    idx = -1

    do i = 1, N_timeslices
        if ((date_time%hour == 3  .and. timeintervals(1,i) == 0.0)  .or. &
            (date_time%hour == 9  .and. timeintervals(1,i) == 0.25) .or. &
            (date_time%hour == 15 .and. timeintervals(1,i) == 0.5)  .or. &
            (date_time%hour == 21 .and. timeintervals(1,i) == 0.75)) then
            idx = i
            exit
        end if
    end do

    if (idx == -1) then
        write(*,*) 'Error: No matching time interval found for HH = ', HH
    else
        write(*,*) 'Index of matching time interval: ', idx
    end if

    ! fill tmp arrays
    N_obs = 0

    do i = 1, N_lon
        do j = 1, N_lat
            if (sm_subd(i,j,idx) .ne. this_obs_param%nodata) then
                N_obs = N_obs + 1
                if (N_obs > max_obs) then
                    write(*,*) 'Error: Number of observations exceeds max_obs'
                    stop
                end if
                tmp1_lon(N_obs) = longitudes(i,j)
                tmp1_lat(N_obs) = latitudes(i,j)
                tmp1_obs(N_obs) = sm_subd(i,j,idx)
                tmp1_err(N_obs) = sigma_subd(i,j,idx)
                tmp1_time(N_obs) = date_time%year
            end if
        end do
    end do

    write(*,*) 'Number of observations: ', N_obs

    allocate(tmp_lon(N_obs))
    allocate(tmp_lat(N_obs))
    allocate(tmp_obs(N_obs))
    allocate(tmp_err(N_obs))
    allocate(tmp_time(N_obs))

    tmp_lon = tmp1_lon(1:N_obs)
    tmp_lat = tmp1_lat(1:N_obs)
    tmp_obs = tmp1_obs(1:N_obs)
    tmp_err = tmp1_err(1:N_obs)
    tmp_time = tmp1_time(1:N_obs)

    ! work out which time interval we want to read based on HH.
    ! 

    ! close the file
    ierr = nf90_close(ncid)

        ! Your subroutine code here
        ! ...

    end subroutine read_obs_sm_CYGNSS        

end program main