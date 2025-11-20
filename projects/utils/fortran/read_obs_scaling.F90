program read_obs_scaling
    
    implicit none       
    
    ! ----------------------------------------------------------
    
    ! local variables
    
    real,    parameter :: no_data_stats = -9999.
    
    real,    parameter :: tol = 1e-2
            
    ! -------------------
        
    logical        :: hpol, scale_mean_only
    
    character(300) :: fname
    
    character( 80) :: tmpstring80, obs_param_descr
    character(  2) :: tmpstring2, orbit_flag
    
    integer        :: i, ind, istat, ind_angle, obs_param_orbit, obs_param_pol, obs_param_N_ang
    
    integer        :: asc_flag, N_data_min, N_sclprm, N_ang
    
    real           :: tmpreal, obs_param_ang(10)
    
    integer, dimension(:), allocatable :: sclprm_tile_id
    
    real,    dimension(:), allocatable :: sclprm_ang
    real,    dimension(:), allocatable :: sclprm_lon,      sclprm_lat 
    real,    dimension(:), allocatable :: sclprm_mean_obs, sclprm_std_obs
    real,    dimension(:), allocatable :: sclprm_mean_mod, sclprm_std_mod
    
    character(len=*), parameter :: Iam = 'scale_obs_Tb_zscore'
    character(len=400) :: err_msg

    ! ------------------------------------------------------------------

       
    scale_mean_only = .true.
       
    ! obs_param_nml(31)%descr          = 'SMAP_L1C_Tbh_A'
    ! obs_param_nml(31)%orbit          = 1
    ! obs_param_nml(31)%pol            = 1
    ! obs_param_nml(31)%N_ang          = 1
    ! obs_param_nml(31)%ang(1)         = 40.

    obs_param_descr = 'SMAP_L1C_Tbh_A'
    obs_param_orbit = 1
    obs_param_pol   = 1
    obs_param_N_ang = 1
    obs_param_ang(1) = 40.


    ! assemble the name of the file with scaling parameters
    !
    ! - first 5 chars of this_obs_param%scalename are NOT part of the file name
    ! - different scaling files for each orbit and pentad
    
    
    fname = '../test_data/OLv7_M36_MULTI_L1C_zscore_stats_A_p40.bin'   
    fname = '../test_data/L4SM_OLv7M36_L1C_zscore_stats_A_p40.bin'                                             

    
    open(10, file=fname, form='unformatted',convert='big_endian',status='old', &
         access='SEQUENTIAL', iostat=istat)

    if (istat/=0) then
         write(*,*) 'Error opening file ', trim(fname)
         stop
     end if
   

    ! read file header 
    !
    ! file format of scaling files mirrors that of pre-processed SMOS obs files

    read(10) asc_flag, N_data_min
    read(10)                 ! start time of interval for stats computation (not used)
    read(10)                 ! end   time of interval for stats computation (not used)
    read(10) N_sclprm, N_ang

    ! minimal consistency checks

    if ( (obs_param_orbit==1 .and. asc_flag/=1) .or.       &
         (obs_param_orbit==2 .and. asc_flag/=0)      ) then
         err_msg = 'orbit flag does not match'
         write(*, *) err_msg
    end if
    
    write (*, *) '  asc_flag   = ', asc_flag 
    write (*, *) '  N_data_min = ', N_data_min
    write (*, *) '  N_sclprm   = ', N_sclprm
    write (*, *) '  N_ang      = ', N_ang
    
    allocate(sclprm_ang(     N_ang   ))
    
    allocate(sclprm_lon(     N_sclprm))
    allocate(sclprm_lat(     N_sclprm))
    allocate(sclprm_tile_id( N_sclprm))
    
    allocate(sclprm_mean_obs(N_sclprm))     
    allocate(sclprm_std_obs( N_sclprm))      
    allocate(sclprm_mean_mod(N_sclprm))     
    allocate(sclprm_std_mod( N_sclprm)) 

    ! read angle and location information
    
    read(10) sclprm_ang
    
    read(10) sclprm_lon      !only valid values where obs were available
    read(10) sclprm_lat      !only valid values where obs were available
    read(10) sclprm_tile_id

    ! Write the max and min values lof the lat and lon values to screen
    write(*,*) 'min lon = ', minval(sclprm_lon)
    write(*,*) 'max lon = ', maxval(sclprm_lon)
    write(*,*) 'min lat = ', minval(sclprm_lat)
    write(*,*) 'max lat = ', maxval(sclprm_lat)
    
    ! find the index for the angle of interest
    ! NOTE: after processing of namelist inputs, each species 
    !       has a unique angle (see subroutine read_ens_upd_inputs())
    
    ind_angle = -9999
    
    do i=1,N_ang
       
       if (abs(sclprm_ang(i)-obs_param_ang(1))<0.01)  ind_angle = i
       
    end do
       
    ! need h-pol or v-pol?
    
    if     (obs_param_pol==1) then
       
       hpol = .true.
       
    elseif (obs_param_pol==2) then
       
       hpol = .false.

    end if
    
    ! read scaling parameters
    !
    ! for each field, loop over all angles, read stats for angle of interest
    !
    ! blocks (1- 5) after header: Tbh stats
    ! blocks (6-10) after header: Tbv stats
    
    if (.not. hpol) then  ! in case of V-pol, skip through H-pol entries
       
       do i=1,N_ang ! block 1 - mean_obs Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 2 - std_obs Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 3 - mean_mod Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 4 - std_mod Tbh           
          read(10) 
       end do
       
       do i=1,N_ang ! block 5 - N_data Tbh          
          read(10)  
       end do
       
    end if
    
    ! from each block, read stats for angle of interest
    
    do i=1,N_ang           ! block 1 (h-pol) or 6 (v-pol) - mean_obs 
       
       if (i==ind_angle) then
          read(10) sclprm_mean_obs
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 2 (h-pol) or 7 (v-pol) - std_obs
       
       if (i==ind_angle) then
          read(10) sclprm_std_obs
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 3 (h-pol) or 8 (v-pol) - mean_mod
       
       if (i==ind_angle) then
          read(10) sclprm_mean_mod
       else
          read(10)
       end if
       
    end do
    
    do i=1,N_ang           ! block 4 (h-pol) or 9 (v-pol) - std_mod
       
       if (i==ind_angle) then
          read(10) sclprm_std_mod
       else
          read(10)
       end if
       
    end do
    
    !do i=1,N_ang           ! block 5 (h-pol) or 10 (v-pol) - N_data
    !   
    !   if (i==ind_angle) then
    !      read(10) sclprm_Ndata
    !   else
    !      read(10)
    !   end if
    !   
    !end do
    
    close(10,status='keep')

    ! ----------------------------------------------------------
    ! Print the max and min values of the mean and std values to the screen
    write(*,*) 'min mean obs = ', minval(sclprm_mean_obs)
    write(*,*) 'max mean obs = ', maxval(sclprm_mean_obs)
    write(*,*) 'min std obs = ', minval(sclprm_std_obs)
    write(*,*) 'max std obs = ', maxval(sclprm_std_obs)
    write(*,*) 'min mean mod = ', minval(sclprm_mean_mod)
    write(*,*) 'max mean mod = ', maxval(sclprm_mean_mod)
    write(*,*) 'min std mod = ', minval(sclprm_std_mod)
    write(*,*) 'max std mod = ', maxval(sclprm_std_mod)

    ! Write the data to a file that can be read by a python script
    ! open(20, file='./OLv7_M36_MULTI_L1C_zscore_stats_A_p40.txt', status='unknown')
    open(20, file='./L4SM_OLv7M36_L1C_zscore_stats_A_p40.txt', status='unknown')

    write(20, *) 'sclprm_lon'
    write(20, *) (sclprm_lon(i), i=1,N_sclprm)
    write(20, *) 'sclprm_lat'
    write(20, *) (sclprm_lat(i), i=1,N_sclprm)
    write(20, *) 'sclprm_tile_id'
    write(20, *) (sclprm_tile_id(i), i=1,N_sclprm)

    write(20, *) 'mean_obs'
    write(20, *) (sclprm_mean_obs(i), i=1,N_sclprm)
    write(20, *) 'std_obs'
    write(20, *) (sclprm_std_obs(i), i=1,N_sclprm)
    write(20, *) 'mean_mod'
    write(20, *) (sclprm_mean_mod(i), i=1,N_sclprm)
    write(20, *) 'std_mod'
    write(20, *) (sclprm_std_mod(i), i=1,N_sclprm)

    close(20)

end program read_obs_scaling