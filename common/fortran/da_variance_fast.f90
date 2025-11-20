program da_variance_fast
  use netcdf
  use, intrinsic :: ieee_arithmetic
  implicit none
!================= Parameters =================
  integer, parameter :: RK = selected_real_kind(12, 100)   ! ~double
  integer, parameter :: DOY = 366
  real(RK), parameter :: FILL_THRESH = 1.0e10_RK
  integer, parameter :: MAXVARS = 4           ! up to 4 anomaly vars
!================= CLI & config ================
  character(len=512) :: ol_list, da_list, outdir, ol_prefix, da_prefix
  character(len=256) :: anom_csv
  real(RK) :: tempK, snow_eps, min_valid_frac
  integer :: clim_window
!================= File lists ==================
  character(len=4096), allocatable :: ol_files(:), da_files(:)
  integer :: nol, nda
!================= Grid & indices ==============
  integer :: ntiles
  real(RK), allocatable :: lat(:), lon(:)
  logical,  allocatable :: keep_tile(:)
  integer,  allocatable :: keep_idx(:)
  integer :: nkeep
!================= Climatology accumulators ====
  real(RK), allocatable :: clim_sum(:,:), clim_cnt(:,:)     ! (DOY, nkeep)
  real(RK), allocatable :: clim(:,:)                        ! (DOY, nkeep)
!================= Variance accumulators =========
  type Welf
     real(RK) :: n=0.0_RK, mean=0.0_RK, M2=0.0_RK
  end type Welf
  type(Welf), allocatable :: w_ol(:), w_da(:)
!================= Var names ===================
  integer :: nvars, iv
  character(len=32) :: varnames(MAXVARS)
!================= Misc ========================
  integer :: i, j, ierr
  logical :: ok
!---------------- Parse args ----------------
  call parse_arg('--ol-list', ol_list, ok);      if(.not.ok) call die('missing --ol-list')
  call parse_arg('--da-list', da_list, ok);      if(.not.ok) call die('missing --da-list')
  call parse_arg('--outdir',  outdir,  ok);      if(.not.ok) outdir = './yearly_outputs'
  call parse_arg('--ol-prefix', ol_prefix, ok);  if(.not.ok) ol_prefix='LS_OLv8_M36'
  call parse_arg('--da-prefix', da_prefix, ok);  if(.not.ok) da_prefix='LS_DAv8_M36'
  call parse_arg('--anom-vars', anom_csv, ok);   if(.not.ok) anom_csv='SFMC,RZMC'
  call parse_real('--temp-K', tempK, ok);        if(.not.ok) tempK=275.15_RK
  call parse_real('--snow-eps', snow_eps, ok);   if(.not.ok) snow_eps=1.0e-2_RK
  call parse_real('--min-valid-frac', min_valid_frac, ok); if(.not.ok) min_valid_frac=0.70_RK
  call parse_int('--clim-window', clim_window, ok); if(.not.ok) clim_window=31

  call split_csv(trim(anom_csv), varnames, nvars)
  if(nvars<1) call die('no anomaly variables parsed')

  call ensure_dir(outdir)

!---------------- Read file lists ----------------
  call read_list(ol_list, ol_files, nol)
  call read_list(da_list, da_files, nda)
  if(nol<1 .or. nda<1) call die('empty file list(s)')
  write(*,'(a,i0,a,i0)') 'INFO: files OL=', nol, ', DA=', nda

!---------------- Discover grid ------------------
  call read_latlon(ol_files(1), lat, lon, ntiles)

  allocate(keep_tile(ntiles))

!---------------- PASS A.1: keep_tile -----------
  write(*,*) 'INFO: PASS A: building keep_tile (valid-day counts BOTH expts)'
  call build_keep_mask(ol_files, nol, varnames, nvars, tempK, snow_eps, ntiles, lat, lon, &
                       da_files, nda, min_valid_frac, keep_tile)

  nkeep = count(keep_tile)
  allocate(keep_idx(nkeep)); j=0
  do i=1,ntiles
    if(keep_tile(i)) then; j=j+1; keep_idx(j)=i; end if
  end do
  write(*,'(a,i0,a,i0)') 'INFO: keep tiles = ', nkeep, '/', ntiles

  call write_keep_nc(trim(outdir)//'/keep_tile.nc', keep_tile, lat, lon, ntiles)

!---------------- PASS A.2: OL climatology ------
  write(*,*) 'INFO: PASS A: accumulating OL climatology (DOY) on kept tiles'
  allocate(clim_sum(DOY,nkeep), clim_cnt(DOY,nkeep), clim(DOY,nkeep))
  clim_sum = 0.0_RK; clim_cnt = 0.0_RK

  do iv=1,nvars
    clim_sum = 0.0_RK; clim_cnt = 0.0_RK
    call accumulate_ol_climatology(ol_files, nol, trim(varnames(iv)), tempK, snow_eps, &
                                   ntiles, keep_idx, nkeep, clim_sum, clim_cnt)

    call fill_missing_366(clim_sum, clim_cnt, DOY, nkeep)
    call finalize_mean_and_smooth(clim_sum, clim_cnt, DOY, nkeep, clim_window, clim)

    call write_clim_nc(trim(outdir)//'/OLv8_climatology_DOY_smooth_kept.nc', &
                       trim(varnames(iv)), clim, DOY, keep_idx, lat, lon, nkeep, iv==1)
  end do

!---------------- PASS B: variance of anomalies -
  write(*,*) 'INFO: PASS B: computing variance of anomalies (OL vs DA)'
  do iv=1,nvars
    call variance_pass(ol_files, nol, da_files, nda, trim(varnames(iv)), tempK, snow_eps, &
                       ntiles, keep_idx, nkeep, trim(outdir)//'/OLv8_climatology_DOY_smooth_kept.nc', &
                       w_ol, w_da)

    call write_variance_summary(trim(outdir)//'/variance_impact_summary.nc', trim(varnames(iv)), &
                                w_ol, w_da, keep_idx, lat, lon, nkeep, iv==1)
    deallocate(w_ol, w_da)
  end do

  write(*,*) 'INFO: Done.'
contains
!======================================================================
subroutine die(s)
  character(len=*), intent(in) :: s
  write(*,*) 'FATAL: ', trim(s)
  stop 2
end subroutine

subroutine ensure_dir(path)
  character(len=*), intent(in) :: path
  logical :: ex
  integer :: ios

  inquire(file=trim(path), exist=ex)
  if (.not. ex) then
     call execute_command_line('mkdir -p "'//trim(path)//'"', exitstat=ios)
     if (ios /= 0) call die('mkdir failed: '//trim(path))
  end if
end subroutine

!---------------- Simple CLI parsers -----------------------------------
subroutine parse_arg(flag, val, ok)
  character(len=*), intent(in) :: flag
  character(len=*), intent(out) :: val
  logical, intent(out) :: ok
  character(len=4096) :: cmd
  integer :: i
  ok=.false.; val=''
  call get_command(cmd)
  i = index(cmd, flag)
  if(i>0) then
     val = adjustl(cmd(i+len_trim(flag):))
     if(len_trim(val)>0) then
        if(val(1:1)==' ') val=val(2:)
        i = index(val, ' ')
        if(i>0) val=val(1:i-1)
        ok=.true.
     end if
  end if
end subroutine

subroutine parse_real(flag, val, ok)
  character(len=*), intent(in) :: flag
  real(RK), intent(out) :: val
  logical, intent(out) :: ok
  character(len=256) :: s
  call parse_arg(flag, s, ok)
  if(ok) read(s,*) val
end subroutine

subroutine parse_int(flag, val, ok)
  character(len=*), intent(in) :: flag
  integer, intent(out) :: val
  logical, intent(out) :: ok
  character(len=256) :: s
  call parse_arg(flag, s, ok)
  if(ok) read(s,*) val
end subroutine

subroutine split_csv(csv, arr, n)
  character(len=*), intent(in) :: csv
  character(len=*), intent(out) :: arr(:)
  integer, intent(out) :: n
  integer :: p, q, m
  character(len=512) :: s, tok
  s=trim(csv); m=0; p=1
  do
     q = index(s(p:), ',')
     if(q==0) then
        tok=adjustl(s(p:))
        if(len_trim(tok)>0) then; m=m+1; arr(m)=tok; end if
        exit
     else
        tok=adjustl(s(p:p+q-2))
        if(len_trim(tok)>0) then; m=m+1; arr(m)=tok; end if
        p=p+q
     end if
  end do
  n=m
end subroutine

!---------------- Read text list -------------------------------
subroutine read_list(listfile, files, n)
  character(len=*), intent(in) :: listfile
  character(len=4096), allocatable, intent(out) :: files(:)
  integer, intent(out) :: n
  character(len=4096) :: line
  integer :: u, cnt, ios
  n=0; cnt=0
  open(newunit=u, file=trim(listfile), status='old', action='read', iostat=ios)
  if(ios/=0) call die('cannot open list: '//trim(listfile))
  do
     read(u,'(A)', iostat=ios) line
     if(ios/=0) exit
     if(len_trim(line)>0) cnt=cnt+1
  end do
  rewind(u)
  allocate(files(cnt)); n=0
  do
     read(u,'(A)', iostat=ios) line
     if(ios/=0) exit
     if(len_trim(line)>0) then; n=n+1; files(n)=adjustl(trim(line)); end if
  end do
  close(u)
end subroutine

!---------------- NetCDF helpers -------------------------------
subroutine nc_check(st, where)
  integer, intent(in) :: st
  character(len=*), intent(in) :: where
  if(st/=nf90_noerr) call die(trim(where)//': '//trim(nf90_strerror(st)))
end subroutine

subroutine read_latlon(path, lat, lon, ntiles)
  character(len=*), intent(in) :: path
  real(RK), allocatable, intent(out) :: lat(:), lon(:)
  integer, intent(out) :: ntiles
  integer :: ncid, vid, did, st
  integer :: dimlen
  st = nf90_open(trim(path), NF90_NOWRITE, ncid); call nc_check(st,'open latlon')
  st = nf90_inq_dimid(ncid, 'tile', did); call nc_check(st,'inq_dim tile')
  st = nf90_inquire_dimension(ncid, did, len=dimlen); call nc_check(st,'inq_dim len')
  ntiles = dimlen
  allocate(lat(ntiles), lon(ntiles))
  st = nf90_inq_varid(ncid, 'lat', vid); call nc_check(st,'inq_var lat')
  st = nf90_get_var(ncid, vid, lat);     call nc_check(st,'get lat')
  st = nf90_inq_varid(ncid, 'lon', vid); call nc_check(st,'inq_var lon')
  st = nf90_get_var(ncid, vid, lon);     call nc_check(st,'get lon')
  st = nf90_close(ncid)
end subroutine

! Read 1D (or time=1, tile) into 1D tile

subroutine read_var_1d_tile(path, vname, arr, ntiles)
  use netcdf
  implicit none
  character(len=*), intent(in) :: path, vname
  real(RK),        intent(out) :: arr(:)     ! size ntiles
  integer,         intent(in)  :: ntiles

  integer :: ncid, vid, st, nd, dims(4)
  character(len=NF90_MAX_NAME) :: dname
  integer :: dimlen, d, tile_dim, time_dim
  integer :: start(4), count(4)
  character(len=32) :: c_nd, c_dimlen, c_ntiles

  tile_dim = -1
  time_dim = -1

  st = nf90_open(trim(path), NF90_NOWRITE, ncid);                         call nc_check(st,'open '//trim(vname))
  st = nf90_inq_varid(ncid, trim(vname), vid);                            call nc_check(st,'inq_var '//trim(vname))
  st = nf90_inquire_variable(ncid, vid, ndims=nd, dimids=dims);           call nc_check(st,'inq dims '//trim(vname))

  if (nd < 1 .or. nd > 2) then
     write(c_nd,'(I0)') nd
     call die('read_var_1d_tile: '//trim(vname)//' has unsupported ndims='//trim(c_nd))
  end if

  ! Identify tile and (optional) time dims by name
  do d = 1, nd
     call nc_check(nf90_inquire_dimension(ncid, dims(d), name=dname, len=dimlen), 'inq dim '//trim(vname))
     if (trim(dname) == 'tile') then
        tile_dim = d
        if (dimlen /= ntiles) then
           write(c_dimlen,'(I0)') dimlen
           write(c_ntiles,'(I0)') ntiles
           call die('tile length mismatch in '//trim(vname)//' file='//trim(path)// &
                    ' (file='//trim(c_dimlen)//' vs expected='//trim(c_ntiles)//')')
        end if
     else if (trim(dname) == 'time') then
        time_dim = d
     end if
  end do

  ! Fallback if no dimension literally named "tile": pick the one whose length==ntiles
  if (tile_dim < 0) then
     do d = 1, nd
        call nc_check(nf90_inquire_dimension(ncid, dims(d), name=dname, len=dimlen), 'inq dim2 '//trim(vname))
        if (dimlen == ntiles) then
           tile_dim = d
           exit
        end if
     end do
  end if
  if (tile_dim < 0) then
     call die('cannot locate tile axis for '//trim(vname)//' in '//trim(path))
  end if

  ! Build start/count with first time slice if time is present
  start = 1; count = 1
  if (nd == 1) then
     ! (tile)
     count(1) = ntiles
  else
     ! 2-D: set time=1 (if present), tile=all
     if (time_dim > 0) then
        start(time_dim) = 1
        count(time_dim) = 1
     else
        ! No explicit "time" name; assume the other dim is time-like → first index
        if (tile_dim == 1) then
           start(2) = 1; count(2) = 1
        else
           start(1) = 1; count(1) = 1
        end if
     end if
     count(tile_dim) = ntiles
  end if

  st = nf90_get_var(ncid, vid, arr, start=start(1:nd), count=count(1:nd)); call nc_check(st,'get '//trim(vname))
  call nc_check(nf90_close(ncid), 'close '//trim(vname))
end subroutine


!---------------- Date parsing from filename -------------------
integer function doy_from_fname(path)
  character(len=*), intent(in) :: path
  integer :: y, m, d
  integer :: pos
  character(len=8) :: ymd
  pos = index(path, '_1200z.nc4')
  if(pos==0) then
     doy_from_fname = -1; return
  end if
  ymd = path(pos-8:pos-1)
  read(ymd(1:4),*) y
  read(ymd(5:6),*) m
  read(ymd(7:8),*) d
  doy_from_fname = doy_from_ymd(y,m,d)
end function

integer function doy_from_ymd(y,m,d)
  integer, intent(in) :: y,m,d
  integer, parameter :: cum(12) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
  logical :: leap
  leap = (mod(y,400)==0) .or. ((mod(y,4)==0) .and. (mod(y,100)/=0))
  doy_from_ymd = cum(m) + d
  if(leap .and. m>2) doy_from_ymd = doy_from_ymd + 1
end function

!---------------- Build keep mask ------------------------------
subroutine build_keep_mask(ol_files, nol, varnames, nvars, tempK, snow_eps, ntiles, lat, lon, &
                           da_files, nda, min_frac, keep_tile)
  character(len=*), intent(in) :: ol_files(:), da_files(:), varnames(:)
  integer, intent(in) :: nol, nda, nvars, ntiles
  real(RK), intent(in) :: tempK, snow_eps, min_frac
  real(RK), intent(in) :: lat(:), lon(:)
  logical, intent(out) :: keep_tile(:)
  integer :: i, v
  real(RK), allocatable :: tmp(:), tsoil(:), frsn(:)
  integer, allocatable :: cnt_ol(:), cnt_da(:)
  logical, allocatable :: ok(:)

  allocate(cnt_ol(ntiles)); cnt_ol=0
  allocate(cnt_da(ntiles)); cnt_da=0
  allocate(tmp(ntiles), tsoil(ntiles), frsn(ntiles))
  allocate(ok(ntiles))

  ! OL counts
  do i=1,nol
     call read_var_1d_tile(ol_files(i), 'TSOIL1', tsoil, ntiles)
     call read_var_1d_tile(ol_files(i), 'FRLANDSNO', frsn, ntiles)
     ok = .true.
     do v=1,nvars
        call read_var_1d_tile(ol_files(i), trim(varnames(v)), tmp, ntiles)
        ok = ok .and. .not.(abs(tmp)>FILL_THRESH) .and. (tsoil>=tempK) .and. (frsn<=snow_eps)
     end do
     cnt_ol = cnt_ol + merge(1,0,ok)
  end do

  ! DA counts
  do i=1,nda
     call read_var_1d_tile(da_files(i), 'TSOIL1', tsoil, ntiles)
     call read_var_1d_tile(da_files(i), 'FRLANDSNO', frsn, ntiles)
     ok = .true.
     do v=1,nvars
        call read_var_1d_tile(da_files(i), trim(varnames(v)), tmp, ntiles)
        ok = ok .and. .not.(abs(tmp)>FILL_THRESH) .and. (tsoil>=tempK) .and. (frsn<=snow_eps)
     end do
     cnt_da = cnt_da + merge(1,0,ok)
  end do

  keep_tile = ( real(cnt_ol,RK) >= min_frac*real(nol,RK) ) .and. &
              ( real(cnt_da,RK) >= min_frac*real(nda,RK) )

  deallocate(tmp, tsoil, frsn, cnt_ol, cnt_da, ok)
end subroutine

!---------------- Accumulate OL climatology --------------------
subroutine accumulate_ol_climatology(files, n, vname, tempK, snow_eps, ntiles, keep_idx, nkeep, sumd, cntd)
  character(len=*), intent(in) :: files(:), vname
  integer, intent(in) :: n, ntiles, keep_idx(:), nkeep
  real(RK), intent(in) :: tempK, snow_eps
  real(RK), intent(inout) :: sumd(:,:), cntd(:,:)      ! (DOY, nkeep)
  integer :: i, doy, kk, t
  real(RK), allocatable :: a(:), tsoil(:), frsn(:)
  real(RK) :: val

  allocate(a(ntiles), tsoil(ntiles), frsn(ntiles))

  do i=1,n
     doy = doy_from_fname(files(i)); if(doy<1) cycle
     if(doy==366) doy = 365
     call read_var_1d_tile(files(i), trim(vname), a,     ntiles)
     call read_var_1d_tile(files(i), 'TSOIL1',   tsoil, ntiles)
     call read_var_1d_tile(files(i), 'FRLANDSNO',frsn,  ntiles)

     do kk=1,nkeep
        t = keep_idx(kk); val = a(t)
        if( (abs(val)<=FILL_THRESH) .and. (tsoil(t)>=tempK) .and. (frsn(t)<=snow_eps) ) then
           sumd(doy,kk) = sumd(doy,kk) + val
           cntd(doy,kk) = cntd(doy,kk) + 1.0_RK
        end if
     end do
  end do

  deallocate(a, tsoil, frsn)
end subroutine

subroutine fill_missing_366(sumd, cntd, ndoy, nkeep)
  real(RK), intent(inout) :: sumd(:,:), cntd(:,:)
  integer, intent(in) :: ndoy, nkeep
  integer :: kk
  do kk=1,nkeep
     if(cntd(ndoy,kk)==0.0_RK) then
        sumd(ndoy,kk) = sumd(ndoy-1,kk)
        cntd(ndoy,kk) = cntd(ndoy-1,kk)
     end if
  end do
end subroutine

subroutine finalize_mean_and_smooth(sumd, cntd, ndoy, nkeep, win, outm)
  real(RK), intent(in) :: sumd(:,:), cntd(:,:)
  integer, intent(in) :: ndoy, nkeep, win
  real(RK), intent(out) :: outm(:,:)
  real(RK), allocatable :: mean0(:,:), pad(:,:), kernel(:)
  integer :: i, kk, left, right, L, j

  allocate(mean0(ndoy,nkeep)); mean0 = 0.0_RK
  do kk=1,nkeep
     do i=1,ndoy
        if(cntd(i,kk)>0.0_RK) mean0(i,kk) = sumd(i,kk)/cntd(i,kk)
     end do
  end do

  if(win<=1) then
     outm = mean0
     deallocate(mean0)
     return
  end if

  left = win/2; right = (win-1)/2; L = left+right+1
  allocate(pad(ndoy+left+right, nkeep)); pad=0.0_RK
  allocate(kernel(L)); kernel = 1.0_RK/real(L,RK)
  outm = 0.0_RK

  ! wrap pad
  pad(1:left,:) = mean0(ndoy-left+1:ndoy,:)
  pad(left+1:left+ndoy,:) = mean0(:,:)
  pad(left+ndoy+1:left+ndoy+right,:) = mean0(1:right,:)

  do kk=1,nkeep
     do j=1,ndoy
        outm(j,kk) = sum( pad(j:j+L-1,kk) * kernel )
     end do
  end do

  deallocate(mean0, pad, kernel)
end subroutine

!---------------- Welford helpers ------------------------------
subroutine welf_init(arr, nkeep)
  type(Welf), allocatable, intent(out) :: arr(:)
  integer, intent(in) :: nkeep
  allocate(arr(nkeep))
end subroutine

subroutine welf_push(a, arr, nkeep)
  real(RK), intent(in) :: a(:)         ! anomalies for kept tiles
  type(Welf), intent(inout) :: arr(:)
  integer, intent(in) :: nkeep
  integer :: k
  real(RK) :: x, nloc, delta
  do k=1,nkeep
     if(.not. ieee_is_nan(a(k))) then
        x = a(k)
        nloc = arr(k)%n + 1.0_RK
        delta = x - arr(k)%mean
        arr(k)%mean = arr(k)%mean + delta/nloc
        arr(k)%M2   = arr(k)%M2   + delta*(x - arr(k)%mean)
        arr(k)%n    = nloc
     end if
  end do
end subroutine

pure real(RK) function qnan()
  qnan = ieee_value(0.0_RK, ieee_quiet_nan)
end function

!---------------- Variance pass -------------------------------
subroutine variance_pass(ol_files, nol, da_files, nda, vname, tempK, snow_eps, ntiles, keep_idx, nkeep, clim_nc, w_ol, w_da)
  character(len=*), intent(in) :: ol_files(:), da_files(:), vname, clim_nc
  integer, intent(in) :: nol, nda, ntiles, keep_idx(:), nkeep
  real(RK), intent(in) :: tempK, snow_eps
  type(Welf), allocatable, intent(out) :: w_ol(:), w_da(:)

  real(RK), allocatable :: clim(:,:), a(:), b(:), tsoil(:), frsn(:)
  real(RK), allocatable :: an_ol(:), an_da(:)
  integer :: i, doy, kk, t
  real(RK) :: val

  call welf_init(w_ol, nkeep)
  call welf_init(w_da, nkeep)

  allocate(clim(DOY,nkeep))
  call read_clim_nc(clim_nc, vname, clim, DOY, nkeep)

  allocate(a(ntiles), b(ntiles), tsoil(ntiles), frsn(ntiles))
  allocate(an_ol(nkeep), an_da(nkeep))

  ! OL
  do i=1,nol
     doy = doy_from_fname(ol_files(i)); if(doy<1) cycle
     if(doy==366) doy=365
     call read_var_1d_tile(ol_files(i), trim(vname), a, ntiles)
     call read_var_1d_tile(ol_files(i), 'TSOIL1', tsoil, ntiles)
     call read_var_1d_tile(ol_files(i), 'FRLANDSNO', frsn, ntiles)

     do kk=1,nkeep
        t = keep_idx(kk); val = a(t)
        if( (abs(val)<=FILL_THRESH) .and. (tsoil(t)>=tempK) .and. (frsn(t)<=snow_eps) ) then
           an_ol(kk) = val - clim(doy,kk)
        else
           an_ol(kk) = qnan()
        end if
     end do
     call welf_push(an_ol, w_ol, nkeep)
  end do

  ! DA
  do i=1,nda
     doy = doy_from_fname(da_files(i)); if(doy<1) cycle
     if(doy==366) doy=365
     call read_var_1d_tile(da_files(i), trim(vname), b, ntiles)
     call read_var_1d_tile(da_files(i), 'TSOIL1', tsoil, ntiles)
     call read_var_1d_tile(da_files(i), 'FRLANDSNO', frsn, ntiles)

     do kk=1,nkeep
        t = keep_idx(kk); val = b(t)
        if( (abs(val)<=FILL_THRESH) .and. (tsoil(t)>=tempK) .and. (frsn(t)<=snow_eps) ) then
           an_da(kk) = val - clim(doy,kk)
        else
           an_da(kk) = qnan()
        end if
     end do
     call welf_push(an_da, w_da, nkeep)
  end do

  deallocate(clim, a, b, tsoil, frsn, an_ol, an_da)
end subroutine

!---------------- NetCDF I/O for keep/clim/summary ------------
subroutine write_keep_nc(path, keep_tile, lat, lon, ntiles)
  character(len=*), intent(in) :: path
  logical, intent(in) :: keep_tile(:)
  real(RK), intent(in) :: lat(:), lon(:)
  integer, intent(in) :: ntiles
  integer :: ncid, did, vid, st
  integer, allocatable :: keep_i(:)

  allocate(keep_i(ntiles)); keep_i = merge(1,0,keep_tile)

  st = nf90_create(trim(path), NF90_CLOBBER, ncid); call nc_check(st,'create keep')
  st = nf90_def_dim(ncid, 'tile', ntiles, did); call nc_check(st,'def dim tile')
  st = nf90_def_var(ncid, 'keep_tile', NF90_INT, (/did/), vid); call nc_check(st,'def keep')
  st = nf90_put_att(ncid, vid, 'long_name', 'tile keep mask (1=true,0=false)')
  st = nf90_def_var(ncid, 'lat', NF90_DOUBLE, (/did/), vid); call nc_check(st,'def lat')
  st = nf90_def_var(ncid, 'lon', NF90_DOUBLE, (/did/), vid); call nc_check(st,'def lon')
  st = nf90_enddef(ncid); call nc_check(st,'enddef')
  st = nf90_inq_varid(ncid, 'keep_tile', vid); call nc_check(st,'inq keep')
  st = nf90_put_var(ncid, vid, keep_i); call nc_check(st,'put keep')
  st = nf90_inq_varid(ncid, 'lat', vid); call nc_check(st,'inq lat')
  st = nf90_put_var(ncid, vid, lat);     call nc_check(st,'put lat')
  st = nf90_inq_varid(ncid, 'lon', vid); call nc_check(st,'inq lon')
  st = nf90_put_var(ncid, vid, lon);     call nc_check(st,'put lon')
  st = nf90_close(ncid)
  deallocate(keep_i)
end subroutine

subroutine write_clim_nc(path, vname, clim, ndoy, keep_idx, lat, lon, nkeep, first)
  character(len=*), intent(in) :: path, vname
  real(RK), intent(in) :: clim(:,:)
  integer, intent(in) :: ndoy, keep_idx(:), nkeep
  real(RK), intent(in) :: lat(:), lon(:)
  logical, intent(in) :: first
  integer :: ncid, dim_doy, dim_tile, vid, st

  if(first) then
     ! ---------- DEFINE MODE ----------
     st = nf90_create(trim(path), NF90_CLOBBER, ncid);           call nc_check(st,'create clim')
     st = nf90_def_dim(ncid, 'dayofyear', ndoy, dim_doy);        call nc_check(st,'def doy')
     st = nf90_def_dim(ncid, 'tile',     nkeep, dim_tile);       call nc_check(st,'def tile')

     st = nf90_def_var(ncid, 'tile', NF90_INT,    (/dim_tile/), vid); call nc_check(st,'def tile idx')
     st = nf90_def_var(ncid, 'lat',  NF90_DOUBLE, (/dim_tile/), vid); call nc_check(st,'def lat')
     st = nf90_def_var(ncid, 'lon',  NF90_DOUBLE, (/dim_tile/), vid); call nc_check(st,'def lon')

     ! >>> INSERTED HERE (before first enddef): define the clim var while in define mode <<<
     st = nf90_def_var(ncid, trim(vname)//'_clim', NF90_DOUBLE, (/dim_doy, dim_tile/), vid)
     call nc_check(st,'def var clim')

     ! Leave define mode
     st = nf90_enddef(ncid);                                      call nc_check(st,'enddef')

     ! ---------- DATA MODE ----------
     st = nf90_inq_varid(ncid, 'tile', vid);                      call nc_check(st,'inq tile idx')
     call nc_check(nf90_put_var(ncid, vid, keep_idx), 'put tile idx')
     st = nf90_inq_varid(ncid, 'lat',  vid);                      call nc_check(st,'inq lat')
     call nc_check(nf90_put_var(ncid, vid, lat(keep_idx)), 'put lat kept')
     st = nf90_inq_varid(ncid, 'lon',  vid);                      call nc_check(st,'inq lon')
     call nc_check(nf90_put_var(ncid, vid, lon(keep_idx)), 'put lon kept')

  else
     ! Appending another variable: open, go to define mode, define, then enddef
     st = nf90_open(trim(path), NF90_WRITE, ncid);                call nc_check(st,'open clim append')

     ! You can inq dims in data mode, but we’ll do it consistently after redef too
     st = nf90_redef(ncid);                                       call nc_check(st,'redef clim append')

     st = nf90_inq_dimid(ncid, 'dayofyear', dim_doy);             call nc_check(st,'inq doy')
     st = nf90_inq_dimid(ncid, 'tile',     dim_tile);             call nc_check(st,'inq tile')

     st = nf90_def_var(ncid, trim(vname)//'_clim', NF90_DOUBLE, (/dim_doy, dim_tile/), vid)
     call nc_check(st,'def var clim (append)')

     st = nf90_enddef(ncid);                                      call nc_check(st,'enddef append')
  end if

  ! Now write the climatology field (data mode)
  st = nf90_inq_varid(ncid, trim(vname)//'_clim', vid);           call nc_check(st,'inq put clim')
  st = nf90_put_var(ncid, vid, clim);                             call nc_check(st,'put clim')

  st = nf90_close(ncid);                                          call nc_check(st,'close clim')
end subroutine

subroutine read_clim_nc(path, vname, clim, ndoy, nkeep)
  character(len=*), intent(in) :: path, vname
  real(RK), intent(out) :: clim(:,:)
  integer, intent(in) :: ndoy, nkeep
  integer :: ncid, vid, st
  st = nf90_open(trim(path), NF90_NOWRITE, ncid); call nc_check(st,'open clim')
  st = nf90_inq_varid(ncid, trim(vname)//'_clim', vid); call nc_check(st,'inq clim var')
  st = nf90_get_var(ncid, vid, clim); call nc_check(st,'get clim')
  st = nf90_close(ncid)
end subroutine

subroutine write_variance_summary(path, vname, w_ol, w_da, keep_idx, lat, lon, nkeep, first)
  character(len=*), intent(in) :: path, vname
  type(Welf), intent(in) :: w_ol(:), w_da(:)
  integer, intent(in) :: keep_idx(:), nkeep
  real(RK), intent(in) :: lat(:), lon(:)
  logical, intent(in) :: first
  integer :: ncid, dim_tile, st, vid
  real(RK), allocatable :: var_ol(:), var_da(:), dvar(:), pct(:), ratio(:), n(:)
  integer :: k

  allocate(var_ol(nkeep), var_da(nkeep), dvar(nkeep), pct(nkeep), ratio(nkeep), n(nkeep))
  call finalize_var(w_ol, var_ol, n, nkeep)
  call finalize_var(w_da, var_da, n, nkeep)   ! n reused as placeholder

  dvar  = var_da - var_ol
  do k=1,nkeep
     if(var_ol(k)>0.0_RK .and. ieee_is_finite(var_ol(k))) then
        pct(k)   = 100.0_RK*dvar(k)/var_ol(k)
        ratio(k) = var_da(k)/var_ol(k)
     else
        pct(k)   = qnan()
        ratio(k) = qnan()
     end if
  end do

  if(first) then
     call nc_check(nf90_create(trim(path), NF90_CLOBBER, ncid), 'create varsum')
     call nc_check(nf90_def_dim(ncid, 'tile', nkeep, dim_tile), 'def dim tile')
     call nc_check(nf90_def_var(ncid, 'tile', NF90_INT, (/dim_tile/), vid), 'def tile idx')
     call nc_check(nf90_def_var(ncid, 'lat',  NF90_DOUBLE, (/dim_tile/), vid), 'def lat')
     call nc_check(nf90_def_var(ncid, 'lon',  NF90_DOUBLE, (/dim_tile/), vid), 'def lon')
     call nc_check(nf90_enddef(ncid), 'enddef init')
     call nc_check(nf90_inq_varid(ncid,'tile',vid),'inq tile')
     call nc_check(nf90_put_var(ncid,vid,keep_idx),'put tile')
     call nc_check(nf90_inq_varid(ncid,'lat',vid),'inq lat')
     call nc_check(nf90_put_var(ncid,vid,lat(keep_idx)),'put lat')
     call nc_check(nf90_inq_varid(ncid,'lon',vid),'inq lon')
     call nc_check(nf90_put_var(ncid,vid,lon(keep_idx)),'put lon')
  else
     call nc_check(nf90_open(trim(path), NF90_WRITE, ncid), 'open varsum')
  end if

  call def_and_put_1d(ncid, trim(vname)//'_var_ol', NF90_DOUBLE, var_ol, nkeep)
  call def_and_put_1d(ncid, trim(vname)//'_var_da', NF90_DOUBLE, var_da, nkeep)
  call def_and_put_1d(ncid, trim(vname)//'_dvar',   NF90_DOUBLE, dvar,   nkeep)
  call def_and_put_1d(ncid, trim(vname)//'_pct',    NF90_DOUBLE, pct,    nkeep)
  call def_and_put_1d(ncid, trim(vname)//'_ratio',  NF90_DOUBLE, ratio,  nkeep)
  call def_and_put_1d(ncid, trim(vname)//'_n',      NF90_DOUBLE, n,      nkeep)

  call nc_check(nf90_close(ncid), 'close varsum')
  deallocate(var_ol, var_da, dvar, pct, ratio, n)
end subroutine

subroutine finalize_var(w, var_out, n_out, nkeep)
  type(Welf), intent(in) :: w(:)
  real(RK), intent(out) :: var_out(:), n_out(:)
  integer, intent(in) :: nkeep
  integer :: k
  do k=1,nkeep
     if(w(k)%n>1.0_RK) then
        var_out(k) = w(k)%M2 / (w(k)%n - 1.0_RK)
     else
        var_out(k) = qnan()
     end if
     n_out(k) = w(k)%n
  end do
end subroutine

subroutine def_and_put_1d(ncid, name, xtype, vec, nkeep)
  integer, intent(in) :: ncid, xtype, nkeep
  character(len=*), intent(in) :: name
  real(RK), intent(in) :: vec(:)
  integer :: vid, st, dim_tile
  call nc_check(nf90_inq_dimid(ncid,'tile',dim_tile),'inq dim tile')
  call nc_check(nf90_redef(ncid), 'redef '//name)
  call nc_check(nf90_def_var(ncid, trim(name), xtype, (/dim_tile/), vid), 'def var '//name)
  call nc_check(nf90_enddef(ncid), 'enddef '//name)
  call nc_check(nf90_put_var(ncid, vid, vec), 'put '//name)
end subroutine

end program da_variance_fast

