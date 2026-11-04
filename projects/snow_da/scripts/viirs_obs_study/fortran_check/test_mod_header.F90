! module header for the extracted reader routines; stubs replace the GEOSldas
! logging and abort routines (build.sh appends the routines and "end module")
module test_mod
  use, intrinsic :: iso_c_binding,   only: c_char
  use, intrinsic :: iso_fortran_env, only: int16
  use io_hdf5,          only: hdf5read
  use hdf4_fortran_api, only: sfstart, sfn2index, sfselect, sfginfo, sfrdata, sfendacc, sfend
  implicit none
  logical, parameter :: logit = .true.
  integer, parameter :: logunit = 6
  integer, parameter :: LDAS_GENERIC_ERROR = 1
contains
  subroutine ldas_abort(code, Iam, msg)
    integer, intent(in) :: code
    character(*), intent(in) :: Iam, msg
    write(*,*) 'ABORT: ', Iam, ': ', trim(msg); stop 1
  end subroutine ldas_abort
