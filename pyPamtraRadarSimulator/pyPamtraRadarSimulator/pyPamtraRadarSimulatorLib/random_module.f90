module random_module

  use kinds
  implicit none
  integer :: counter = 0

  contains
  subroutine get_random(errorstatus, n, seedval, x_noise)
    ! Description:
    ! returns n random numbers
    ! if seedval is not zero, always the same numbers are returned
    !
    !initiation of random nummer generator taken from
    !taken from http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.1   1/12/2012    creation of file 


    use report_module
    use kinds
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: seedval
    real(kind=dbl), intent(out), dimension(n) :: x_noise
    integer :: i, m, clock
    integer, dimension(:), allocatable :: seed

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=10) :: nameOfRoutine = 'get_random'      
    
      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
      !get the required size of the seed
      call random_seed(size = m)

      allocate(seed(m))
      if (seedval > 0) then
        !if we want always the same random numbers, but still a difference for each profile
        clock=seedval -1 + counter! minus one is for historical reasons
        counter = counter + 1
      else if (seedval < 0) then
        !if we want always the same random numbers, for all profiles the same!!
        clock=seedval -1 ! minus one is for historical reasons
      else
        !get 'real' random numbers
        call system_clock(count=clock) 
        clock=clock +  counter 
        counter = counter + 1
      end if
    
      seed = clock + 37 * (/ (i - 1, i = 1, m) /)
      !seed the seed
      call random_seed(put = seed)
      deallocate(seed)
      !get the random numbers
      call RANDOM_NUMBER(x_noise)

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
  end subroutine get_random

end module random_module
