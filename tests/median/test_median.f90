!>Test driver for day_median (DigiBog_Hydro_procs.f90, Section 8.0).
!>Reads the case lengths from lengths.txt and the test values from
!>values.bin, both written by calcs/median_test_vectors.R, and writes one
!>median per case to medians.bin. The R script calcs/check_median.R
!>compares these with R's median().
!>
!>Values are read and written as raw doubles (stream access), not text,
!>so the comparison with R is exact. Run from tests/median/.
program test_median

  use hydro_procedures, only: q, day_median

  implicit none

  real(kind=q), dimension(:), allocatable :: values
  integer :: n_cases, n, c

  open(unit=10, file="lengths.txt", status="old", action="read")
  open(unit=11, file="values.bin", status="old", action="read", &
       access="stream", form="unformatted")
  open(unit=20, file="medians.bin", status="replace", action="write", &
       access="stream", form="unformatted")

  read (10, *) n_cases
  do c = 1, n_cases
    read (10, *) n
    allocate(values(n))
    read (11) values
    write (20) day_median(values)
    deallocate(values)
  end do

  close(10)
  close(11)
  close(20)
  write (*, '(A, I0, A)') "Wrote ", n_cases, " medians to medians.bin"

end program test_median
