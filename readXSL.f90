
module precisions
implicit none
integer, parameter :: dp=kind(1.0d0)    
end module precisions
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

PROGRAM read_rewrite
USE precisions
IMPLICIT NONE
INTEGER :: I, Nx_max, iu, Nx, Ny, entier,  entierH, nfile, j, k
INTEGER,parameter :: ncolumns = 12, nr_lines_max = 100

  integer, parameter :: max_csv_lines = 100 ! max number of lines in CSV
  integer, parameter :: csv_columns   = 12   ! number of columns in CSV
  integer stat, line_nr, nr, j, computed_number
  character :: separator = ' '
  character(80) :: line
  type(csv_header_line) :: csv_header
  type(csv_data_line), dimension(max_csv_lines) :: csv_data  




!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

OPEN(19,FILE = '509_alongtrack_profiles.xls', status = 'old', action = 'read')

if (stat .ne. 0) then
    write(*,*) 'File cannot be opened !'
    go to 99
end if

OPEN(50 ,FILE = "509_alongtrack_profiles.xls.txt")
    
  write(*,*) 'Reading xsl-file...' 
  ! process file
  line_nr = 0
  do while (.true.) 
            read(19,'(8(E16.8,1x))') line
            line_nr = line_nr + 1
            WRITE(50,'(8(E16.8,1x))') line 
            nr = line_nr - 1
  ENDDO

       close(19)
       close(50) 


  ! close file
  99 continue
  close(1)
  write(*,*) 'Done.' 

  write(*,*)  'Number of all lines found in CSV  = ', line_nr
  write(*,*)  'Number of data lines found in CSV = ', nr

  ! write the data
  write(*, '(A)') '****************************************'
  write(*,'(A10, A10, A10)') csv_header%c01, &
                             adjustr(csv_header%c02), & 
                             '3*NUMCOL'
  write(*, '(A)') '****************************************'
  do j = 1, nr
     computed_number = 3 * csv_data(j)%c02
     write (*,'(A10, I10, I10)') csv_data(j)%c01, &
                                csv_data(j)%c02, &
                                computed_number

  end do
  write(*, '(A)') '****************************************'  

!//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  

END PROGRAM read_rewrite