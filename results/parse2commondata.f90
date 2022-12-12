PROGRAM parse2commondata

   IMPLICIT NONE
   
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   CHARACTER(LEN=120) :: filename,outname
   REAL(KIND=dp),ALLOCATABLE :: dat(:,:)
   INTEGER(KIND=4) :: nlines,narg,i,j,nsys_err
   INTEGER(KIND=4),DIMENSION(5) :: k
   REAL(KIND=dp),PARAMETER :: M=0.938
   CHARACTER(LEN=250) :: line

   k = (/ 3,6,9,10,12 /)
   narg = iargc()
   IF (narg .lt. 2) THEN 
      WRITE(0,"(A)") "Please give two arguments. Input filename and output filename."
      CALL exit(1)
   END IF
   CALL getarg(1,filename)
   CALL getarg(2,outname)
   CALL get_n_lines(filename,nlines)

   ALLOCATE(dat(nlines-2,13))

   OPEN(10,FILE=filename)
   READ(10,"(A)",END=20) line

   READ(10,"(A)",END=20) line
   DO i=1,UBOUND(dat,1)
      READ(10,"(A)",END=20) line
      DO j=LBOUND(dat,2),UBOUND(dat,2)
         CALL read_variable_from_line(line,dat(i,j))
! Calculate y and replace Enu avg by it
         dat(i,9) = dat(i,6) / (2*dat(i,3)*M*dat(i,9))
      END DO
   END DO
   20 CLOSE(10)
   OPEN(60,FILE=outname)
! Need to think about how to handle systematic errors
   nsys_err = 0
   WRITE(60,"(A,I3,I3)") "FASERnu2",nsys_err,UBOUND(dat,1)
   DO j=1,UBOUND(dat,1)
      WRITE(60,"(I3,A4,5ES16.8)") j,"DIS",(dat(j,k(i)),i=1,UBOUND(k,1))
   END DO
   CLOSE(60)
   DEALLOCATE(dat)

END PROGRAM parse2commondata

SUBROUTINE read_variable_from_line(line,var)

   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   CHARACTER(LEN=250),INTENT(INOUT) :: line
   CHARACTER(LEN=250) :: substr
   INTEGER(KIND=4) :: ind
   REAL(KIND=dp),INTENT(OUT) :: var

   line = ADJUSTL(line)
   ind = INDEX(line," ")
   IF (ind .ne. 0) THEN
      substr = line(1:ind)
      READ(substr,*) var
      line = line(ind+1:)
   ELSE
      READ(line,*) var
   END IF

END SUBROUTINE read_variable_from_line

SUBROUTINE get_n_lines(filename, nlines)

   IMPLICIT NONE
   
   CHARACTER(len=120),INTENT(IN) :: filename
   INTEGER(KIND=4),INTENT(OUT) :: nlines
   nlines = 0 
   OPEN(40,FILE=filename)
   DO
      READ(40,*,END=50)
      nlines = nlines + 1
   END DO
   50 CLOSE(40)
END SUBROUTINE get_n_lines

