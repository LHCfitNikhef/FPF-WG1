PROGRAM parse2commondata
   IMPLICIT NONE
   
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   CHARACTER(LEN=220) :: filename1,filename2,outname
   REAL(KIND=dp),ALLOCATABLE :: dat(:,:)
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: fluctuated_central_vals,&
      &central_vals
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: fluctuated_central_vals_os
   INTEGER(KIND=4) :: nlines,narg,i,j,nsys_err
   INTEGER(KIND=4),DIMENSION(7) :: k
   REAL(KIND=dp),PARAMETER :: M=0.938
   CHARACTER(LEN=500) :: line
   REAL(KIND=dp),DIMENSION(:),ALLOCATABLE :: syserr,staterr,err
   REAL(KIND=dp) :: fluctuate
   REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE :: centr
   INTEGER,DIMENSION(2) :: lenfile
   CHARACTER(LEN=40) :: proc,ptype,procstr,pname,experiment

   narg = iargc()
   IF (narg.EQ.1) THEN
      CALL getarg(1,filename1)
      IF (filename1.EQ."--help") THEN
         WRITE(6,"(A)") "This program fluctuates the data for the total &
            & cross-section and writes the result to the commondata format. "
         WRITE(6,"(A)") "p2c PATH-TO-SYSTEMATIC PATH-TO-XSEC CHARGE PROCESS-TYPE"
         WRITE(6,"(A,/,A,/,A,/,A)") "Allowed CHARGES are:","nu","nub","sum"
         WRITE(6,"(A,/,A,/,A)") "Allowed PROCESS-TYPE are:","charm","inclusive"
      END IF
      CALL exit(1)
   END IF
   IF (narg .lt. 6) THEN 
      WRITE(0,"(A,/,A)")  "Please give four arguments",&
  & "Input filenames and output filename followed by process, type and experiment"
      CALL exit(1)
   END IF
   CALL getarg(1,filename1)
   CALL getarg(2,filename2)
   CALL getarg(3,outname)
   CALL getarg(4,proc)
   CALL getarg(5,ptype)
   CALL getarg(6,experiment)
   CALL get_n_lines(filename1,lenfile(1))
   CALL get_n_lines(filename2,lenfile(2))
   ! The file containing the central values has one header
   ! line. The file with the uncertainties has two.
   IF (lenfile(1)-2.NE.lenfile(2)-1) THEN
      WRITE(*,*) "Did not find an equal number of central values and uncertainties."
      CALL exit(0)
   END IF
   nlines=lenfile(1)
   ! File containing the systematics has the following columns:
   ! 1: x_lower
   ! 2: x_upper
   ! 3: x_avg
   ! 4: Q2_lower
   ! 5: Q2_upper
   ! 6: Q2_avg
   ! 7: E_nu_lower
   ! 8: E_nu_upper
   ! 9: E_nu_avg
   !10: d2sigma/dxdQ2
   !11: N_events
   !12: N_events_errs = SQRT(N_events)
   !13: N_sys_errs
   !14: Percentage_error_theta
   !15: Percentage_error_Elepton
   !16: Percentage_error_Ehadron
   !17: MC_Samples
   k = (/ 3,6,9,10,12,13,14 /)

   ! The central values are found in column 
   ! 4: nu
   ! 5: nub

   ALLOCATE(dat(nlines-2,17))
   ALLOCATE(centr(nlines-2,5))
   ALLOCATE(fluctuated_central_vals(nlines-2))
   ALLOCATE(fluctuated_central_vals_os(nlines-2))
   ALLOCATE(central_vals(nlines-2))
   ALLOCATE(syserr(nlines-2))
   ALLOCATE(staterr(nlines-2))
   ALLOCATE(err(nlines-2))

   OPEN(10,FILE=filename1)
   READ(10,"(A)",END=20) line

   READ(10,"(A)",END=20) line
   DO i=1,UBOUND(dat,1)
      READ(10,"(A)",END=20) line
      DO j=LBOUND(dat,2),UBOUND(dat,2)
         CALL read_variable_from_line(line,dat(i,j))
      END DO
      ! Calculate y and replace Enu avg by it
      dat(i,9) = dat(i,6) / (2*dat(i,3)*M*dat(i,9))
   END DO
   20 CLOSE(10)
   OPEN(30,FILE=filename2)
   READ(30,"(A)",END=40) line
   DO i=1,UBOUND(centr,1)
      READ(30,"(A)",END=40) line
      DO j=1,UBOUND(centr,2)
         CALL read_variable_from_line(line,centr(i,j))
      END DO
   END DO
   40 CLOSE(30)
   IF (TRIM(proc).EQ."nu") THEN
      dat(:,10)=centr(:,4)
      IF(TRIM(ptype).EQ."charm") THEN
         procstr="DIS_SNU_C"
         pname=TRIM(experiment)//"NU_DSIGMA_CHARM"
      ELSE IF(TRIM(ptype).EQ."inclusive") THEN
         procstr="DIS_SNU"
         pname=TRIM(experiment)//"NU_DSIGMA_INCLUSIVE"
      END IF
   ELSE IF (TRIM(proc).EQ."nub") THEN
      dat(:,10)=centr(:,5)
      IF(TRIM(ptype).EQ."charm") THEN
         procstr="DIS_SNB_C"
         pname=TRIM(experiment)//"NB_DSIGMA_CHARM"
      ELSE IF(TRIM(ptype).EQ."inclusive") THEN
         procstr="DIS_SNB"
         pname=TRIM(experiment)//"NB_DSIGMA_INCLUSIVE"
      END IF
   ELSE IF (TRIM(proc).EQ."sum") THEN
      dat(:,10)=centr(:,4)+centr(:,5)
      procstr="DIS"
      IF (TRIM(ptype).EQ."charm") pname=TRIM(experiment)//"_SUM_DSIGMA_CHARM"
      IF (TRIM(ptype).EQ."inclusive") pname=TRIM(experiment)//"_SUM_DSIGMA_INCLUSIVE"
   ELSE
      WRITE(*,*) "Unknown process type, expect nu or nub."
      CALL exit(0)
   END IF
   OPEN(60,FILE=outname)
   nsys_err = 1
   staterr(:)=1.0_dp /dat(:,12)*dat(:,10)
   DO i=LBOUND(syserr,1),UBOUND(syserr,1)
      syserr(i)=DSQRT((dat(i,10)*dat(i,14))**2 &
      & +(dat(i,10)*dat(i,15))**2 &
      & +(dat(i,10)*dat(i,16))**2)
      err(i)=DSQRT(staterr(i)**2 +syserr(i)**2)
      fluctuated_central_vals(i) = fluctuate(dat(i,10),err(i))
      fluctuated_central_vals_os(i) = fluctuate(dat(i,10),staterr(i))
   END DO
   central_vals=dat(:,10)
   dat(:,14)=syserr(:)/dat(:,10)*100
   dat(:,10)=fluctuated_central_vals(:)
   dat(:,12)=staterr(:)
   dat(:,13)=syserr(:)
   WRITE(60,"(A,I3,I3)") TRIM(pname),nsys_err,UBOUND(dat,1)
   DO j=1,UBOUND(dat,1)
      WRITE(60,"(I3,X,A,7ES16.8)") j,TRIM(procstr),(dat(j,k(i)),i=1,UBOUND(k,1))
   END DO
   CLOSE(60)
   dat(:,10)=central_vals(:)
   dat(:,14)=staterr(:)/dat(:,10)*100
   dat(:,10)=fluctuated_central_vals(:)
   dat(:,12)=0.0_dp
   dat(:,13)=staterr(:)
   dat(:,10)=fluctuated_central_vals_os(:)
   OPEN(70,FILE=outname(1:LEN_TRIM(outname)-4)//"_OS.dat")
   WRITE(70,"(A,I3,I3)") TRIM(pname)//"_OS",nsys_err,UBOUND(dat,1)
   DO j=1,UBOUND(dat,1)
      !WRITE(70,"(I3,X,A,7ES16.8)") j,TRIM(procstr),(dat(j,k(i)),i=1,UBOUND(k,1)-2),0.0_dp,0.0_dp
      WRITE(70,"(I3,X,A,7ES16.8)") j,TRIM(procstr),(dat(j,k(i)),i=1,UBOUND(k,1))
   END DO
   CLOSE(70)

   DEALLOCATE(dat)

END PROGRAM parse2commondata

SUBROUTINE read_variable_from_line(line,var)

   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)

   CHARACTER(LEN=500),INTENT(INOUT) :: line
   CHARACTER(LEN=500) :: substr
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

FUNCTION fluctuate(loc,sc) RESULT(fluc)
   IMPLICIT NONE
   INTEGER,PARAMETER :: dp=SELECTED_REAL_KIND(14,200)
   REAL,PARAMETER :: pi=3.14159265
   REAL(KIND=dp) :: loc,sc,fluc,x,r1,r2,u1,u2

   CALL random_number(r1)
   CALL random_number(r2)
   u1=1-r1
   u2=1-r2
   x = DSQRT(-2*DLOG(u1))*DCOS(2*pi*u2)
   fluc=loc+sc*x
END FUNCTION fluctuate
