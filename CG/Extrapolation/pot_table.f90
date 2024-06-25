      PROGRAM potential

!     Kuntal Ghosh
!     Takes an extrapolated pair potential file and derives the
!     corresponding .table file

      IMPLICIT NONE

      INTEGER*8 :: ios,nlines,i,a,j
      REAL*8 :: aa,dx,term,au
      REAL*8, ALLOCATABLE :: u(:),f(:),d(:)

      OPEN (1, FILE='u_extrap.dat', STATUS='OLD')
      OPEN (2, FILE='forces_extrap.dat', STATUS='OLD')
      OPEN (3, FILE='1_1_bumper.table', STATUS='UNKNOWN')

      dx = 0.1d0
         
      nlines = 0
      read_loop: DO
         READ (1,*,IOSTAT=ios)
         IF (ios/=0) EXIT read_loop
         nlines = nlines + 1
      END DO read_loop

      ALLOCATE (u(nlines),f(nlines),d(nlines))

      REWIND (1)

      DO i = 1,nlines
         READ (1,*) d(i),u(i)
         READ (2,*) d(i),f(i)
      END DO

!     Generate .table file

      WRITE (3,'(A)') "# Header information on force file"
      WRITE (3,*)
      WRITE (3,'(A)') "1_1"
      WRITE (3,'(A1,I5,2X,A1,2F16.6)') "N",nlines-1,"R",d(2),d(nlines)
      WRITE (3,*)

      DO i = 1,nlines-1
         IF (d(i)<10.0d0) THEN
            WRITE (3,'(I5,3F16.8)') i,d(i+1),u(i+1),f(i+1)
         ELSE
            WRITE (3,'(I5,3F16.8)') i,d(i+1),0.0d0,0.0d0
         END IF
      END DO 

      DEALLOCATE (u,f,d)

      END PROGRAM potential
