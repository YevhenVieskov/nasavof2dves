SUBROUTINE VTK(IOpst,filename1,imax, jmax,delx,dely,FS,U,V,P)
use nrtype
implicit none
integer,intent(in):: IOpst
character(200),intent(in)::filename1
REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: FS,U,V,P
REAL(RP), INTENT(IN) ::delx,dely
INTEGER, INTENT(IN) :: imax, jmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(RP), DIMENSION(:,:), ALLOCATABLE :: UTEMP,VTEMP
REAL(RP), DIMENSION(:), ALLOCATABLE :: X,Y
integer i,j,ios
ALLOCATE (UTEMP(0:imax+1,0:jmax+1))
ALLOCATE (VTEMP(0:imax+1,0:jmax+1))
ALLOCATE (X(0:imax+1))
ALLOCATE (Y(0:jmax+1))
 DO i = 0,imax
  DO j = 0,jmax
    X(I) = i*delx
    Y(J) = j*dely
  END DO
 END DO

!header vtk file
open(IOpst,file=filename1,status='new',iostat=ios)
write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
write(IOpst,'(A,A,A)')   'vof result output'
write(IOpst,'(A)')      'ASCII'

!grid data

write(IOpst,'(A)')      'DATASET RECTILINEAR_GRID' 
write(IOpst,'(A,1X,I5,1X,I5,1X,I5)')      'DIMENSIONS',IMAX+1,JMAX+1,1
write(IOpst,'(A,1X,I5,1X,A)')      'X_COORDINATES',IMAX+1,'double'  
write(IOpst,*)   X(0:IMAX)  
write(IOpst,'(A,1X,I5,1X,A)')      'Y_COORDINATES',JMAX+1,'double' 
write(IOpst,*)   Y(0:JMAX)  
write(IOpst,'(A,1X,I5,1X,A)')      'Z_COORDINATES',1,'double'
write(IOpst,*)   0.0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
write(IOpst,'(A,1X,I5)')   'POINT_DATA', (IMAX+1)*(JMAX+1)
write(IOpst,'(A,1X,A,1X,A)') 'VECTORS', 'velocity_point', 'double'
utemp=0d0
vtemp=0d0

do i=0,IMAX
  do j=0,JMAX !!!
     utemp(i,j)=0.5d0*(u(i,j+1)+u(i,j))
  end do
end do

do i=0,IMAX !!!
  do j=0,JMAX
     vtemp(i,j)=0.5d0*(v(i+1,j)+v(i,j))
  end do
end do

do j=0,JMAX
  do i=0,IMAX 
    write(IOpst,'(f8.3,1X,f8.3,1X,f8.3)') utemp(i,j),vtemp(i,j),0.0  
  end do
end do

write(IOpst,'(A,1X,I5)')   'CELL_DATA', (IMAX)*(JMAX)

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'pressure', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=1,JMAX
  do i=1,IMAX
    write(IOpst,'(f20.3)') p(i,j)
  end do
end do

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'vof_fraction', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=1,JMAX
  do i=1,IMAX  
    write(IOpst,'(f8.3)') fs(i,j)
  end do
end do 

write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'divergence', 'double',1
write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
do j=1,JMAX
  do i=1,IMAX  
    write(IOpst,'(f8.3)') DIV(i,j)
  end do
end do 
!write(IOpst,'(A,1X,A,1X,A,1X,I1)') 'SCALARS', 'divergence', 'double',1
!write(IOpst,'(A,1X,A)') 'LOOKUP_TABLE', 'default'
!do j=2,JMAX
!  do i=2,IMAX  
!    write(IOpst,'(f8.3)') divs(i,j)
!  end do
!end do   
DEALLOCATE (UTEMP, VTEMP,X,Y)
END SUBROUTINE VTK