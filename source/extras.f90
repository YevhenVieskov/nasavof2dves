!
!%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT_ONE_real %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE OUTPUT_ONE_real (FIELD, imin, imax, jmin, jmax)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Outputs a single real field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FIELD
!
! local variables
!
      INTEGER :: i,j
      CHARACTER(LEN=30) filename

      filename = 'contour_data'
      OPEN(1,file=filename,status='UNKNOWN',form='FORMATTED')
      DO j = jmin, jmax
         DO i = imin, imax
            write(1,'(1x,i5,1x,i5,1x,f12.6)') i, j, FIELD(i,j)
         END DO
      END DO
      CLOSE(1)

      RETURN

   END SUBROUTINE OUTPUT_ONE_real
!
!%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT_ONE_integer %%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE OUTPUT_ONE_integer (IFIELD, imin, imax, jmin, jmax, &
                                  filename)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Outputs a single integer field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
      INTEGER, DIMENSION(0:,0:), INTENT(INOUT) :: IFIELD
      CHARACTER(LEN=30) filename
!
! local variables
!
      INTEGER :: i,j

      OPEN(1,file=filename,status='UNKNOWN',form='FORMATTED')
      DO j = jmin, jmax
         DO i = imin, imax
            write(1,'(1x,i5,1x,i5,1x,i6)') i, j, IFIELD(i,j)
         END DO
      END DO
      CLOSE(1)

      RETURN

   END SUBROUTINE OUTPUT_ONE_integer
!
!%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT_THREE_real %%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE OUTPUT_THREE_real (U,V,P, imin,imax, jmin,jmax, filename)
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P
      CHARACTER(LEN=30) filename
!
! local variables
!
      INTEGER :: i,j

      OPEN(1,file=filename,status='UNKNOWN',form='FORMATTED')
      DO j = jmin, jmax
         DO i = imin, imax
            write(1,'(1x,i5,1x,i5,3(1x,g12.6))') i, j, &
                      U(i,j), V(i,j), P(i,j)
         END DO
      END DO
      CLOSE(1)

      RETURN

   END SUBROUTINE OUTPUT_THREE_real
