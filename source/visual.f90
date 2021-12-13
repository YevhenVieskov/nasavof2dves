!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUTVEC_bin %%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE OUTPUTVEC_bin (U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, &
                             WORK1, xlength,ylength, imax,jmax, vecfile)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Writes U, V, P, PSI, and ZETA into "vecfile" for visualization
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: xlength, ylength
      CHARACTER(LEN=30), INTENT(IN) :: vecfile
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, P, TEMP, PSI, HEAT
      REAL(RP), DIMENSION(1:,1:), INTENT(IN) :: ZETA
      REAL(RP), DIMENSION(1:,1:), INTENT(OUT) :: WORK1
!
!  local variables
!
      INCLUDE 'defs.h'
      INTEGER :: i,j

      OPEN (1,file=vecfile,status='UNKNOWN',form='UNFORMATTED') 
      WRITE(1) xlength, ylength
      WRITE(1) imax, jmax

      DO j = 1,jmax
         DO i = 1,imax
            IF ((IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j)<C_E)) THEN
               WORK1(i,j) = (U(i,j) + U(i-1,j)) / 2.0
            ELSE
               WORK1(i,j) = 0.0
            END IF
         END DO
      END DO
      WRITE(1) WORK1

      DO j = 1,jmax
         DO i = 1,imax
            IF ((IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j)<C_E)) THEN
               WORK1(i,j) = (V(i,j) + V(i,j-1)) / 2.0
            ELSE
               WORK1(i,j) = 0.0
            END IF
         END DO
      END DO
      WRITE(1) WORK1

      DO j = 1,jmax
         DO i = 1,imax
            IF ((IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j)<C_E)) THEN
               WORK1(i,j) = P(i,j)
            ELSE
               WORK1(i,j) = 0.0
            END IF
         END DO
      END DO
      WRITE(1) WORK1

      DO j = 1,jmax
         DO i = 1,imax
            IF ((IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j)<C_E)) THEN
               WORK1(i,j) = TEMP(i,j)
            ELSE
               WORK1(i,j) = -0.5
            END IF
         END DO
      END DO
      WRITE(1) WORK1

      WRITE(1) ZETA
      WRITE(1) PSI
      WRITE(1) HEAT

      CLOSE (1)
      RETURN

   END SUBROUTINE OUTPUTVEC_bin
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_PSI_ZETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_PSI_ZETA (U, V, PSI, ZETA, FLAG, imax,jmax,delx,dely)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Computation of stream function and vorticity
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delx, dely
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: PSI
      REAL(RP), DIMENSION(1:,1:), INTENT(OUT) :: ZETA
!
! local variables
!
      INCLUDE 'defs.h'
      INTEGER :: i, j
!
!  Computation of the vorticity zeta at the upper right corner 
!  of cell (i,j) (only if the corner is surrounded by fluid cells)
!
      DO i = 1,imax-1   
         DO j = 1,jmax-1

            IF(((IAND(FLAG(i,  j),C_F)/=0).AND.(FLAG(i,j)<C_E)).AND.& 
               ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E)).AND.&
               ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)< C_E)).AND.&
               ((IAND(FLAG(i+1,j+1),C_F)/=0).AND.(FLAG(i+1,j+1)<C_E)))&
               THEN

               ZETA(i,j) = (U(i,j+1)-U(i,j))/dely-(V(i+1,j)-V(i,j))/delx
 
            ELSE

               ZETA(i,j) = 0.0

            END IF

         END DO
      END DO
!
! Computation of the stream function at the upper right corner 
! of cell (i,j) (only if both lower cells are fluid cells) 
!
      DO i = 0,imax
         PSI(i,0) = 0.0
         DO j = 1,jmax

            IF (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,  j)<C_E)).OR. &
                ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E))) THEN
               PSI(i,j) = PSI(i,j-1) + U(i,j)*dely
            ELSE
               PSI(i,j) = PSI(i,j-1);
            END IF

         END DO
      END DO

      RETURN

   END SUBROUTINE COMP_PSI_ZETA
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_HEAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_HEAT (U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, &
                         delx, dely)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Computation of the heat function 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delx, dely, Re, Pr
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, TEMP
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: HEAT
!
! local variables
!
      INCLUDE 'defs.h'
      INTEGER :: i, j

      DO i = 0,imax
         HEAT(i,0) = 0.0
         DO j = 1,jmax
            IF  (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,j)  <C_E)).OR.&
                 ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E)))THEN
               HEAT(i,j) = HEAT(i,j-1) &
                   + dely*(U(i,j)*0.5*(1.0+TEMP(i+1,j)+TEMP(i,j))*Re*Pr&
                   -                      (TEMP(i+1,j)-TEMP(i,j))/delx )
            END IF       
         END DO
      END DO

      RETURN

   END SUBROUTINE COMP_HEAT
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET_PARTICLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE SET_PARTICLES (N, pos1x, pos1y, pos2x, pos2y, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Set initual coordinates where particles are injected. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      REAL(RP), INTENT(IN) :: pos1x, pos1y, pos2x, pos2y
      TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT) :: Partlines 
!
! local variables
!
      INTEGER :: i, status_var
      REAL(RP) :: hx, hy, x, y

      IF (N >= 2) THEN
         hx  = (pos2x-pos1x) / (N-1)
         hy  = (pos2y-pos1y) / (N-1)
      END IF

      DO i = 1,N
         x = pos1x + hx*(i-1)
         y = pos1y + hy*(i-1)
         Partlines(i)%Particles%x = x
         Partlines(i)%Particles%y = y
         Partlines(i)%length = 0
      END DO

      RETURN

   END SUBROUTINE SET_PARTICLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%% ADVANCE_PARTICLES %%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE ADVANCE_PARTICLES (imax, jmax, delx, dely, delt, &
                                 U, V, FLAG, N, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Advance particles by Euler method
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax, N
      REAL(RP), INTENT(IN) :: delx, dely, delt
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      TYPE(PARTICLELINE),DIMENSION(:),INTENT(INOUT),TARGET :: Partlines
!
! local variables
! 
      INCLUDE 'defs.h'
      INCLUDE 'interfaces.h'
      INTEGER :: i, j, k
      REAL(RP) :: x, y, x1, y1, x2, y2, uu, vv
      TYPE(PARTICLE), POINTER :: part, temp, help

      DO k = 1,N

       part => Partlines(k)%Particles 

       DO

         IF (.NOT.(ASSOCIATED(part%next))) EXIT
         temp => part%next
!
! advance all elements except the first element, which contains the injection
! point for new particles
!
         x = temp%x
         y = temp%y
!
! Computation of new x-coordinates by discretizing dx/dt=u 
!
         i = INT(x/delx)+1
         j = INT((y+0.5*dely)/dely)+1

         x1 = (i-1)*delx
         y1 = ((j-1)-0.5)*dely
         x2 = i*delx
         y2 = (j-0.5)*dely
!
! bilinear interpolation 
!
         uu= ((x2-x)*(y2-y)*U(i-1,j-1) &
           +  (x-x1)*(y2-y)*U(i,j-1)   &
           +  (x2-x)*(y-y1)*U(i-1,j)   &
           +  (x-x1)*(y-y1)*U(i,j))/delx/dely
!
! Computation of new y-coordinates by discretizing dy/dt=v 
!
         i = INT((x+0.5*delx)/delx)+1
         j = INT(y/dely)+1
 
         x1 = ((i-1)-0.5)*delx
         y1 = (j-1)*dely
         x2 = (i-0.5)*delx
         y2 = j*dely
!
! bilinear interpolation
!
         vv= ((x2-x)*(y2-y)*V(i-1,j-1) &
           +  (x-x1)*(y2-y)*V(i,j-1)   &
           +  (x2-x)*(y-y1)*V(i-1,j)   &
           +  (x-x1)*(y-y1)*V(i,j))/delx/dely
!
! velocity updates
!
         x = x + delt*uu
         y = y + delt*vv
!
! determine new cell for the particle
!
         i = INT(x/delx)+1
         j = INT(y/dely)+1
!
! if the particle exits the fluid domain, delete it
!
         IF ( (x >= imax*delx) .OR. (y >= jmax*dely) .OR. &
              (x <= 0.0)       .OR. (y <= 0.0) ) THEN
            IF (.NOT.(ASSOCIATED(temp%next))) THEN
               !
               ! temp is last particle
               !
               DEALLOCATE(temp)
               NULLIFY(part%next)
            ELSE
               !
               ! temp is NOT last particle
               !
               help => temp%next
               DEALLOCATE(temp)
               part%next => help
            END IF
            Partlines(k)%length = Partlines(k)%length - 1
         ELSE
!
! special treatment if particle would be in an inner obstacle cell
!
            IF (FLAG(i,j) < C_F) THEN
               CALL ADVANCE_AT_BOUND (i, j, x, y, uu, vv, U, V, FLAG, &
                                      delx, dely, delt)
            END IF

            temp%x = x
            temp%y = y

         END IF

         part => part%next

       END DO
      END DO

      RETURN

   END SUBROUTINE ADVANCE_PARTICLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%% ADVANCE_AT_BOUND %%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE ADVANCE_AT_BOUND (i, j, x, y, uu, vv, U, V, FLAG, &
                                delx, dely, delt)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Computation of new particle location of a particle near a no-slip
! wall, guaranteeing, that the new position is not in the obstacle
! cell.  Here a modIFied interpolation algorithm is applied, using the
! fact that at no-skip walls, the velocity is not only given at the
! midpoint of the edge but on the whole edge 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i, j
      REAL(RP), INTENT(IN) :: delx, dely, delt
      REAL(RP), INTENT(INOUT) :: x, y, uu, vv
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
!
! local variables
! 
      INCLUDE 'defs.h'
      INTEGER :: iold, jold
      REAL(RP) :: xold, yold
      REAL(RP) :: ul, ur, vo, vu      
      REAL(RP) :: x1, x2, y1, y2
!
! get old particle position
!
      xold = x - delt*uu
      yold = y - delt*vv
      iold = INT(xold/delx)+1
      jold = INT(yold/dely)+1
!
!-------------------- compute new (x,y) --------------------------------
!
      IF (i /= iold) THEN
         !
         !  compute new x
         !
         IF (FLAG(iold+1,jold) < C_F) THEN
            ur = 0.0
         ELSE
            IF (yold>= (jold-0.5)*dely) THEN
               IF (FLAG(iold+1,jold+1) < C_F) THEN
                  y2 = jold *dely
               ELSE
                  y1 = (jold-0.5)*dely
                  y2 = (jold+0.5)*dely
                  ur = (U(iold,jold)  *(y2-yold) &
                     +  U(iold,jold+1)*(yold-y1))/dely
               END IF
            ELSE      
               IF (FLAG(iold+1,jold-1) < C_F) THEN
                  y1 = (jold-1.0)*dely
                  ur = U(iold,jold)*(yold-y1)*2.0/dely
               ELSE
                  y1 = (jold-1.5)*dely
                  y2 = (jold-0.5)*dely
               END IF
            END IF
         END IF

         IF (FLAG(iold-1,jold) < C_F) THEN
            ul = 0.0  
         ELSE
            IF (yold >= (jold-0.5)*dely) THEN
               IF (FLAG(iold-1,jold+1) < C_F) THEN
                  y2 = jold *dely
                  ul = U(iold-1,jold)*(y2-yold)*2.0/dely
               ELSE   
                  y1 = (jold-0.5)*dely
                  y2 = (jold+0.5)*dely
                  ul = (U(iold-1,jold)  *(y2-yold) &
                     +  U(iold-1,jold+1)*(yold-y1))/dely
	       END IF 
            ELSE       
               IF (FLAG(iold-1,jold-1) < C_F) THEN
                  y1 = (jold-1.0)*dely
                  ul = U(iold-1,jold)*(yold-y1)*2.0/dely
               ELSE
                  y1 = (jold-1.5)*dely
                  y2 = (jold-0.5)*dely
                  ul = (U(iold-1,jold-1)*(y2-yold) &
                     +  U(iold-1,jold)  *(yold-y1))/dely
               END IF
	    END IF 
         END IF

         uu = (ul*(iold*delx-xold)+ur*(xold-(iold-1)*delx))/delx
         x = xold + uu*delt

      END IF  ! new x is finished

      IF (j /= jold) THEN
!
!  compute new y
!
         IF (FLAG(iold,jold+1) < C_F) THEN
            vo = 0.0   
         ELSE
            IF (xold >= (iold-0.5)*delx) THEN
               IF (FLAG(iold+1,jold+1) < C_F) THEN
                  x2 = iold*delx
                  vo = V(iold,jold)*(x2-xold)*2.0/delx
               ELSE  
                  x1 = (iold-0.5)*delx
                  x2 = (iold+0.5)*delx
                  vo = (V(iold,jold)  *(x2-xold) &
                     +  V(iold+1,jold)*(xold-x1))/delx
               END IF
            ELSE      
               IF (FLAG(iold-1,jold+1) < C_F) THEN
                  x1 = (iold-1.0)*delx
                  vo = V(iold,jold)*(xold-x1)*2.0/delx
               ELSE
                  x1 = (iold-1.5)*delx
                  x2 = (iold-0.5)*delx
                  vo = (V(iold-1,jold)*(x2-xold) &
                     +  V(iold,jold)  *(xold-x1))/delx
               END IF
            END IF
         END IF

         IF (FLAG(iold,jold-1) < C_F) THEN
            vu = 0.0  
         ELSE
            IF (xold>= (iold-0.5)*delx) THEN
               IF (FLAG(iold+1,jold-1) < C_F) THEN
                  x2 = iold*delx
                  vu = V(iold,jold-1)*(x2-xold)*2.0/delx
               ELSE   
                  x1 = (iold-0.5)*delx
                  x2 = (iold+0.5)*delx
                  vu = (V(iold,jold-1)  *(x2-xold) &
                     +  V(iold+1,jold-1)*(xold-x1))/delx
               END IF
            ELSE       
               IF (FLAG(iold-1,jold-1) < C_F) THEN
                  x1 = (iold-1.0)*delx
                  vu = V(iold,jold-1)*(xold-x1)*2.0/delx
               ELSE 
                  x1 = (iold-1.5)*delx
                  x2 = (iold-0.5)*delx
                  vu = (V(iold-1,jold-1)*(x2-xold) &
                     + V(iold,jold-1)*(xold-x1))/delx
               END IF
            END IF
         END IF

         vv = (vu*(jold*dely-yold)+vo*(yold-(jold-1)*dely))/dely
         y = yold +vv*delt
  
      END IF  ! new y is finished

      RETURN

  END SUBROUTINE ADVANCE_AT_BOUND
!
!%%%%%%%%%%%%%%%%%%%%%%%%%% WRITE_PARTICLES_ASCII %%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE WRITE_PARTICLES_ASCII (tracefile, itype, t, N, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Append particle positions to file "partfile" in ascii format
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      CHARACTER(LEN=30) :: tracefile
      INTEGER, INTENT(IN) :: N, itype
      REAL(RP), INTENT(IN) :: t
      TYPE(PARTICLELINE),DIMENSION(:),INTENT(IN),TARGET :: Partlines
!
! local variables
!
      TYPE(PARTICLE), POINTER :: part, temp
      INTEGER :: k, length, count 
      
      OPEN (1,file=tracefile,position='APPEND', &
           status='OLD',form='FORMATTED')
!
!  write time stamp only if streakfile (not particle trace)
!
      IF (itype == 1) WRITE (1,*) t, ' -1.0    -1    -1' 
      DO k = 1, N
         count = 0
         part => Partlines(k)%Particles
         length = Partlines(k)%length
         DO WHILE (length > 0) 
            count = count + 1
            part => part%next
            WRITE(1,'(1x,g12.6,1x,g12.6,1x,i4,1x,i4)') part%x, part%y, &
               count, length
            IF (.NOT.(ASSOCIATED(part%next))) EXIT
         END DO
      END DO

      CLOSE(1)
      RETURN

   END SUBROUTINE WRITE_PARTICLES_ASCII
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTICLE_TRACING %%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE PARTICLE_TRACING (tracefile, t, imax, jmax, &
                                delx, dely, delt, U, V, FLAG, &
                                N, Partlines, iwrite)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Move particle positions and append them to a file if desired
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE

      CHARACTER(LEN=30), INTENT(IN) :: tracefile
      INTEGER, INTENT(IN) :: imax, jmax, N, iwrite
      REAL(RP), INTENT(IN) :: t, delx, dely, delt
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT) :: Partlines

      INCLUDE 'interfaces.h'

      IF (t <= 0) THEN
         OPEN(1,file=tracefile,status='REPLACE',form='FORMATTED')        
         CLOSE(1)
         CALL WRITE_PARTICLES_ASCII (tracefile, 0, t, N, Partlines)
      END IF

      CALL ADVANCE_PARTICLES (imax, jmax, delx, dely, delt, &
                              U, V, FLAG, N, Partlines)

      IF ( IAND(iwrite,1) /= 0 ) THEN
         CALL WRITE_PARTICLES_ASCII (tracefile, 0, t+delt, N, Partlines)
      END IF
 
      RETURN 
    
   END SUBROUTINE PARTICLE_TRACING
!
!%%%%%%%%%%%%%%%%%%%%%%% INJECT_PARTICLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE INJECT_PARTICLES (N, Partlines)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Injection of new particles for streaklines
!
!  NOTE:  most recent particle is place at head of list, just after
!         the injection point
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      TYPE(PARTICLELINE), DIMENSION(:),INTENT(INOUT),TARGET :: Partlines
!
!  local variables
!
      INTEGER :: k, status_var
      TYPE(PARTICLE), POINTER :: part, predecessor

      DO k = 1,N
         ALLOCATE (part,STAT=status_var)
         IF (status_var == 0) THEN
            part%x = Partlines(k)%Particles%x
            part%y = Partlines(k)%Particles%y
            predecessor => Partlines(k)%Particles
            part%next => predecessor%next
            predecessor%next => part
            Partlines(k)%length = Partlines(k)%length + 1  
         ELSE
            WRITE (6,*) ' memory allocation failure in INJECT_PARTICLES'
            STOP 'inject'
         END IF
      END DO

      RETURN

   END SUBROUTINE INJECT_PARTICLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% STREAKLINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE STREAKLINES (streakfile, iwrite, imax, jmax, &
                          delx, dely, delt, t, U, V, FLAG, &
                          N, Partlines)
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Move particles for streaklines, inject and write particle positions
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
      IMPLICIT NONE
 
      CHARACTER(LEN=30), INTENT(IN) :: streakfile
      INTEGER, INTENT(IN) :: imax, jmax, iwrite, N
      REAL(RP), INTENT(IN) :: delx, dely, delt, t
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      TYPE(PARTICLELINE), DIMENSION(:), INTENT(OUT) :: Partlines
 
      INCLUDE 'interfaces.h'

      IF (t <= 0) THEN
         OPEN(1,file=streakfile,status='REPLACE',form='FORMATTED')
         CLOSE(1)
         CALL WRITE_PARTICLES_ASCII (streakfile, 1, t, N, Partlines)
      END IF

      CALL ADVANCE_PARTICLES (imax, jmax, delx, dely, delt, &
                              U, V, FLAG, N, Partlines)

      IF ( IAND(IWRITE,2) /=0 ) THEN
         CALL INJECT_PARTICLES (N, Partlines)
      END IF

      IF ( IAND(IWRITE,4) /= 0) THEN
        CALL WRITE_PARTICLES_ASCII (streakfile, 1, t+delt, N, Partlines)
      END IF

      RETURN

   END SUBROUTINE STREAKLINES
