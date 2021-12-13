!
!%%%%%%%%%%%%%%%%%%%%%%%% READ_PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE READ_PARAMETERS (problem, &
                              inputfile, infile, outfile, &
                              xlength, ylength, imax, jmax, &
                              t_end, delt, tau, &
                              delx, dely, del_vec, &
                              del_trace, del_streak, del_inj, &
                              vecfile, tracefile, streakfile, &
                              N, pos1x, pos1y, pos2x, pos2y, &
                              itermax, eps, omega, gamma, p_bound, &
                              Re, Pr, beta, GX, GY, UI, VI, TI, &
                              wW, wE, wN, wS )
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  reads parameter data and returns to main program
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
     USE nrtype

     IMPLICIT NONE

     INTEGER, INTENT(OUT) :: imax, jmax, &
                             N, itermax, &
                             p_bound, &
                             wW, wE, wN, wS

     REAL(RP), INTENT(OUT) :: xlength, ylength, &
                              delx, dely, &
                              t_end, delt,tau, &
                              del_trace, del_inj, &
                              del_streak, del_vec, &
                              pos1x, pos1y, pos2x, pos2y, &
                              eps, omega, gamma, Re, Pr, beta, &
                              GX, GY, UI, VI, TI

     CHARACTER (LEN=30), INTENT(OUT) :: problem, &
                                        inputfile, &
                                        vecfile, tracefile, streakfile, &
                                        infile, outfile
!
!-----------------------------------------------------------------------
!
     inputfile="none" !main
     open(5,file="dam.par",STATUS="OLD") !!!main
     !open(6,file="backstep.rez",STATUS="new")
     read (5,*) problem
     read (5,'(1a1)') 

     read (5,*) infile
     read (5,*) outfile
     read (5,'(1a1)') 

     read (5,*) xlength
     read (5,*) ylength
     read (5,*) imax
     read (5,*) jmax
        delx = xlength / imax
        dely = ylength / jmax
     read (5,'(1a1)') 

     read (5,*) t_end
     read (5,*) delt
     read (5,*) tau
     read (5,'(1a1)') 

     read (5,*) del_trace
     read (5,*) del_inj
     read (5,*) del_streak
     read (5,*) del_vec
     read (5,'(1a1)') 
  
     read (5,*) vecfile
     read (5,*) tracefile
     read (5,*) streakfile
     read (5,'(1a1)') 

     read (5,*) N
     read (5,*) pos1x
     read (5,*) pos1y
     read (5,*) pos2x
     read (5,*) pos2y
     read (5,'(1a1)') 

     read (5,*) itermax
     read (5,*) eps 
     read (5,*) omega
     read (5,*) gamma
     read (5,*) p_bound
     read (5,'(1a1)') 
     read (5,'(1a1)') 
     read (5,'(1a1)') 

     read (5,*) Re
     read (5,*) Pr
     read (5,*) beta
     read (5,*) GX
     read (5,*) GY
     read (5,*) UI
     read (5,*) VI
     read (5,*) TI
     read (5,'(1a1)') 

     read (5,*) wW
     read (5,*) wE
     read (5,*) wN
     read (5,*) wS
!
!-----------------------------------------------------------------------
!
     write (6,*) " Problem: ", problem

     write (6,*) " xlength = ", xlength
     write (6,*) " ylength = ", ylength

     write (6,*) " imax = ", imax
     write (6,*) " jmax = ", jmax

     write (6,*) " delx = ", delx
     write (6,*) " dely = ", dely

     write (6,*) " delt = ", delt
     write (6,*) " t_end = ", t_end
     write (6,*) " tau = ", tau

     write (6,*) " del_trace = ", del_trace
     write (6,*) " del_inj = ", del_inj
     write (6,*) " del_streak = ", del_streak
     write (6,*) " del_vec = ", del_vec

     write (6,*) " vecfile: ", vecfile
     write (6,*) " tracefile: ", tracefile
     write (6,*) " streakfile: ", streakfile
     write (6,*) " infile: ", infile
     write (6,*) " outfile: ", outfile

     write (6,*) " N = ", N
     write (6,*) " pos1x = ", pos1x
     write (6,*) " pos1y = ", pos1y
     write (6,*) " pos2x = ", pos2x
     write (6,*) " pos2y = ", pos2y

     write (6,*) " itermax = ", itermax
     write (6,*) " eps = ", eps
     write (6,*) " omega = ", omega
     write (6,*) " gamma = ", gamma
     IF ((p_bound /= 1) .AND. (p_bound /= 2)) THEN
        write (6,*) " p_bound must be 1 or 2"
        stop ' init'
     ELSE
        write (6,*) " p_bound = ", p_bound
     END IF

     write (6,*) " Re = ", Re
     write (6,*) " Pr = ", Pr
     write (6,*) " beta = ", beta
     write (6,*) " GX = ", GX
     write (6,*) " GY = ", GY
     write (6,*) " UI = ", UI
     write (6,*) " VI = ", VI
     write (6,*) " TI = ", TI

     IF ((wW > 4).OR.(wW < 1)) THEN
        write (6,*) "wW must be 1,2,3, or 4"
        stop ' wW'
     ELSE
        write (6,*) " wW = ", wW
     END IF

     IF ((wE > 4).OR.(wE < 1)) THEN
        write (6,*) "wE must be 1,2,3, or 4"
        stop ' wE'
     ELSE
        write (6,*) " wE = ", wE
     END IF

     IF ((wN > 4).OR.(wN < 1)) THEN
        write (6,*) "wN must be 1,2,3, or 4"
        stop ' wN'
     ELSE
        write (6,*) " wN = ", wN
     END IF

     IF ((wS > 4).OR.(wS < 1)) THEN
        write (6,*) "wS must be 1,2,3, or 4"
        stop ' wS'
     ELSE
        write (6,*) " wS = ", wS
     END IF

     IF (((wW == 4).AND.(wE /= 4)) .OR. (wW /= 4).AND.(wE == 4)) THEN
        write (6,*) "Periodic boundary conditions need wW=wE=4"
        stop ' wW or wE'
     END IF

     IF (((wS == 4).AND.(wN /= 4)) .OR. (wS /= 4).AND.(wN == 4)) THEN
        write (6,*) "Periodic boundary conditions need wS=wN=4"
        stop ' wN or wS'
     END IF

     RETURN

   END SUBROUTINE READ_PARAMETERS
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT_UVP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE INIT_UVP (problem, U, V, P, TEMP, imax, jmax, UI, VI, TI)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  initializes U, V, P, and TEMP fields
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER :: i,j

      REAL(RP), INTENT(IN) :: UI, VI, TI
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P, TEMP

      CHARACTER (LEN=30), INTENT(IN) :: problem

      DO i = 0,imax+1
         DO j = 0,jmax+1

            U(i,j) = UI
            V(i,j) = VI
            P(i,j) = 0.
            TEMP(i,j) = TI

         END DO
      END DO

      IF (problem == "backstep") THEN

         DO i = 0,imax+1
            DO j = 0,jmax/2
               U(i,j) = 0.0
            END DO
         END DO

      END IF

      RETURN

   END SUBROUTINE INIT_UVP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ_BIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE READ_BIN (filename, U,V,P,TEMP,FLAG, imax,jmax, check) 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Binary read of fields U, V, P, and TEMP for use in restarts
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: check
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V, P, TEMP
      CHARACTER (LEN=30), INTENT(IN) :: filename
!
! local variables
!
      INTEGER :: i_in, j_in

      OPEN (1, file=filename, status='OLD', form='UNFORMATTED', &
               IOSTAT=check)

      IF (check == 0) then

         READ(1,IOSTAT=check) i_in, j_in
         IF ( (i_in /= imax) .OR. (j_in /= jmax) ) THEN
            WRITE (6,*) ' i_in = ', i_in, ' imax = ', imax
            WRITE (6,*) ' j_in = ', j_in, ' jmax = ', jmax
            WRITE (6,*) ' dimensional mismatch upon restart: STOP'
            STOP 'READ_BIN'
         END IF
         READ(1,IOSTAT=check) U
         READ(1,IOSTAT=check) V
         READ(1,IOSTAT=check) P
         READ(1,IOSTAT=check) TEMP
         READ(1,IOSTAT=check) FLAG

         if (check /= 0) check = 2
      
      ELSE

         check = 1

      END IF

      CLOSE(1)

      RETURN

   END SUBROUTINE READ_BIN
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%% WRITE_BIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE WRITE_BIN (filename, U, V, P, TEMP, FLAG, imax, jmax)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Binary write of fields U, V, P, and TEMP for use in restarts
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, P, TEMP
      CHARACTER (LEN=30), INTENT(IN) :: filename

      OPEN(1,file=filename,status='UNKNOWN',form='UNFORMATTED')

      WRITE(1) imax, jmax

      WRITE(1) U
      WRITE(1) V
      WRITE(1) P
      WRITE(1) TEMP
      WRITE(1) FLAG

      CLOSE(1)

      RETURN

   END SUBROUTINE WRITE_BIN
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INIT_FLAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE INIT_FLAG (problem, FLAG, imax, jmax, delx, dely, ibound,FS)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Defines and prints the geometry for the selected case
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: FLAG
      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: ibound
      REAL(RP), INTENT(IN) :: delx, dely
      CHARACTER (LEN=30), INTENT(IN) :: problem
      CHARACTER (LEN=1), DIMENSION(0:200) :: one_line
	  REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: FS
!	  
!  local variables
!
      INCLUDE 'defs.h'

      INTEGER :: i, j, low, up
      REAL(RP) :: mx, my, x, y, rad1

      IF (imax > 200) THEN
	  !IF (imax > 600) THEN
          write (6,*) ' imax too big in INIT_FLAG'
          stop ' imax'
      END IF
!
!  initialize boundary cells
!
      DO i = 0,imax+1
         FLAG(i,0) = C_B
         FLAG(i,jmax+1) = C_B
      END DO

      DO j = 1,jmax
         FLAG(0,j) = C_B
         FLAG(imax+1,j) = C_B
      END DO
!
!  initialize fluid cells
!
      DO i = 1,imax
         DO j = 1,jmax
            FLAG(i,j) = C_F
         END DO
      END DO
!
!  specialize to the problem
!
      IF (problem == "convection") THEN

         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "rayleigh") THEN

         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "dcavity") THEN

         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "dam") THEN

			DO J=0,JMAX+1
			  DO I=0,IMAX+1
			    FS(I,J)=0.0
			  END DO
			END DO
			
		DO i = 1,imax
         DO j = 1,jmax
            
			if(i<=20.and.j<=20)then
			!if(i<=40.and.j<=80)then
			!if(i*delx<=1.0.and.j*dely<=2.0) then
			  FS(I,J) =1.0
			  FLAG(I,J)=C_F
			else
              FLAG(I,J)=C_E
			end if
         END DO
        END DO


         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "drop") THEN

         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "fluidtrap") THEN

         DO i = 9*imax/22+1,13*imax/22
            DO j = 1,4*jmax/11
               FLAG(i,j) = C_B
            END DO
            DO j = 8*jmax/11+1,jmax
               FLAG(i,j) = C_B
            END DO
         END DO
         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "plate") THEN
         !
         ! flow past an inclined plate
         !
         low = 2*jmax/5
         up = 3*jmax/5
         FLAG(low,low) = C_B
         FLAG(low,low+1) = C_B
         FLAG(up,up-1) = C_B
         FLAG(up,up) = C_B
         DO i = low+1,up-1
            DO j = i-1,i+1
               FLAG(i,j) = C_B
            END DO
         END DO
         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "backstep" .OR. problem == "wave") THEN

         DO i = 1,jmax
            DO j = 1,jmax/2
               FLAG(i,j) = C_B
            END DO
         END DO
         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "circle") THEN

         mx = 20.0/41.0*jmax*dely
         my = mx
         rad1 = 5.0/41.0*jmax*dely
         DO i = 1,imax
            DO j = 1,jmax
               x = (i-0.5)*delx
               y = (j-0.5)*dely
               IF ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1) THEN
                  FLAG(i,j) = C_B
               END IF
            END DO
         END DO
         WRITE (6,*) ' problem = ', problem

      ELSE IF (problem == "molding") THEN
         !
         ! circular obstacle
         ! 
         mx = jmax*dely/2
         my = jmax*dely/2
         rad1 = jmax*dely/6
         DO i = 1,imax
            DO j = 1,jmax
               x = (i-0.5)*delx
               y = (j-0.5)*dely
               IF ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1) THEN
                  FLAG(i,j) = C_B
               END IF
            END DO
         END DO
         WRITE (6,*) ' problem = ', problem

      ELSE

         WRITE (6,*) ' undefined problem: STOP'
         STOP 'init'

      END IF
!
!  graphical output of geometry
!
      WRITE (6,*) ' '
      WRITE (6,*) ' Geometry of the fluid domain: '
      WRITE (6,*) ' '

 !     DO j = jmax+1,0,-1
!         DO i = 0,imax+1
 !           IF ( IAND(FLAG(i,j),C_F) /= C_F) THEN
 !              WRITE (one_line(i),'(a1)') '*'
!            ELSE
!               WRITE (one_line(i),'(a1)') ' '
!            END IF
!         END DO
!         WRITE (6,*) one_line(0:imax+1)
!      END DO
      WRITE (6,*) ' '
      !
      ! flags for boundary cells
      !
      ibound = 0
      DO i = 1, imax
         DO j = 1, jmax

            IF ( IAND(FLAG(i,j),C_F) /= C_F) THEN
               ibound = ibound + 1
            END IF

            FLAG(i,j) = FLAG(i,j) + ( IAND(FLAG(i-1,j),C_F) * B_W  &
                                  +   IAND(FLAG(i+1,j),C_F) * B_E  &
                                  +   IAND(FLAG(i,j-1),C_F) * B_S  &
                                  +   IAND(FLAG(i,j+1),C_F) * B_N ) / C_F

         END DO
      END DO

      WRITE (6,*) ' ibound = ', ibound
      WRITE (6,*) ' '
      DO j = jmax+1,0,-1
         WRITE (6,'(1x,200i2)') (FLAG(i,j), i=0,imax+1)
      END DO

      RETURN

  END SUBROUTINE INIT_FLAG
