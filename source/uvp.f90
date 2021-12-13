!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_TEMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_TEMP (U, V, TEMP, FLAG, WORK, imax, jmax, delt, &
                         delx, dely, gamma, Re, Pr)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes temperature field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delt, delx, dely, gamma, Re, Pr
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: WORK
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: TEMP
!
!  local variables
!
      INCLUDE 'defs.h'

      REAL(RP) :: laplt, dutdx, dvtdy, indelx2, indely2
      INTEGER :: i, j
      
      indelx2 = 1. / (delx * delx)
      indely2 = 1. / (dely * dely)

      DO i = 1, imax
         DO j = 1, jmax

            IF ( (IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j) < C_E) ) THEN

               laplt = (TEMP(i+1,j)-2.*TEMP(i,j)+TEMP(i-1,j))*indelx2 &
                     + (TEMP(i,j+1)-2.*TEMP(i,j)+TEMP(i,j-1))*indely2

               dutdx = ( ( U(i,j)  *0.5*( TEMP(i,j)   + TEMP(i+1,j) ) &
                     -     U(i-1,j)*0.5*( TEMP(i-1,j) + TEMP(i,j)   ))&
                   + gamma*(ABS(U(i,j))  *0.5*(TEMP(i,j)  -TEMP(i+1,j)) &
                   -        ABS(U(i-1,j))*0.5*(TEMP(i-1,j)-TEMP(i,j)))&
                     ) / delx

               dvtdy = ( ( V(i,j)  *0.5*( TEMP(i,j)   + TEMP(i,j+1) ) &
                     -     V(i,j-1)*0.5*( TEMP(i,j-1) + TEMP(i,j)   ))&
                   + gamma*(ABS(V(i,j))  *0.5*(TEMP(i,j)  -TEMP(i,j+1)) &
                   -        ABS(V(i,j-1))*0.5*(TEMP(i,j-1)-TEMP(i,j)))&
                     ) / dely

               WORK(i,j) = TEMP(i,j) + delt*(laplt/Re/Pr -dutdx -dvtdy)

            END IF

         END DO
      END DO

      DO i = 1, imax
         DO j = 1, jmax

            IF ( (IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j) < C_E) ) THEN

               TEMP(i,j) = WORK(i,j)

            END IF

         END DO
      END DO

      RETURN

   END SUBROUTINE COMP_TEMP
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_FG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_FG (U, V, TEMP, F, G, FLAG, imax, jmax, & 
                       delt, delx, dely, GX, GY, gamma, Re, beta)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes F, G fields
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delt, delx, dely, GX, GY, gamma, Re, beta
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, TEMP
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: F, G
      !
      !  local variables
      !
      INCLUDE 'defs.h'

      INTEGER :: i, j
      REAL(RP) :: DU2DX, DUVDY, DUVDX, DV2DY, LAPLU, LAPLV
!
!  compute flux field F
!
      DO i = 1, imax-1
         DO j = 1, jmax
            ! 
            !  only if both adjacent cells are fluid cells
            !
            IF (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,j)  <C_E)).AND.&
                ((IAND(FLAG(i+1,j),C_F)/=0).AND.(FLAG(i+1,j)<C_E))) THEN

               DU2DX = ((U(i,j)+U(i+1,j))*(U(i,j)+U(i+1,j)) &
                     + gamma*ABS(U(i,j)+U(i+1,j))*(U(i,j)-U(i+1,j)) &
                     - (U(i-1,j)+U(i,j))*(U(i-1,j)+U(i,j)) &
                     - gamma*ABS(U(i-1,j)+U(i,j))*(U(i-1,j)-U(i,j))) &
                     / (4.0*delx)

               DUVDY = ((V(i,j)+V(i+1,j))*(U(i,j)+U(i,j+1)) &
                     + gamma*ABS(V(i,j)+V(i+1,j))*(U(i,j)-U(i,j+1)) &
                     - (V(i,j-1)+V(i+1,j-1))*(U(i,j-1)+U(i,j)) &
                     - gamma*ABS(V(i,j-1)+V(i+1,j-1))*(U(i,j-1)-U(i,j))) &
                     / (4.0*dely)

               LAPLU = (U(i+1,j)-2.0*U(i,j)+U(i-1,j))/delx/delx &
                     + (U(i,j+1)-2.0*U(i,j)+U(i,j-1))/dely/dely
   
               F(i,j) = U(i,j) + delt * (LAPLU/Re - DU2DX - DUVDY + GX) &
		      - delt*beta*GX*(TEMP(i,j)+TEMP(i+1,j))/2

            ELSE

               F(i,j) = U(i,j)

            END IF

         END DO
      END DO
!
!  compute flux field G
!
      DO i = 1, imax
         DO j = 1, jmax-1
            ! 
            !  only if both adjacent cells are fluid cells
            !
            IF (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,j)<C_E)).AND.&
                ((IAND(FLAG(i,j+1),C_F)/=0).AND.(FLAG(i,j+1)<C_E)))THEN

               DUVDX = ((U(i,j)+U(i,j+1))*(V(i,j)+V(i+1,j)) &
                     + gamma*ABS(U(i,j)+U(i,j+1))*(V(i,j)-V(i+1,j)) &
                     - (U(i-1,j)+U(i-1,j+1))*(V(i-1,j)+V(i,j)) &
                     - gamma*ABS(U(i-1,j)+U(i-1,j+1))*(V(i-1,j)-V(i,j))) &
	             / (4.0*delx)

               DV2DY = ((V(i,j)+V(i,j+1))*(V(i,j)+V(i,j+1)) &
                     + gamma*ABS(V(i,j)+V(i,j+1))*(V(i,j)-V(i,j+1)) &
                     - (V(i,j-1)+V(i,j))*(V(i,j-1)+V(i,j)) &
                     - gamma*ABS(V(i,j-1)+V(i,j))*(V(i,j-1)-V(i,j))) &
	             / (4.0*dely)

               LAPLV = (V(i+1,j)-2.0*V(i,j)+V(i-1,j))/delx/delx &
                     + (V(i,j+1)-2.0*V(i,j)+V(i,j-1))/dely/dely

               G(i,j) = V(i,j)+delt*(LAPLV/Re-DUVDX-DV2DY+GY) &
                      - delt*beta*GY*(TEMP(i,j)+TEMP(i,j+1))/2    

            ELSE
 
               G(i,j) = V(i,j)

            END IF

         END DO
      END DO
!
! F and G at external boundary
!
      DO j = 1, jmax
         F(0,j)    = U(0,j)
         F(imax,j) = U(imax,j)
      END DO
 
      DO i = 1, imax
         G(i,0)    = V(i,0)
         G(i,jmax) = V(i,jmax)
      END DO

      RETURN

   END SUBROUTINE COMP_FG
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_RHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_RHS (F, G, RHS, FLAG, imax, jmax, delt, delx, dely)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  computes right-hand side field
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delt, delx, dely
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: RHS
!
!     local variables
!
      INCLUDE 'defs.h'

      INTEGER :: i,j

      DO i = 1,imax
         DO j = 1,jmax
!
!           only for fluid and non-surface cells
!
            IF ((IAND(FLAG(i,j),C_F)/=0).AND.(FLAG(i,j)<C_O)) THEN  
               RHS(i,j) = ( ( F(i,j) - F(i-1,j) ) / delx &
                        +   ( G(i,j) - G(i,j-1) ) / dely) / delt
            END IF

         END DO
      END DO

      RETURN

   END SUBROUTINE COMP_RHS
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POISSON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE POISSON (P, RHS, FLAG, imax, jmax, delx, dely, &
                       eps, iter, itermax, omega, res, ifull, p_bound)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  iteration for Poisson's equation
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype
     
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax, itermax, ifull, p_bound 
      INTEGER, INTENT(INOUT) :: iter
      REAL(RP), INTENT(IN) :: delx, dely, eps, omega
      REAL(RP), INTENT(OUT) :: res
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: RHS
!
!     local variables
!
      INCLUDE 'defs.h'

      INTEGER i, j    !, iter
      INTEGER eps_E, eps_W, eps_N, eps_S
      REAL(RP) :: rdx2, rdy2
      REAL(RP) :: add, beta_2, beta_mod
      REAL(RP) :: p0

      rdx2 = 1./delx/delx
      rdy2 = 1./dely/dely
      beta_2 = -omega/(2.0*(rdx2+rdy2))

      p0 = 0.0
      DO i = 1, imax
         DO j = 1, jmax

            IF ( IAND(FLAG(i,j),C_F) /= 0 ) THEN
               p0 = p0 + P(i,j)*P(i,j)
            END IF

         END DO
      END DO

      p0 = sqrt(p0/ifull)
      IF (p0 < 0.0001) THEN
         p0 = 1.0
      END IF
!
!  SOR iteration
!
      DO iter = 1, itermax
        !
        ! p_bound is the boundary-condition type for pressure
        ! 
        IF (p_bound == 1) THEN  ! modify the equation at the boundary
          !
          ! relaxation for fluid cells 
          !
          DO i = 1, imax         !for (i=1;i<=imax;i+=1)
            DO j = 1, jmax       !for (j=1;j<=jmax;j+=1)   
              !
              ! five point star for interior fluid cells
              !
              IF (FLAG(i,j) == C_A) THEN 

                P(i,j) = (1.-omega)*P(i,j) & 
                       - beta_2*( (P(i+1,j) + P(i-1,j) )*rdx2 &
                       +          (P(i,j+1) + P(i,j-1) )*rdy2 &
                       - RHS(i,j) )
              !
              ! modified star near boundary 
              !
              ELSE IF((IAND(FLAG(i,j),C_F)/=0).AND.(FLAG(i,j)<C_O)) THEN

                INCLUDE 'eps.h'
                beta_mod = -omega / ( &
                         + ( eps_E + eps_W) * rdx2 &
                         + ( eps_N + eps_S) * rdy2 &
                         )
                P(i,j) = (1. - omega) * P(i,j) &
                       - beta_mod * ( &
                       ( eps_E*P(i+1,j) &
                       + eps_W*P(i-1,j) &
                       ) * rdx2 &
                     + ( eps_N*P(i,j+1) &
                       + eps_S*P(i,j-1) &
                       ) * rdy2 &
                       - RHS(i,j) )

              END IF

            END DO
          END DO
!
!  compute the residual
!
          res = 0.0
          DO i = 1,imax
            DO j = 1,jmax
               !
               ! only fluid cells
               !
               INCLUDE 'eps.h'
               IF ((IAND(FLAG(i,j),C_F)/=0) .AND. (FLAG(i,j)<C_O )) THEN  
                  add = ( eps_E*( P(i+1,j) - P(i,j) ) &
                      -   eps_W*( P(i,j)   - P(i-1,j) ) ) * rdx2 &
                      + ( eps_N*( P(i,j+1) - P(i,j) ) &
                      -   eps_S*( P(i,j)   - P(i,j-1) ) ) * rdy2 &
                      -  RHS(i,j)
                  res = res + add*add
 	       END IF 

            END DO
         END DO

         res = sqrt(res/ifull)/p0
         IF (res < eps) RETURN

        ELSE IF (p_bound == 2) THEN
          !
          ! copy values at external boundary   
          !
          DO i = 1,imax
            P(i,0)      = P(i,1)
	    P(i,jmax+1) = P(i,jmax)
          END DO
          DO j = 1,jmax
	     P(0,j)      = P(1,j)
	     P(imax+1,j) = P(imax,j)
          END DO
          !
          ! and at interior boundary cells
          !
          DO i = 1,imax
            DO j = 1,jmax

	      IF ( (FLAG(i,j) >= B_N) .AND. (FLAG(i,j) <= B_SE) ) THEN
	        IF ( FLAG(i,j) == B_N ) THEN
                  P(i,j) = P(i,j+1)
                ELSE IF ( FLAG(i,j) == B_E ) THEN
                  P(i,j) = P(i+1,j)
                ELSE IF ( FLAG(i,j) == B_S ) THEN
                  P(i,j) = P(i,j-1)
                ELSE IF ( FLAG(i,j) == B_W ) THEN
                  P(i,j) = P(i-1,j)
                ELSE IF ( FLAG(i,j) == B_NE ) THEN
                  P(i,j) = 0.5*(P(i,j+1) + P(i+1,j))
                ELSE IF ( FLAG(i,j) == B_SE ) THEN
                  P(i,j) = 0.5*(P(i,j-1) + P(i+1,j))
                ELSE IF ( FLAG(i,j) == B_SW ) THEN
                  P(i,j) = 0.5*(P(i,j-1) + P(i-1,j))
                ELSE IF ( FLAG(i,j) == B_NW ) THEN
                  P(i,j) = 0.5*(P(i,j+1) + P(i-1,j))
                END IF
              END IF

            END DO
          END DO
          !
          !  relaxation method for fluid cells
          !
          DO i = 1,imax 
            DO j = 1,jmax

	      IF ((IAND(FLAG(i,j),C_F)/=0).AND.(FLAG(i,j)<C_O)) THEN
	        P(i,j) = (1. - omega) * P(i,j) &
                       - beta_2 * ( (P(i+1,j) + P(i-1,j)) * rdx2 &
                       +            (P(i,j+1) + P(i,j-1)) * rdy2 &
                       - RHS(i,j) )
              END IF
       
            END DO
          END DO
          !
          ! computation of residual 
          !
          res = 0.0

          DO i = 1,imax
            DO j = 1,jmax
              !
              ! only fluid cells
              !
	      IF ((IAND(FLAG(i,j),C_F)/=0).AND.(FLAG(i,j)< C_O)) THEN
                add = (P(i+1,j)-2.*P(i,j)+P(i-1,j))*rdx2 &
                    + (P(i,j+1)-2.*P(i,j)+P(i,j-1))*rdy2 &
                    - RHS(i,j)
                res = res + add*add
              END IF

            END DO
          END DO
	
          res = sqrt(res/ifull)/p0
          IF (res < eps) RETURN
      
        END IF
      
      END DO

      RETURN

   END SUBROUTINE POISSON
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADAP_UV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE ADAP_UV (U, V, F, G, P, FLAG, imax, jmax, delt,delx,dely)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  compute the U, V fields if adjacent cells are fluid cells
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: imax, jmax
      REAL(RP), INTENT(IN) :: delt, delx, dely
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G, P
      REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V
!
!     local variables
!
      INCLUDE 'defs.h'

      INTEGER :: i, j

      DO i = 1, imax - 1
         DO j = 1, jmax
            !
            ! only if both adjacent cells are fluid cells
            !
            IF (((IAND(FLAG(i,j),  C_F)>0).AND.(FLAG(i,  j)<C_E)).AND.&
                ((IAND(FLAG(i+1,j),C_F)>0).AND.(FLAG(i+1,j)<C_E))) THEN

               U(i,j) = F(i,j) - ( P(i+1,j) - P(i,j) )*delt/delx

            END IF

         END DO
      END DO

      DO i = 1, imax
         DO j = 1, jmax - 1
            !
            ! only if both adjacent cells are fluid cells
            !
            IF (((IAND(FLAG(i,j),  C_F)>0).AND.(FLAG(i,  j)<C_E)).AND.&
                ((IAND(FLAG(i,j+1),C_F)>0).AND.(FLAG(i,j+1)<C_E))) THEN

               V(i,j) = G(i,j) - ( P(i,j+1) - P(i,j) )*delt/dely

            END IF

         END DO
      END DO


 DO j = 1,jmax 
  DO i = 1,imax
  IF(FLAG(I,J)>C_B.AND.FLAG(I,J)>C_E) THEN
    DIV(I,J)=(U(I,J)-U(I-1,J))/DX+(V(I,J)-V(I,J-1))/DY
  ELSE
    DIV(I,J)=0.0
  END IF
  
  END DO
 END DO

      RETURN

   END SUBROUTINE ADAP_UV
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMP_delt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
   SUBROUTINE COMP_delt (delt, t, imax, jmax, delx, dely, &
                         U, V, Re, Pr, tau, iwrite, &
                         del_trace, del_inj, del_streak, del_vec)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Computation of adaptive time stepsize satisfying 
! the CFL stability criteria and set the flag "write" if some data
! has to be written into a file. 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
      USE nrtype

      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: iwrite
      REAL(RP), INTENT(IN) :: delx, dely, Re, Pr, t, tau, &
                              del_trace, del_inj, del_streak, del_vec
      REAL(RP), INTENT(OUT) :: delt
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
!
!     local variables
!
      INTEGER :: i, j
      REAL(RP) :: umax, vmax, deltu, deltv, deltRePr, t_new
!
! delt satisfying CFL conditions
!
      IF ( tau >= 1.0e-10 ) THEN ! else no time stepsize control

         umax = 1.0e-10
         vmax = 1.0e-10 

         DO i = 0,imax+1
            DO j = 1,jmax+1
               IF (ABS(U(i,j)) > umax) umax = ABS(U(i,j))
            END DO
         END DO

         DO i = 1,imax+1
            DO j = 0,jmax+1
               IF (ABS(V(i,j)) > vmax) vmax = ABS(V(i,j))
            END DO
         END DO

         deltu = delx / umax
         deltv = dely / vmax 

         IF (Pr < 1) THEN
            deltRePr = 1./(1./(delx*delx)+1./(dely*dely))*Re*Pr/2.
         ELSE
            deltRePr = 1./(1./(delx*delx)+1./(dely*dely))*Re/2.
         END IF

         IF (deltu < deltv) THEN
            IF (deltu < deltRePr) THEN
               delt = deltu
            ELSE
               delt = deltRePr
            END IF
         ELSE
            IF (deltv < deltRePr) THEN
               delt = deltv
            ELSE
               delt = deltRePr
            END IF
         END IF

         delt = tau * delt   ! multiply by safety factor

      END IF
!
! look if some data has to be written to a file in the next time step
!
      iwrite = 0
      t_new = t + delt

      IF ( INT(t/del_trace) /= INT(t_new/del_trace) ) THEN
         iwrite = iwrite + 1
      END IF

      IF ( INT(t/del_inj) /= INT(t_new/del_inj) ) THEN
         iwrite = iwrite + 2
      END IF

      IF ( INT(t/del_streak) /= INT(t_new/del_streak) ) THEN
         iwrite = iwrite + 4
      END IF

      IF ( INT(t/del_vec) /= INT(t_new/del_vec) ) THEN
         iwrite = iwrite + 8
      END IF

      RETURN

   END SUBROUTINE COMP_delt
