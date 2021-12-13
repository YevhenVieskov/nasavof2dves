!
! init.f90 interfaces
!
   INTERFACE
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
      END SUBROUTINE READ_PARAMETERS
   END INTERFACE

   INTERFACE
      SUBROUTINE INIT_UVP (problem, U,V,P,TEMP, imax, jmax, UI, VI, TI)
         USE nrtype
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: UI, VI, TI
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P, TEMP
         CHARACTER (LEN=30), INTENT(IN) :: problem
      END SUBROUTINE INIT_UVP
   END INTERFACE

   INTERFACE
      SUBROUTINE READ_BIN (filename, U,V,P,TEMP,FLAG, imax,jmax, check)
         USE nrtype
         INTEGER, INTENT(IN) :: imax, jmax
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
         INTEGER, INTENT(OUT) :: check
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V, P, TEMP
         CHARACTER (LEN=30), INTENT(IN) :: filename
      END SUBROUTINE READ_BIN
   END INTERFACE

   INTERFACE
      SUBROUTINE WRITE_BIN (filename, U, V, P, TEMP, FLAG, imax, jmax)
         USE nrtype
         INTEGER, INTENT(IN) :: imax, jmax
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, P, TEMP
         CHARACTER (LEN=30), INTENT(IN) :: filename
      END SUBROUTINE WRITE_BIN
   END INTERFACE

   INTERFACE
   SUBROUTINE INIT_FLAG (problem, FLAG, imax, jmax, delx, dely, ibound,FS)

      USE nrtype

      IMPLICIT NONE
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: FLAG
      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: ibound
      REAL(RP), INTENT(IN) :: delx, dely
      CHARACTER (LEN=30), INTENT(IN) :: problem
      CHARACTER (LEN=1), DIMENSION(0:200) :: one_line
	  REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: FS
   END SUBROUTINE INIT_FLAG
 END INTERFACE
!
!  uvp.f90 interfaces
!   
   INTERFACE
      SUBROUTINE COMP_TEMP (U, V, TEMP, FLAG, WORK, imax, jmax, delt, &
                            delx, dely, gamma, Re, Pr)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delt, delx, dely, gamma, Re, Pr
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: TEMP
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: WORK
      END SUBROUTINE COMP_TEMP
   END INTERFACE

   INTERFACE
      SUBROUTINE COMP_FG (U, V, TEMP, F, G, FLAG, imax, jmax, & 
                          delt, delx, dely, GX, GY, gamma, Re, beta)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delt, delx, dely, GX, GY, gamma,Re,beta
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, TEMP
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: F, G
      END SUBROUTINE COMP_FG
   END INTERFACE

   INTERFACE
      SUBROUTINE COMP_RHS (F, G, RHS, FLAG, imax, jmax, delt,delx,dely)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delt, delx, dely
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: RHS
      END SUBROUTINE COMP_RHS
   END INTERFACE

   INTERFACE
      SUBROUTINE POISSON (P, RHS, FLAG, imax, jmax, delx, dely, &
                       eps, iter, itermax, omega, res, ifull, p_bound)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax, itermax, ifull, p_bound 
         INTEGER, INTENT(INOUT) :: iter
         REAL(RP), INTENT(IN) :: delx, dely, eps, omega
         REAL(RP), INTENT(OUT) :: res
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: P
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: RHS
      END SUBROUTINE POISSON
   END INTERFACE

   INTERFACE
      SUBROUTINE ADAP_UV (U,V, F,G, P, FLAG, imax, jmax, delt,delx,dely)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delt, delx, dely
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: F, G, P
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V
      END SUBROUTINE ADAP_UV
   END INTERFACE

   INTERFACE
      SUBROUTINE COMP_delt (delt, t, imax, jmax, delx, dely, &
                            U, V, Re, Pr, tau, iwrite, &
                            del_trace, del_inj, del_streak, del_vec)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         INTEGER, INTENT(OUT) :: iwrite
         REAL(RP), INTENT(IN) :: delx, dely, Re, Pr, t, tau, &
                                 del_trace, del_inj, del_streak, del_vec
         REAL(RP), INTENT(OUT) :: delt
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      END SUBROUTINE COMP_delt
   END INTERFACE
!
!  boundary.f90 interfaces
!
   INTERFACE
      SUBROUTINE SETBCOND (U, V, P, TEMP, FLAG, imax, jmax, wW,wE,wN,wS,FS)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax, wW, wE, wN, wS
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P, TEMP
		 REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
      END SUBROUTINE SETBCOND
   END INTERFACE

   INTERFACE
      SUBROUTINE SETSPECBCOND (problem, U, V, TEMP, imax, jmax, UI, VI)
         USE nrtype
         IMPLICIT NONE
         CHARACTER (LEN=30) :: problem
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, TEMP
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: UI, VI
      END SUBROUTINE SETSPECBCOND
   END INTERFACE
!
! intefaces for visual.f90
!
   INTERFACE
      SUBROUTINE OUTPUTVEC_bin (U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, &
                                WORK1, xlength,ylength, imax,jmax, vecfile)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: xlength, ylength
         CHARACTER(LEN=30), INTENT(IN) :: vecfile
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, P, TEMP, PSI, HEAT
         REAL(RP), DIMENSION(1:,1:), INTENT(IN) :: ZETA
         REAL(RP), DIMENSION(1:,1:), INTENT(OUT) :: WORK1
      END SUBROUTINE OUTPUTVEC_bin
   END INTERFACE

   INTERFACE
      SUBROUTINE COMP_PSI_ZETA (U,V,PSI,ZETA,FLAG,imax,jmax,delx,dely)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delx, dely
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: PSI
         REAL(RP), DIMENSION(1:,1:), INTENT(OUT) :: ZETA
      END SUBROUTINE COMP_PSI_ZETA
   END INTERFACE

   INTERFACE
      SUBROUTINE COMP_HEAT (U, V, TEMP, HEAT, FLAG, Re, Pr, &
                            imax, jmax, delx, dely)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: delx, dely, Re, Pr
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V, TEMP
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: HEAT
      END SUBROUTINE COMP_HEAT
   END INTERFACE

   INTERFACE
      SUBROUTINE SET_PARTICLES (N, pos1x,pos1y,pos2x,pos2y, Partlines)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N
         REAL(RP), INTENT(IN) :: pos1x, pos1y, pos2x, pos2y
         TYPE (PARTICLELINE), DIMENSION(:), INTENT(INOUT) :: Partlines 
      END SUBROUTINE SET_PARTICLES
   END INTERFACE

   INTERFACE
      SUBROUTINE ADVANCE_PARTICLES (imax, jmax, delx, dely, delt, &
                                    U, V, FLAG, N, Partlines)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax, N
         REAL(RP), INTENT(IN) :: delx, dely, delt
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
         TYPE(PARTICLELINE),DIMENSION(:),INTENT(INOUT), &
                            TARGET :: Partlines
      END SUBROUTINE ADVANCE_PARTICLES
   END INTERFACE

   INTERFACE
      SUBROUTINE WRITE_PARTICLES_ASCII (tracefile,itype,t,N,Partlines)
         USE nrtype
         IMPLICIT NONE
         CHARACTER(LEN=30) :: tracefile
         INTEGER, INTENT(IN) :: N, itype
         REAL(RP), INTENT(IN) :: t
         TYPE(PARTICLELINE),DIMENSION(:),INTENT(IN),TARGET :: Partlines
      END SUBROUTINE WRITE_PARTICLES_ASCII
   END INTERFACE

   INTERFACE
      SUBROUTINE PARTICLE_TRACING (tracefile, t, imax, jmax, &
                                   delx, dely, delt, U, V, FLAG, &
                                   N, Partlines, iwrite)
         USE nrtype
         IMPLICIT NONE
         CHARACTER(LEN=30), INTENT(IN) :: tracefile
         INTEGER, INTENT(IN) :: imax, jmax, N, iwrite
         REAL(RP), INTENT(IN) :: t, delx, dely, delt
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
         TYPE(PARTICLELINE), DIMENSION(:), INTENT(IN) :: Partlines
      END SUBROUTINE PARTICLE_TRACING
   END INTERFACE

   INTERFACE
      SUBROUTINE ADVANCE_AT_BOUND (i, j, x, y, uu, vv, U, V, FLAG, &
                                   delx, dely, delt)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: i, j
         REAL(RP), INTENT(IN) :: delx, dely, delt
         REAL(RP), INTENT(INOUT) :: x, y, uu, vv
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      END SUBROUTINE ADVANCE_AT_BOUND
   END INTERFACE

   INTERFACE
      SUBROUTINE INJECT_PARTICLES (N, Partlines)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: N
         TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT),TARGET :: &
             Partlines
      END SUBROUTINE INJECT_PARTICLES


	  
   END INTERFACE

   INTERFACE
      SUBROUTINE STREAKLINES (streakfile, iwrite, imax, jmax, &
                             delx, dely, delt, t, U, V, FLAG, &
                             N, Partlines)
         USE nrtype
         IMPLICIT NONE
         CHARACTER(LEN=30), INTENT(IN) :: streakfile
         INTEGER, INTENT(IN) :: imax, jmax, iwrite, N
         REAL(RP), INTENT(IN) :: delx, dely, delt, t
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
         TYPE(PARTICLELINE),DIMENSION(:),INTENT(OUT),TARGET :: Partlines
      END SUBROUTINE STREAKLINES
   END INTERFACE
!
! intefaces for surface.f90
!
   INTERFACE
      SUBROUTINE INIT_PARTICLES (N, imax, jmax, delx, dely, ppc, problem, &
                                 U, V, Partlines)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: N
         INTEGER, INTENT(IN) :: imax, jmax, ppc
         REAL(RP), INTENT(IN) :: delx, dely
         REAL(RP), DIMENSION(0:,0:), INTENT(OUT) :: U, V
         TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT) :: Partlines
         CHARACTER(LEN=30), INTENT(IN) :: problem
      END SUBROUTINE INIT_PARTICLES
   END INTERFACE

   INTERFACE
      SUBROUTINE SET_PART (Partline, x, y)
         USE nrtype
         IMPLICIT NONE
         REAL(RP), INTENT(IN) :: x, y
         TYPE(PARTICLELINE), INTENT(INOUT), TARGET :: Partline
      END SUBROUTINE SET_PART
   END INTERFACE

!   INTERFACE
!      SUBROUTINE MARK_CELLS (FLAG, imax, jmax, delx, dely, &
!                             ifull, isurf, N, Partlines)
!         USE nrtype
!         IMPLICIT NONE
!         INTEGER, INTENT(IN) :: imax, jmax
!         INTEGER, INTENT(OUT) :: ifull, isurf, N
!         REAL(RP), INTENT(IN) :: delx, dely
!         INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
!         TYPE(PARTICLELINE), DIMENSION(:), INTENT(INOUT), TARGET :: &
!                             Partlines
!      END SUBROUTINE MARK_CELLS
!   END INTERFACE

	INTERFACE
	SUBROUTINE MARK_CELLS (FLAG,FS, imax, jmax, delx, dely, ifull, isurf)

      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: ifull, isurf
      REAL(RP), INTENT(IN) :: delx, dely
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
	  REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: FS
	END SUBROUTINE MARK_CELLS
	END INTERFACE

	

   INTERFACE
      SUBROUTINE SET_UVP_SURFACE (U, V, P, FLAG, GX, GY, imax, jmax, &
                                  Re, delx, dely, delt)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imax, jmax
         REAL(RP), INTENT(IN) :: GX, GY, Re, delx, dely, delt
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P
      END SUBROUTINE SET_UVP_SURFACE
   END INTERFACE
!
! diagnostics
!
   INTERFACE
      SUBROUTINE OUTPUT_ONE_real(FIELD, imin, imax, jmin, jmax)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FIELD
      END SUBROUTINE OUTPUT_ONE_real
   END INTERFACE

   INTERFACE
      SUBROUTINE OUTPUT_ONE_integer(IFIELD, imin, imax, jmin, jmax, filename)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
         INTEGER(I2B), DIMENSION(0:,0:), INTENT(INOUT) :: IFIELD
         CHARACTER(LEN=30) filename
      END SUBROUTINE OUTPUT_ONE_integer
   END INTERFACE

   INTERFACE
      SUBROUTINE OUTPUT_THREE_real (U,V,P,imin,imax,jmin,jmax,filename)
         USE nrtype
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: imin, imax, jmin, jmax
         REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: U, V, P
         CHARACTER(LEN=30) filename
      END SUBROUTINE OUTPUT_THREE_real
   END INTERFACE


   INTERFACE
     SUBROUTINE VFCONV(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)

      use nrtype
      implicit none
      INTEGER, INTENT(IN) :: imax, jmax,cycles,wW, wE, wN, wS
      REAL(RP), INTENT(IN) :: delx, dely,delt
      REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
      REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
	END SUBROUTINE VFCONV
   END INTERFACE


   INTERFACE
     SUBROUTINE VTK(IOpst,filename1,imax, jmax,delx,dely,FS,U,V,P)
		use nrtype
		implicit none
		integer,intent(in):: IOpst
		character(200),intent(in)::filename1
		REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: FS,U,V,P
		REAL(RP), INTENT(IN) ::delx,dely
		INTEGER, INTENT(IN) :: imax, jmax
	END SUBROUTINE VTK
   END INTERFACE

   INTERFACE
     subroutine transport(c,itype, tanalfa, tanbeta, u1,u3,u4,u2,delx,dely,delt,f1,f3,f4,f2)	
		use nrtype
		implicit none	
		
		REAL(RP), INTENT(IN):: u1,u3,u4,u2,delx,dely,delt,c,tanalfa, tanbeta
		REAL(RP), INTENT(INOUT):: f1,f3,f4,f2
		INTEGER, INTENT(IN)::itype
      end subroutine transport
   END INTERFACE

   INTERFACE
     SUBROUTINE VFY(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)
	use nrtype
	implicit none

	INTEGER, INTENT(IN) :: imax, jmax,cycles,wW, wE, wN, wS
	REAL(RP), INTENT(IN) :: delx, dely,delt
	REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
	REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
	INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
	END SUBROUTINE VFY
   END INTERFACE

   
  INTERFACE 
   SUBROUTINE VFLAIR(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)
	use nrtype
	implicit none

	INTEGER, INTENT(IN) :: imax, jmax,cycles,wW, wE, wN, wS
	REAL(RP), INTENT(IN) :: delx, dely,delt
	REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
	REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
	INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
   END SUBROUTINE VFLAIR
 END INTERFACE


  INTERFACE
	SUBROUTINE VFCT(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)

    use nrtype
    implicit none
!INCLUDE 'interfaces.h'
!INCLUDE 'defs.h'
!global variables
   INTEGER, INTENT(IN) :: imax, jmax,cycles,wW, wE, wN, wS
   REAL(RP), INTENT(IN) :: delx, dely,delt
   REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
   REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
   INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG
   END SUBROUTINE VFCT


  END INTERFACE