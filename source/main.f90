!
!%%%%%%%%%%%%%%%%%%%%%%%%% PROGRAM MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
PROGRAM MAIN
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Fortran 90 re-write of C code to accompany text:
!
!     Numerical Simulation in Fluid Dynamics (SIAM 1998)
!       by Griebel, Dornseifer, and Neunhoeffer
!
!  Translated by:
!
!     Dr. C. David Pruett
!     Dept. of Mathematics & Statistics
!     MSC 7803
!     James Madison University
!     Harrisonburg, VA 22807 USA
!     dpruett@math.jmu.edu
!     www.math.jmu.edu/~dpruett
!     540-568-6227
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  USE nrtype
 
  IMPLICIT NONE

  INTEGER :: N, &
             imax, jmax, &
             wW, wE, wN, wS, &
             itermax, itersor, iwrite, p_bound 

  INTEGER :: ppc, &
             ifull, &
             isurf, &
             ibound, &
             init_case, &
             cycles, i

  INTEGER(I2B), DIMENSION(:,:), ALLOCATABLE :: FLAG

  REAL(RP) :: t, &
              xlength, ylength, &
              delx, dely, &
              t_end, delt, tau, &
              del_trace, del_inj, del_streak, del_vec, &
              pos1x, pos2x, pos1y, pos2y, &
              Re, Pr, &
              GX, GY, &
              UI, VI, TI, &
              beta, eps, omega, gamma, res


  REAL(RP), DIMENSION(:,:), ALLOCATABLE :: U, V, P, PSI, ZETA, &
                                        RHS, F, G, TEMP, HEAT, &
                                        WORK, WORK1,FS

  CHARACTER (LEN=30) :: problem, &
                        vecfile, tracefile, streakfile, &
                        inputfile, infile, outfile

  TYPE (PARTICLELINE), DIMENSION(:), ALLOCATABLE :: Partlines

  INTEGER :: status_var
  INTEGER(I4B) :: idim, jdim, kdim, kout

!
!--------------------------VTK OUTPUT-----------------------------------
!

character(200) filename1
character(50) szNumber
integer IOpst


INCLUDE 'interfaces.h'

IOpst=1


!
!------------------------- INITIALIZE ----------------------------------
!
  CALL READ_PARAMETERS (problem, &
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

  itersor = 0
  ifull = 0
  isurf = 0
  ibound = 0

  DX=DELX
  DY=DELY
!
! allocate arrays 
!
  ALLOCATE (U(0:imax+1,0:jmax+1))
  ALLOCATE (V(0:imax+1,0:jmax+1))
  ALLOCATE (P(0:imax+1,0:jmax+1))
  ALLOCATE (TEMP(0:imax+1,0:jmax+1))

  ALLOCATE (WORK(0:imax+1,0:jmax+1))
  ALLOCATE (WORK1(1:imax,1:jmax))

  ALLOCATE (F(0:imax+1,0:jmax+1))
  ALLOCATE (G(0:imax+1,0:jmax+1))
  ALLOCATE (RHS(0:imax+1,0:jmax+1))

  ALLOCATE (PSI(0:imax,0:jmax))
  ALLOCATE (HEAT(0:imax,0:jmax))
  ALLOCATE (ZETA(1:imax-1,1:jmax-1))

  ALLOCATE (FLAG(0:imax+1,0:jmax+1)) 

  ALLOCATE (FS(0:imax+1,0:jmax+1)) 
  ALLOCATE (DIV(0:imax+1,0:jmax+1))
!
! Read initial values from file "infile" if it exists
!
  CALL READ_BIN (infile, U, V, P, TEMP, FLAG, imax, jmax, init_case)
  IF ( init_case == 1 ) THEN
     CALL INIT_UVP (problem, U, V, P, TEMP, imax, jmax, UI, VI, TI)
     CALL INIT_FLAG (problem, FLAG, imax, jmax, delx, dely, ibound,FS)
  ELSE IF ( init_case == 2 ) THEN
     write (6,*) ' read error in READ_bin:  STOP'
     stop ' read_bin'
  ENDIF
!
! Initialize particles for streaklines or particle tracing
! (make an additional copy; one to hold injection points)
! 
!  IF ( ( streakfile /= "none" ) .OR. ( tracefile /= "none" ) ) THEN
!     ALLOCATE (Partlines(1:N))
!     DO i = 1,N
!        NULLIFY (Partlines(i)%Particles%next)
!     END DO
!     CALL SET_PARTICLES (N, pos1x, pos1y, pos2x, pos2y, Partlines)
!     CALL INJECT_PARTICLES (N, Partlines)  ! makes an extra copy
!  END IF
!
! Initialize particles for free-boundary problems
!
!  IF ( (problem == "drop"    ) .OR. (problem == "dam" ) ) THEN
!     DEALLOCATE(Partlines)
!     ppc = 4
!     IF (problem == "drop") THEN
!        N = 2
!     ELSE IF (problem == "dam") THEN
!        N = 1
!     END IF
!     ALLOCATE (Partlines(1:N))
!     DO i = 1,N
!        NULLIFY (Partlines(i)%Particles%next)
!     END DO
!     CALL INIT_PARTICLES (N, imax, jmax, delx, dely, ppc, problem, &
!                          U, V, Partlines)
!     CALL INJECT_PARTICLES (N, Partlines)  ! makes an dummy copy
!  END IF
!
! initialize boundary conditions
!
  CALL SETBCOND (U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS,FS)
  CALL SETSPECBCOND (problem, U, V, TEMP, imax, jmax, UI, VI)

  t = 0.0
  cycles = 0
!
!---------------------- TIME ADVANCEMENT LOOP --------------------------
!
  DO WHILE (t < t_end)
     
    CALL COMP_delt (delt, t, imax, jmax, delx, dely, &
                    U, V, Re, Pr, tau, iwrite, &
                    del_trace, del_inj, del_streak, del_vec)
    !
    ! determine fluid cells for free-boundary problems and set
    !   set boundary conditions at free surface
    !
    IF ( (problem == "drop"   ) .OR. (problem == "dam" ) .OR. & 
         (problem == "molding") .OR. (problem == "wave") ) THEN
       !CALL MARK_CELLS (FLAG, imax, jmax, delx, dely, ifull, isurf, &
                        !N, Partlines)
	CALL MARK_CELLS (FLAG,FS, imax, jmax, delx, dely, ifull, isurf)

       CALL SET_UVP_SURFACE (U, V, P, FLAG, GX, GY, imax, jmax, &
                             Re, delx, dely, delt)
    ELSE
       ifull = imax*jmax - ibound
    END IF
    !
    ! compute the new temperature
    ! 
    CALL COMP_TEMP (U, V, TEMP, FLAG, WORK, imax, jmax, delt, &
                    delx, dely, gamma, Re, Pr)
    !
    ! compute tentative velocity field (F,G)
    !
    CALL COMP_FG (U, V, TEMP, F, G, FLAG, imax, jmax, & 
                  delt, delx, dely, GX, GY, gamma, Re, beta)
    !
    ! compute right-hand side for pressure equation
    !
    CALL COMP_RHS (F, G, RHS, FLAG, imax, jmax, delt, delx, dely)
    !
    ! solve the pressure equation by successive over relaxation (SOR)
    !
    IF (ifull > 0) THEN
       CALL POISSON (P, RHS, FLAG, imax, jmax, delx, dely, &
                     eps, itersor, itermax, omega, res, ifull, p_bound)
       WRITE (6,'(1x,6ht_end=,g10.4,&
                  &1x,2ht=,g10.4,&
                  &1x,5hdelt=,g10.4,&
                  &1x,11hiterations=,i4,&
                  &1x,4hres=,g12.6,&
                  &1x,8hF-cells=,i4,&
                  &1x,8hS-cells=,i4,&
                  &1x,8hB-cells=,i4)') &
                  t_end, t+delt, delt, itersor, &
                  res, ifull, isurf, ibound
    END IF
    !
    ! compute the new velocity field
    !
    CALL ADAP_UV (U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely)
    !
    ! set the boundary conditions
    !
    CALL SETBCOND (U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS,FS)
    CALL SETSPECBCOND (problem, U, V, TEMP, imax, jmax, UI, VI)

    IF ( (problem == "drop"   ) .OR. (problem == "dam" ) .OR. & 
         (problem == "molding") .OR. (problem == "wave") ) THEN
       CALL SET_UVP_SURFACE (U, V, P, FLAG, GX, GY, imax, jmax, &
                             Re, delx, dely, delt)
!Для включения метода решения уравнения переноса закомментировать ненужный
!раскомментировать нужный. 
!VFCONV- стандартный метод Херта-Николса, VFY-метод Янга,
!VFLAIR- FLAIR метод, VFCT-flux-corrected transport 
!	   CALL	VFCONV(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
!                              FLAG,U, V,FS)

!	   CALL VFY(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
!                              FLAG,U, V,FS)

!	   CALL VFLAIR(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
!                              FLAG,U, V,FS)
		CALL VFCT(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)

    END IF
     
   ! IF ( (IAND(iwrite,8 )/=0) .AND. (vecfile  /= "none") ) THEN
   !    CALL COMP_PSI_ZETA (U, V, PSI, ZETA, FLAG, imax,jmax, delx,dely)
   !    CALL COMP_HEAT (U, V, TEMP, HEAT, FLAG, Re, Pr, &
   !                    imax, jmax, delx, dely)
   !    CALL OUTPUTVEC_bin (U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, &
    !                  WORK1, xlength, ylength, imax, jmax, vecfile)
   ! END IF

    !IF ( (IAND(iwrite,8)/=0) .AND. (outfile /= "none") ) THEN
    !   CALL WRITE_BIN (outfile, U, V, P, TEMP, FLAG, imax, jmax)
    !END IF

    !IF ( tracefile /= "none") THEN
    !    CALL PARTICLE_TRACING (tracefile, t, imax, jmax, &
    !                           delx, dely, delt, U, V, FLAG, &
    !                           N, Partlines, iwrite)
   ! END IF

    !IF ( streakfile /= "none") THEN
      ! CALL STREAKLINES (streakfile, iwrite, imax, jmax, &
     !                    delx, dely, delt, t, U, V, FLAG, &
     !                    N, Partlines)
   ! END IF
    !
    ! advance the time
    !
    t = t + delt
    cycles = cycles + 1

if(mod(cycles,10)==0) then

	Write(szNumber,'(i6.6)') cycles
   !filename1=Trim("../nast2dvof/vtk/")//Trim("uvpf")//Trim(szNumber)//Trim(".vtk")
   	filename1=Trim("D:/cfdprog/nast2dvof/vtk/")//Trim("uvpf")//Trim(szNumber)//Trim(".vtk")
    call VTK(IOpst,filename1,imax, jmax,delx,dely,FS,U,V,P)
end if	

  END DO

!  IF ( vecfile /= "none" ) THEN
!     CALL COMP_PSI_ZETA (U, V, PSI, ZETA, FLAG, imax,jmax, delx,dely)
!     CALL COMP_HEAT (U, V, TEMP, HEAT, FLAG, Re, Pr, &
!                     imax, jmax, delx, dely)
!     CALL OUTPUTVEC_bin (U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, &
!                    WORK1, xlength,ylength, imax,jmax, vecfile)
!  END IF

!  IF ( outfile /= "none" ) THEN
!     CALL WRITE_BIN (outfile, U, V, P, TEMP, FLAG, imax, jmax)
!  END IF
!
! FLOW VIZ
!     CALL OUTPUT_ONE_real (PSI, 0, imax, 0, jmax)
! FLOW VIZ
!
  !
  ! free memory
  !
  DEALLOCATE (U, V, P, TEMP, WORK, F, G, RHS, PSI, HEAT, ZETA, FLAG,FS,DIV)

  STOP ' normal'

END PROGRAM MAIN
