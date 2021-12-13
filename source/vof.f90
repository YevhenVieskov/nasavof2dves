 !SUBROUTINE MARK_CELLS (FLAG,FS, imax, jmax, delx, dely, ifull, isurf)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
! Mark the cells of the fluid domain                            
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
!      USE nrtype
!      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: imax, jmax
!      INTEGER, INTENT(OUT) :: ifull, isurf
!      REAL(RP), INTENT(IN) :: delx, dely
!      INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
!	  REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: FS
      !TYPE(PARTICLELINE),DIMENSION(:),INTENT(INOUT),TARGET :: Partlines
!
!  local variables
!
!      INCLUDE 'defs.h'
!      INTEGER :: i, j, k
!      REAL(RP) :: x, y
      !TYPE(PARTICLE), POINTER :: part, temp, help
!	  REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
!
!  Set all cells which are not obstacle cells to empty cells 
!
!      DO i = 0,imax+1
!         DO j = 0,jmax+1
!            IF (FLAG(i,j) >= C_F.and.FS(i,j)<=EMF) THEN
!               FLAG(i,j) = IAND(IOR(FLAG(i,j),C_E),NOT(C_NSWO))
!            END IF
!         END DO
!      END DO
!
!  Mark cells containing LIQUID as fluid cells (loop over particles) 
!


!DO j = 1,jmax
! DO i = 1,imax
	
!    IF(FS(I,J)>EMF1) FLAG(I,J)=C_F
! END DO
!END DO



!
! Mark surface cells 
!
!      ifull = 0
!      isurf = 0

!      DO j = 1,jmax
!         DO i = 1,imax

!            IF ( (IAND(FLAG(i,j),C_F) /= 0) .AND. (FLAG(i,j) < C_E) ) THEN
         
!               IF (IAND(FLAG(i-1,j),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_W)
!               END IF

!               IF (IAND(FLAG(i+1,j),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_O)
!               END IF

!               IF (IAND(FLAG(i,j-1),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_S)
!               END IF

!               IF (IAND(FLAG(i,j+1),C_E) /= 0) THEN
!                  FLAG(i,j) = IOR(FLAG(i,j),C_N)
!               END IF

!               IF (FLAG(i,j) < C_O) THEN  
!                  ifull = ifull + 1
!               ELSE
!                  isurf = isurf + 1
!               END IF

!            END IF
!      
!         END DO
!      END DO
!
!  DIAGNOSTIC:  Output geometry of the fluid domain
!
!     WRITE (6,*) ' '
!     WRITE (6,*) ' Geometry of the fluid domain'
!     WRITE (6,*) ' '
!     DO j = jmax+1,0,-1
!        WRITE (6,'(1x,200i5)') (FLAG(i,j), i=0,imax+1)
!     END DO

!      RETURN

!   END SUBROUTINE MARK_CELLS





SUBROUTINE MARK_CELLS (FLAG,FS, imax, jmax, delx, dely, ifull, isurf)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
! Mark the cells of the fluid domain                            
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
      USE nrtype
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imax, jmax
      INTEGER, INTENT(OUT) :: ifull, isurf
      REAL(RP), INTENT(IN) :: delx, dely
      INTEGER(I2B), DIMENSION(0:,0:), INTENT(OUT) :: FLAG
	  REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: FS
      !TYPE(PARTICLELINE),DIMENSION(:),INTENT(INOUT),TARGET :: Partlines
!
!  local variables
!
      INCLUDE 'defs.h'
      INTEGER :: i, j, k
      REAL(RP) :: x, y
      !TYPE(PARTICLE), POINTER :: part, temp, help
	  REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
!
!  Set all cells which are not obstacle cells to empty cells 
!
      DO i = 0,imax+1
         DO j = 0,jmax+1
            IF (FLAG(i,j) >= C_F.and.FS(i,j)>EMF) THEN
               FLAG(i,j) = IAND(IOR(FLAG(i,j),C_E),NOT(C_NSWO))
            END IF
         END DO
      END DO
!
!  Mark cells containing LIQUID as fluid cells (loop over particles) 
!


DO j = 1,jmax
 DO i = 1,imax
	
    IF(FS(I,J)>EMF) FLAG(I,J)=C_F
 END DO
END DO



!
! Mark surface cells 
!
      ifull = 0
      isurf = 0

      DO j = 1,jmax
         DO i = 1,imax

            IF ( (IAND(FLAG(i,j),C_F) /= 0) .AND. (FLAG(i,j) < C_E) ) THEN
         
               IF (IAND(FLAG(i-1,j),C_E) /= 0) THEN
                  FLAG(i,j) = IOR(FLAG(i,j),C_W)
               END IF

               IF (IAND(FLAG(i+1,j),C_E) /= 0) THEN
                  FLAG(i,j) = IOR(FLAG(i,j),C_O)
               END IF

               IF (IAND(FLAG(i,j-1),C_E) /= 0) THEN
                  FLAG(i,j) = IOR(FLAG(i,j),C_S)
               END IF

               IF (IAND(FLAG(i,j+1),C_E) /= 0) THEN
                  FLAG(i,j) = IOR(FLAG(i,j),C_N)
               END IF

               IF (FLAG(i,j) < C_O) THEN  
                  ifull = ifull + 1
               ELSE
                  isurf = isurf + 1
               END IF

            END IF
      
         END DO
      END DO
!
!  DIAGNOSTIC:  Output geometry of the fluid domain
!
!     WRITE (6,*) ' '
!     WRITE (6,*) ' Geometry of the fluid domain'
!     WRITE (6,*) ' '
!     DO j = jmax+1,0,-1
!        WRITE (6,'(1x,200i5)') (FLAG(i,j), i=0,imax+1)
!     END DO

      RETURN

   END SUBROUTINE MARK_CELLS








SUBROUTINE VFCONV(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
                              FLAG,U, V,FS)
!
!moves fluid i.e. adjusts VOF function, and re-labels

use nrtype
implicit none
INCLUDE 'interfaces.h'
INCLUDE 'defs.h'
!global variables
INTEGER, INTENT(IN) :: imax, jmax,cycles,wW, wE, wN, wS
REAL(RP), INTENT(IN) :: delx, dely,delt
REAL(RP), DIMENSION(0:,0:), INTENT(IN) :: U, V
REAL(RP), DIMENSION(0:,0:), INTENT(INOUT) :: FS
INTEGER(I2B), DIMENSION(0:,0:), INTENT(IN) :: FLAG

!local variables
real DXR,DXL,DYT,DYB,RXDEN,RYDEN,FL,FC,FR,FB,FT
real AVFT,AVFB,AVFCY,AVFCX,AVFL,AVFR,PFX,PFY,ABPFX,ABPFY
real VX,VY,RB,RA,RD,FDM,FX1,FX,FY,FY1,ABVX,ABVY
integer I,J,IA,ID,IAD,IDM,JAD,JA,JD,JDM,MK
real RDX,RDY
REAL FSN(0:IMAX+1,0:JMAX+1)
INTEGER, DIMENSION(:,:), ALLOCATABLE ::IHV
REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
ALLOCATE (IHV(0:IMAX+1,0:JMAX+1))
RDX=1./DELX
RDY=1./DELY

DO J=0,jmax
 DO I=0,imax
   FSN(I,J)=FS(I,J)
 END DO
END DO

DO J=0,jmax+1
 DO I=0,imax+1
   ihv(I,J)=0
 END DO
END DO

!IF(cycles.ge.1) THEN
DO J=1,jmax
 DO I=1,imax
   IF(.NOT.(FS(I,J).LT.EMF.OR.FS(I,J).GT.EMF1)) THEN

    DXR=0.5*(DELX+DELX)
    DXL=0.5*(DELX+DELX)
    DYT=0.5*(DELY+DELY)
	DYB=0.5*(DELY+DELY)
	RXDEN=1.0/(DXR*DXL*(DXR+DXL))
	RYDEN=1.0/(DYT*DYB*(DYT+DYB))
	FL=FS(I-1,J+1)
	IF(IAND(FLAG(I-1,J+1),C_B)/=0.OR.(I.EQ.1.AND.WW.LT.3)) FL=1.0	
	FC=FS(I,J+1)	
	IF(IAND(FLAG(I,J+1),C_B)/=0) FC=1.0
	FR=FS(I+1,J+1)
	IF(IAND(FLAG(I+1,J+1),C_B)/=0.OR.(I.EQ.IMAX.AND.WE.LT.3)) FR=1.0
	AVFT=FL*DELX+FC*DELX+FR*DELX
	FL=FS(I-1,J-1)
	
	IF(IAND(FLAG(I-1,J-1),C_B)/=0.OR.(I.EQ.1.AND.WW.LT.3)) FL=1.0
	FC=FS(I,J-1)	
	IF(IAND(FLAG(I,J-1),C_B)/=0) FC=1.0
	FR=FS(I+1,J-1)	
	IF(IAND(FLAG(I+1,J-1),C_B)/=0.OR.(I.EQ.IMAX.AND.WE.LT.3)) FR=1.0
	AVFB=FL*DELX+FC*DELX+FR*DELX
	FL=FS(I-1,J)	
	IF(IAND(FLAG(I-1,J-1),C_B)/=0.OR.(I.EQ.1.AND.WW.LT.3)) FL=1.0
	FR=FS(I+1,J)
	
	IF(IAND(FLAG(I+1,J),C_B)/=0.OR.(I.EQ.IMAX.AND.WE.LT.3)) FR=1.0
	AVFCY=FL*DELX+FS(I,J)*DELX+FR*DELX
	FB=FS(I,J-1)	
	IF(IAND(FLAG(I+1,J),C_B)/=0.OR.(J.EQ.1.AND.WS.LT.3)) FB=1.0
	FT=FS(I,J+1)
	IF(IAND(FLAG(I,J+1),C_B)/=0.OR.(J.EQ.JMAX.AND.WN.LT.3)) FT=1.0
	AVFCX=FB*DELY+FS(I,J)*DELY+FT*DELY
	FB=FS(I-1,J-1)	
	IF(IAND(FLAG(I-1,J-1),C_B)/=0.OR.(J.EQ.1.AND.WS.LT.3)) FB=1.0
	FC=FS(I-1,J)
	IF(IAND(FLAG(I-1,J),C_B)/=0) FC=1.0
	FT=FS(I-1,J+1)	
	IF(IAND(FLAG(I-1,J+1),C_B)/=0.OR.(J.EQ.JMAX.AND.WN.LT.3)) FT=1.0
	AVFL=FB*DELY+FC*DELY+FT*DELY
	FB=FS(I+1,J-1)
	IF(IAND(FLAG(I+1,J-1),C_B)/=0.OR.(J.EQ.2.AND.WS.LT.3)) FB=1.0
	FC=FS(I+1,J)	
	IF(IAND(FLAG(I+1,J),C_B)/=0) FC=1.0
	FT=FS(I+1,J+1)
	IF(IAND(FLAG(I+1,J+1),C_B)/=0.OR.(J.EQ.JMAX.AND.WN.LT.3)) FT=1.0
	AVFR=FB*DELY+FC*DELY+FT*DELY

	PFX=RXDEN*((AVFR-AVFCX)*DXL**2+(AVFCX-AVFL)*DXR**2)
	PFY=RYDEN*((AVFT-AVFCY)*DYB**2+(AVFCY-AVFB)*DYT**2)

	ABPFX=ABS(PFX)
	ABPFY=ABS(PFY)
	IF(ABPFX.LE.ABPFY) THEN
	 ihv(i,j)=0
	ELSE
	 ihv(i,j)=1
	END IF
   END IF   
 END DO
END DO



DO j=1,JMAX
 DO i=1,IMAX
 !X-DIRECTION
 IF(IOR(FLAG(I,J),C_B)/=0) THEN ! IF (((IAND(FLAG(i,j),  C_F)/=0).AND.(FLAG(i,j)<C_E))
   VX=U(i,j)*DelT
   ABVX=ABS(VX)
   IF(VX.GE.0) THEN
	 IA=I+1
	 ID=I
	 IDM=MAX0(I-1,1)
!	 RB=X(I)
!	 RA=XI(I+1)
!	 RD=XI(I)
   ELSE
     IA=I
	 ID=I+1
	 IDM=MIN0(I+2,IMAX)
!	 RA=XI(I)
!     RD=XI(I+1)
   END IF

   IF(IHV(ID,J).EQ.1) THEN
	 IAD=IA
   ELSE
     IAD=ID
   END IF
   IF(FSN(IA,J).LT.EMF.OR.FSN(IDM,J).LT.EMF) IAD=IA
      FDM=MAX(FSN(IDM,J),FSN(ID,J))
      FX1=FSN(IAD,J)*ABS(VX)+MAX((FDM-FSN(IAD,J))*ABS(VX)-(FDM-FSN(ID,J))&
      *DELX,0.0)
      FX=MIN(FX1,FSN(ID,J)*DELX)
      FS(ID,J)=FS(ID,J)-FX*RDX
      FS(IA,J)=FS(IA,J)+FX*RDX
	  !Y-DIRECTION
	  VY=V(i,j)*DelT
      ABVY=ABS(VY)
	  IF(VY.GE.0.0) THEN
	    JA=J+1
		JD=J
		JDM=MAX0(J-1,1)
	  ELSE
	    JA=J
		JD=J+1
		JDM=MIN0(J+2,JMAX)
	  END IF

	  IF(IHV(I,JD).EQ.1) THEN
	   JAD=JD
      ELSE
       JAD=JA
      END IF

	  IF(FSN(I,JA).LT.EMF.OR.FSN(I,JDM).LT.EMF) JAD=JA
      FDM=MAX(FSN(I,JDM),FSN(I,JD))
      FY1=FSN(I,JAD)*ABS(VY)+MAX((FDM-FSN(I,JAD))*ABS(VY)-(FDM-FSN(I,JD))&
      *DELY,0.0)
      FY=MIN(FY1,FSN(I,JD)*DELY)
      FS(I,JD)=FS(I,JD)-FY*RDY
      FS(I,JA)=FS(I,JA)+FY*RDY
 END IF
 END DO
END DO

!END IF

!DO  J=1,JMAX
!  DO I=1,IMAX
!  IF(FS(I,J).GE.EMF1)THEN
!  IF(FS(I+1,J).LT.EMF.OR.FS(I-1,J).LT.EMF.OR.FS(I,J+1).LT.EMF.OR.FS(I,J-1).LT.EMF) THEN     
!   FS(I,J)=FS(I,J)-1.1*EMF
!  END IF
!  END IF
!  END DO
!END DO

DO J=0,jmax
 DO I=0,imax
   FSN(I,J)=FS(I,J)
 END DO
END DO

do i=1,imax 
  do j=1,jmax   
    fs(i,j)=min(1.0,max(fsn(i,j),0.0)) 
   enddo 
enddo    

DEALLOCATE (IHV)

END SUBROUTINE VFCONV