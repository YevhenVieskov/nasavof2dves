!	calculate the transportation from neighbour cell
	
	
subroutine transport(c,itype, tanalfa, tanbeta, u1,u3,u4,u2,delx,dely,delt,f1,f3,f4,f2)	
use nrtype
implicit none
INCLUDE 'interfaces.h'
INCLUDE 'defs.h'
!global variables
REAL(RP), INTENT(IN):: u1,u3,u4,u2,delx,dely,delt,c,tanalfa, tanbeta
REAL(RP), INTENT(INOUT):: f1,f3,f4,f2
INTEGER, INTENT(IN)::itype
!local variables
real cotalfa, cotbeta,s1,s2,s3,s4,temp1 
 cotalfa=1.0/tanalfa 
 cotbeta=1.0/tanbeta 
       
      f1=0 
      f3=0 
      f4=0 
      f2=0       
      if(itype.eq.1) then 
        s1=0 
        s3=sqrt(2.0*c*cotalfa) 
        s4=0 
        s2=sqrt(2.0*c*tanalfa) 
        if(u1.gt.0) then 
          if(u1*delt.le.(1.0-s2)*dely) then 
            f1=0 
          else                           
            temp1=u1*delt-(1.0-s2)*dely 
            f1=0.5*temp1*temp1*cotbeta 
          endif 
        endif   
        if(u2.gt.0) then 
          if(u2*delt.ge.s3*delx) then 
            f2=c*delx*dely 
          else 
            f2=0.5*u2*delt*(2.0-u2*delt/(s3*delx))*s2*dely 
          endif 
        endif   
        if(u3.lt.0) then 
          if(abs(u3)*delt.ge.s2*dely) then 
            f3=c*delx*dely 
          else 
            f3=0.5*abs(u3)*delt*(2.0-abs(u3)*delt/(s2*dely))*s3*delx 
          endif 
        endif   
        if(u4.lt.0) then 
          if(abs(u4)*delt.le.(1.0-s3)*delx) then 
            f4=0 
          else    
            temp1=abs(u4)*delt-(1.0-s3)*delx 
            f4=0.5*temp1*temp1*tanbeta 
          endif 
        endif 
      else if(itype.eq.2) then 
        s1=0 
        s3=1.0 
        s4=c-0.5*tanalfa 
        s2=c+0.5*tanalfa 
        if(u1.gt.0) then 
          if(u1*delt.le.(1.0-s2)*dely) then       
            f1=0 
          else if(u1*delt.le.(1.0-s4)*dely) then 
            temp1=u1*delt-(1.0-s2)*dely 
            f1=0.5*temp1*temp1*cotbeta 
          else 
            f1=u1*delt*delx-(1.0-c)*delx*dely 
          endif 
        endif   
        if(u2.gt.0) then 
          f2=u2*delt*(s2*dely-0.5*u2*delt*tanbeta) 
        endif   
        if(u3.lt.0) then 
          if(abs(u3)*delt.le.s4*dely) then 
            f3=abs(u3)*delt*delx 
          else if(abs(u3)*delt.le.s2*dely) then  
            temp1=abs(u3)*delt-s4*dely 
            f3=abs(u3)*delt*delx-0.5*temp1*temp1*cotbeta 
          else 
            f3=c*delx*dely 
          endif 
        endif 
        if(u4.lt.0) then 
          f4=abs(u4)*delt*(s4*dely+0.5*abs(u4)*delt*tanbeta) 
        endif 
      else if(itype.eq.3) then 
        s1=c-0.5*cotalfa 
        s3=c+0.5*cotalfa 
        s4=0 
        s2=1.0 
        if(u1.gt.0) then 
          f1=u1*delt*(s1*delx+0.5*u1*delt*cotbeta) 
        endif 
        if(u2.gt.0) then 
          if(u2*delt.le.s1*delx) then 
            f2=u2*delt*dely 
          else if(u2*delt.le.s3*delx) then 
            temp1=u2*delt-s1*delx 
            f2=u2*delt*dely-0.5*temp1*temp1*tanbeta 
          else 
            f2=c*delx*dely 
          endif 
        endif 
        if(u3.lt.0) then 
          f3=abs(u3)*delt*(s3*delx-0.5*abs(u3)*delt*cotbeta) 
        endif 
        if(u4.lt.0) then 
          if(abs(u4)*delt.le.(1.0-s3)*delx) then 
            f4=0 
          else if(abs(u4)*delt.le.(1.0-s1)*delx) then 
            temp1=abs(u4)*delt-(1.0-s3)*delx 
            f4=0.5*temp1*temp1*tanbeta 
          else 
            f4=abs(u4)*delt*dely-(1.0-c)*delx*dely 
          endif 
        endif 
      else if(itype.eq.4) then 
        s1=1.0-sqrt(2.0*(1.0-c)*cotalfa) 
        s3=1.0 
        s4=1.0-sqrt(2.0*(1.0-c)*tanalfa) 
        s2=1.0 
        if(u1.gt.0) then 
          if(u1*delt.ge.(1.0-s4)*dely) then 
            f1=u1*delt*delx-(1.0-c)*delx*dely 
          else   
            f1=u1*delt*(s1*delx+0.5*u1*delt*cotbeta) 
          endif 
        endif 
        if(u2.gt.0) then 
          if(u2*delt.le.s1*delx) then 
            f2=u2*delt*dely 
          else                        
            temp1=u2*delt-s1*delx 
            f2=u2*delt*dely-0.5*tanbeta*temp1*temp1 
          endif 
        endif 
        if(u3.lt.0) then 
          if(abs(u3)*delt.le.s4*dely) then 
            f3=abs(u3)*delt*delx 
          else                                    
            temp1=abs(u3)*delt-s4*dely 
            f3=abs(u3)*delt*delx-0.5*temp1*temp1*cotbeta 
          endif 
        endif 
        if(u4.lt.0) then 
          if(abs(u4)*delt.ge.(1.0-s1)*delx) then 
            f4=abs(u4)*delt*dely-(1.0-c)*delx*dely 
          else 
            f4=abs(u4)*delt*(s4*dely+0.5*abs(u4)*delt*tanbeta) 
          endif 
        endif 
      else 
        write(*,*)'error in subroutine transport' 
        stop 
      endif 


	
end subroutine transport	
	

SUBROUTINE VFY(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
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
REAL(RP) emk2,c,ft,fl,fr,fb,ut,ul,ur,ub,rnx,rny
REAL(RP) u1,u2,u3,u4,f1,f2,f3,f4,rnx1,rny1
REAL(RP) tanbeta,cotbeta,tanalfa,cotalfa
integer I,J,itype
real RDX,RDY
REAL FSN(0:IMAX+1,0:JMAX+1)
REAL(RP), DIMENSION(:,:), ALLOCATABLE ::FS1,FS0
REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
ALLOCATE (FS1(0:IMAX+1,0:JMAX+1))	
ALLOCATE (FS0(0:IMAX+1,0:JMAX+1))

DO J=0,jmax
 DO I=0,imax
   FS0(I,J)=0.0
   FS1(I,J)=0.0
 END DO
END DO

DO J=0,jmax
 DO I=0,imax
   FS1(I,J)=FS(I,J)   
 END DO
END DO

emk2=0

	   do i=1,imax
		 do j=1,jmax
		  IF(IOR(FLAG(I,J),C_B)/=0) THEN
			ft=0 
            fb=0 
            fl=0 
            fr=0  
            c=FS(i,j)     
            ut=v(i,j) 
            ub=v(i,j-1) 
            ul=u(i-1,j) 
            ur=u(i,j) 
            if(c.ge.1.0-emf) then 
              if(ut.gt.0) then 
                ft=ut*delt*delx*c 
              endif   
              if(ub.lt.0) then 
                fb=abs(ub)*delt*delx*c 
              endif 
              if(ul.lt.0) then 
                fl=abs(ul)*delt*dely*c 
              endif 
              if(ur.gt.0) then 
                fr=ur*delt*dely*c 
              endif   
            else if(c.gt.0) then 
              rnx=(FS(i+1,j+1)+2.0*FS(i+1,j)+FS(i+1,j-1) &
                 -FS(i-1,j+1)-2.0*FS(i-1,j)-FS(i-1,j-1))/delx 
              rny=(FS(i+1,j+1)+2.0*FS(i,j+1)+FS(i-1,j+1) &
                 -FS(i+1,j-1)-2.0*FS(i,j-1)-FS(i-1,j-1))/dely  
              if(abs(rny).le.emk2) then 
                if(rnx.ge.0) then 
                  u1=ut 
                  u2=ur 
                  u3=ub 
                  u4=ul 
                else 
                  u1=ut 
                  u2=-ul 
                  u3=ub 
                  u4=-ur 
                endif  
                f1=0 
                f2=0 
                f3=0 
                f4=0    
                if(u1.gt.0) then 
                  f1=abs(u1)*delt*c*delx 
                endif 
                if(u3.lt.0) then 
                  f3=abs(u3)*delt*c*delx 
                endif 
                if(u2.gt.0) then 
                  if(abs(u2)*delt.le.c*delx) then 
                    f2=abs(u2)*delt*dely 
                  else 
                    f2=c*delx*dely 
                  endif 
                endif 
                if(u4.lt.0) then 
                  if(abs(u4)*delt.le.(1.0-c)*delx) then 
                    f4=0 
                  else 
                    f4=(abs(u4)*delt-(1.0-c)*delx)*dely 
                  endif 
                endif 
                if(rnx.ge.0) then 
                  ft=f1 
                  fb=f3 
                  fr=f2 
                  fl=f4 
                else 
                  ft=f1 
                  fb=f3 
                  fr=f4 
                  fl=f2 
                endif  
              else if(abs(rnx).le.emk2) then 
                if(rny.ge.0) then 
                  u1=ut 
                  u2=ur 
                  u3=ub 
                  u4=ul 
                else 
                  u1=-ub 
                  u2=ur 
                  u3=-ut 
                  u4=ul 
                endif 
                f1=0 
                f2=0 
                f3=0 
                f4=0     
                if(u1.gt.0) then 
                  if(u1*delt.le.c*dely) then 
                    f1=u1*delt*delx 
                  else 
                    f1=c*delx*dely 
                  endif 
                endif 
                if(u3.lt.0) then 
                  if(abs(u3)*delt.le.(1.0-c)*dely) then 
                    f3=0 
                  else 
                    f3=(abs(u3)*delt-(1.0-c)*dely)*delx 
                  endif   
                endif 
                if(u2.gt.0) then 
                  f2=u2*delt*c*dely 
                endif 
                if(u4.lt.0) then 
                  f4=abs(u4)*delt*c*dely 
                endif   
                if(rny.ge.0) then 
                  ft=f1 
                  fb=f3 
                  fr=f2 
                  fl=f4 
                else 
                  ft=f3 
                  fb=f1 
                  fr=f2 
                  fl=f4 
                endif 
              else   
                if(rnx.gt.0.and.rny.lt.0) then 
                  u1=ut 
                  u2=ur 
                  u3=ub 
                  u4=ul 
                else if(rnx.lt.0.and.rny.gt.0) then 
                  u1=-ub 
                  u2=-ul 
                  u3=-ut 
                  u4=-ur 
                else if(rnx.lt.0.and.rny.lt.0) then 
                  u1=ut 
                  u2=-ul 
                  u3=ub 
                  u4=-ur                     
                else 
                  u1=-ub 
                  u2=ur 
                  u3=-ut 
                  u4=ul 
                endif 
                rnx1=abs(rnx) 
                rny1=-abs(rny)        
                tanbeta=-rnx1/rny1 
                cotbeta=1.0/tanbeta 
                tanalfa=delx/dely*tanbeta 
                cotalfa=1.0/tanalfa 
                if(tanalfa.le.1.0) then 
                  if(c.le.0.5*tanalfa) then 
                    itype=1 
                  else if(c.le.1.0-0.5*tanalfa) then 
                    itype=2 
                  else 
                    itype=4 
                  endif 
                else 
                  if(c.le.0.5*cotalfa) then 
                    itype=1 
                  else if(c.le.1.0-0.5*cotalfa) then 
                    itype=3 
                  else 
                    itype=4 
                  endif 
                endif                           
                call transport(c,itype,tanalfa,tanbeta,u1,u3,u4,u2, &
                          delx,dely,delt,f1,f3,f4,f2) 
                if(rnx.gt.0.and.rny.lt.0) then 
                  ft=f1 
                  fb=f3 
                  fl=f4 
                  fr=f2 
                else if(rnx.lt.0.and.rny.gt.0) then 
                  ft=f3 
                  fb=f1 
                  fl=f2 
                  fr=f4            
                else if(rnx.lt.0.and.rny.lt.0) then 
                  ft=f1 
                  fb=f3 
                  fl=f2 
                  fr=f4            
                else 
                  ft=f3 
                  fb=f1 
                  fl=f4 
                  fr=f2            
                endif 
              endif 
            endif    
            FS1(i,j+1)=FS1(i,j+1)+ft/delx/dely 
            FS1(i,j-1)=FS1(i,j-1)+fb/delx/dely 
            FS1(i+1,j)=FS1(i+1,j)+fr/delx/dely 
            FS1(i-1,j)=FS1(i-1,j)+fl/delx/dely 
            FS1(i,j)=FS1(i,j)-(ft+fb+fl+fr)/delx/dely 
				end if !C_B
				enddo
			enddo
			do i=1,imax
				do j=1,jmax
					FS(i,j)=max(min(FS1(i,j),1.0),0.0)
				enddo
			enddo


			
	
	DEALLOCATE (FS0,FS1)	
END SUBROUTINE VFY	
	
	
	




