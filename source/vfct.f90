!distance of point(x1,y1) and point(x2,y2) 
!      real function dist(x1,y1,x2,y2)  
!      implicit real(a-h,o-z) 
!      dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) 
!      end
	    
SUBROUTINE VFCT(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
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


integer I,J
REAL(RP) fl(0:imax+1),fh(0:imax+1),af(0:imax+1),q(0:imax+1)
REAL(RP) rplus(0:imax+1),rminus(0:imax+1)
REAL(RP) w1,w2,w3,w4,h,pplus,qplus,pminus,qminus,dv,cmax,cmin
REAL(RP) isize
REAL FSN(0:IMAX+1,0:JMAX+1)
REAL(RP), DIMENSION(:,:), ALLOCATABLE ::FS1,FS2
REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
ALLOCATE (FS1(0:IMAX+1,0:JMAX+1))
ALLOCATE (FS2(0:IMAX+1,0:JMAX+1))
h=delx
isize=imax

DO J=0,jmax
 DO I=0,imax
   FS1(I,J)=FS(I,J)   
 END DO
END DO
if(mod(cycles,2)==0)then
!X--direction 
        do i=0,imax 
          do j=0,jmax 
            FS(i,j)=FS(i,j)*delx*delx 
          enddo 
        enddo     
        do j=0,jmax 
          do i=0,imax-1 
            if(u(i,j).ge.0) then 
               fl(i)=u(i,j)*delt*FS(i,j) 
            else 
               fl(i)=u(i,j)*delt*FS(i+1,j) 
            endif   
           enddo     !!!!!  low-order numerical flux 
           do i=1,imax-1 
             FS1(i,j)=FS(i,j)-(fl(i)-fl(i-1))/delx 
           enddo 
           do i=0,imax-1 
             if(u(i,j).ge.0) then 
               fh(i)=u(i,j)*delt*FS(i+1,j) 
             else  
               fh(i)=u(i,j)*delt*FS(i,j) 
             endif 
           enddo     !!!!!  high-order numerical flux 
           do i=0,imax-1 
             af(i)=fh(i)-fl(i)    
           enddo     !!!!!  anti-disperse numerical flux                
           do i=1,imax-2 
             if(af(i)*(FS1(i+1,j)-FS1(i,j)).lt.0.and. &
                 (af(i)*(FS1(i+2,j)-FS1(i+1,j)).lt.0.or. &
                 af(i)*(FS1(i,j)-FS1(i-1,j)).lt.0)) then 
               af(i)=0 
             endif 
           enddo   
           do i=1,imax-1 
             w1=max(FS(i-1,j),FS1(i-1,j)) 
             w2=max(FS(i,j),FS1(i,j))  
             w3=max(FS(i+1,j),FS1(i+1,j)) 
             cmax=max(w1,w2,w3) 
             w1=min(FS(i-1,j),FS1(i-1,j)) 
             w2=min(FS(i,j),FS1(i,j))  
             w3=min(FS(i+1,j),FS1(i+1,j)) 
             cmin=min(w1,w2,w3) 
             pplus=max(0.0,af(i-1))-min(0.0,af(i)) 
             qplus=(cmax-FS1(i,j))*delx 
             if(pplus.gt.0) then 
               rplus(i)=min(1.0,qplus/pplus) 
             else 
               rplus(i)=0 
             endif     
             pminus=max(0.0,af(i))-min(0.0,af(i-1)) 
             qminus=(FS1(i,j)-cmin)*delx 
             if(pminus.gt.0.0) then 
               rminus(i)=min(1.0,qminus/pminus) 
             else 
               rminus(i)=0 
             endif        
           enddo 
           do i=1,imax-2                      
             if(af(i).ge.0) then 
               q(i)=min(rplus(i+1),rminus(i)) 
             else 
               q(i)=min(rplus(i),rminus(i+1))  
             endif 
           enddo   !!!! limit the anti-disperse numerical flux 
           do i=2,imax-2 
             FS2(i,j)=FS1(i,j)-(q(i)*af(i)-q(i-1)*af(i-1))/delx 
             dv=delx*delx-delt*delx*(u(i,j)-u(i-1,j)) 
             FS(i,j)=FS2(i,j)/dv 
           enddo                
        enddo     
!c       Y--direction 
        do i=0,imax 
          do j=0,jmax 
            FS(i,j)=FS(i,j)*dely*dely 
          enddo 
        enddo     
        do i=0,imax 
           do j=0,jmax-1 
            if(v(i,j).ge.0) then 
               fl(j)=v(i,j)*delt*FS(i,j) 
            else 
               fl(j)=v(i,j)*delt*FS(i,j+1) 
            endif   
           enddo 
           do j=1,jmax-1 
             FS1(i,j)=FS(i,j)-(fl(j)-fl(j-1))/dely 
           enddo    
           do j=1,jmax-2 
             if(af(j)*(FS1(i,j+1)-FS1(i,j)).lt.0.and. & 
                 (af(j)*(FS1(i,j+2)-FS1(i,j+1)).lt.0.or. &
                 af(j)*(FS1(i,j)-FS1(i,j-1)).lt.0)) then 
               af(j)=0 
             endif 
           enddo   
           do j=1,jmax-1 
             if(v(i,j).ge.0) then 
               fh(j)=v(i,j)*delt*FS(i,j+1) 
             else  
               fh(j)=v(i,j)*delt*FS(i,j) 
             endif 
           enddo 
           do j=0,jmax-1 
             af(j)=fh(j)-fl(j) 
           enddo          
           do j=1,jmax-1 
             w1=max(FS(i,j-1),FS1(i,j-1)) 
             w2=max(FS(i,j),FS1(i,j))  
             w3=max(FS(i,j+1),FS1(i,j+1)) 
             cmax=max(w1,w2,w3) 
             w1=min(FS(i,j-1),FS1(i,j-1)) 
             w2=min(FS(i,j),FS1(i,j))  
             w3=min(FS(i,j+1),FS1(i,j+1)) 
             cmin=min(w1,w2,w3) 
             pplus=max(0.0,af(j-1))-min(0.0,af(j)) 
             qplus=(cmax-FS1(i,j))*dely 
             if(pplus.gt.0) then 
               rplus(j)=min(1.0,qplus/pplus) 
             else 
               rplus(j)=0 
             endif     
             pminus=max(0.0,af(j))-min(0.0,af(j-1)) 
             qminus=(FS1(i,j)-cmin)*dely 
             if(pminus.gt.0.0) then 
               rminus(j)=min(1.0,qminus/pminus) 
             else 
               rminus(j)=0.0 
             endif        
           enddo 
           do j=1,jmax-2                      
             if(af(j).ge.0.0) then 
               q(j)=min(rplus(j+1),rminus(j)) 
             else 
               q(j)=min(rplus(j),rminus(j+1))  
             endif 
           enddo 
           do j=2,jmax-2 
             FS2(i,j)=FS1(i,j)-(q(j)*af(j)-q(j-1)*af(j-1))/dely 
             FS(i,j)=max(min(FS2(i,j)/dely/dely,1.0),0.0) 
           enddo                
        enddo  
        !t=t+delt 
       ! it=it+1  
!c     Y--direction !X and Y--direction solve alternately, 
!c     the goal is to avoid systemic error.
    else 

        do i=0,imax 
          do j=0,jmax 
            FS(i,j)=FS(i,j)*dely*dely 
          enddo 
       enddo     
        do i=0,imax 
           do j=0,jmax-1 
            if(v(i,j).ge.0) then 
               fl(j)=v(i,j)*delt*FS(i,j) 
            else 
               fl(j)=v(i,j)*delt*FS(i,j+1) 
            endif   
           enddo 
           do j=1,jmax-1 
             FS1(i,j)=FS(i,j)-(fl(j)-fl(j-1))/dely 
           enddo 
           do j=1,jmax-1 
             if(v(i,j).ge.0) then 
               fh(j)=v(i,j)*delt*FS(i,j+1) 
             else  
               fh(j)=v(i,j)*delt*FS(i,j) 
             endif 
           enddo 
           do j=0,jmax-1 
             af(j)=fh(j)-fl(j) 
           enddo  
           do j=1,jmax-2 
             if(af(j)*(FS1(i,j+1)-FS1(i,j)).lt.0.and. & 
                 (af(j)*(FS1(i,j+2)-FS1(i,j+1)).lt.0.or.& 
                 af(j)*(FS1(i,j)-FS1(i,j-1)).lt.0)) then 
               af(j)=0 
             endif 
           enddo         
           do j=1,jmax-1 
             w1=max(FS(i,j-1),FS1(i,j-1)) 
             w2=max(FS(i,j),FS1(i,j))  
             w3=max(FS(i,j+1),FS1(i,j+1)) 
             cmax=max(w1,w2,w3) 
             w1=min(FS(i,j-1),FS1(i,j-1)) 
             w2=min(FS(i,j),FS1(i,j))  
             w3=min(FS(i,j+1),FS1(i,j+1)) 
             cmin=min(w1,w2,w3) 
             pplus=max(0.0,af(j-1))-min(0.0,af(j)) 
             qplus=(cmax-FS1(i,j))*dely 
             if(pplus.gt.0) then 
               rplus(j)=min(1.0,qplus/pplus) 
             else 
               rplus(j)=0.0 
             endif     
             pminus=max(0.0,af(j))-min(0.0,af(j-1)) 
             qminus=(FS1(i,j)-cmin)*dely 
             if(pminus.gt.0) then 
               rminus(j)=min(1.0,qminus/pminus) 
             else 
               rminus(j)=0.0
             endif        
           enddo 
           do j=1,jmax-2                      
             if(af(j).ge.0.0) then 
               q(j)=min(rplus(j+1),rminus(j)) 
             else 
               q(j)=min(rplus(j),rminus(j+1))  
             endif 
           enddo 
           do j=2,jmax-2 
             FS2(i,j)=FS1(i,j)-(q(j)*af(j)-q(j-1)*af(j-1))/dely 
             dv=dely*dely-delt*dely*(v(i,j)-v(i,j-1)) 
             FS(i,j)=FS2(i,j)/dv 
           enddo                
        enddo  
!c        X--direction 
        do i=0,imax 
          do j=0,jmax 
            FS(i,j)=FS(i,j)*delx*delx 
          enddo 
        enddo     
        do j=0,jmax 
          do i=0,imax-1 
            if(u(i,j).ge.0) then 
               fl(i)=u(i,j)*delt*FS(i,j) 
            else 
               fl(i)=u(i,j)*delt*FS(i+1,j) 
            endif   
           enddo 
           do i=1,imax-1 
             FS1(i,j)=FS(i,j)-(fl(i)-fl(i-1))/delx 
           enddo   
           do i=1,imax-2 
             if(af(i)*(FS1(i+1,j)-FS1(i,j)).lt.0.and. &
                 (af(i)*(FS1(i+2,j)-FS1(i+1,j)).lt.0.or.& 
                 af(i)*(FS1(i,j)-FS1(i-1,j)).lt.0)) then 
               af(i)=0 
             endif 
           enddo 
           do i=0,imax-1 
             if(u(i,j).ge.0) then 
               fh(i)=u(i,j)*delt*FS(i+1,j) 
             else  
               fh(i)=u(i,j)*delt*FS(i,j) 
             endif 
           enddo 
           do i=0,imax-1 
             af(i)=fh(i)-fl(i) 
           enddo          
           do i=1,imax-1 
             w1=max(FS(i-1,j),FS1(i-1,j)) 
             w2=max(FS(i,j),FS1(i,j))  
             w3=max(FS(i+1,j),FS1(i+1,j)) 
             cmax=max(w1,w2,w3) 
             w1=min(FS(i-1,j),FS1(i-1,j)) 
             w2=min(FS(i,j),FS1(i,j))  
             w3=min(FS(i+1,j),FS1(i+1,j)) 
             cmin=min(w1,w2,w3) 
             pplus=max(0.0,af(i-1))-min(0.0,af(i)) 
             qplus=(cmax-FS1(i,j))*delx
             if(pplus.gt.0) then 
               rplus(i)=min(1.0,qplus/pplus) 
             else 
               rplus(i)=0 
             endif     
             pminus=max(0.0,af(i))-min(0.0,af(i-1)) 
             qminus=(FS1(i,j)-cmin)*delx
             if(pminus.gt.0) then 
               rminus(i)=min(1.0,qminus/pminus) 
             else 
               rminus(i)=0 
             endif        
           enddo 
           do i=1,imax-2                      
             if(af(i).ge.0) then 
               q(i)=min(rplus(i+1),rminus(i)) 
             else 
               q(i)=min(rplus(i),rminus(i+1))  
             endif 
           enddo 
           do i=2,imax-2 
             FS2(i,j)=FS1(i,j)-(q(i)*af(i)-q(i-1)*af(i-1))/delx 
             FS(i,j)=max(0.0,min(FS2(i,j)/delx/delx,1.0)) 
           enddo                
        enddo     
  end if       


DEALLOCATE(FS1,FS2)

END SUBROUTINE VFCT