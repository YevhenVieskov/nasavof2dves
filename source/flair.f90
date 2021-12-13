real function dist(x1,y1,x2,y2)   
      implicit real (a-h,o-z) 
      dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) 
      end 
!c    approximate the slope a of interface according  
!c    to two neighbour cell's fluid volume function 
!c    fa and fb 
      subroutine calslope(fa,fb,a) 
      implicit real (a-h,o-z) 
         
      temp1=sqrt((1.0-fa)*(2.0-fa-fb)) 
      if(3.0*fb.ge.fa.and.3.0*fa-fb.le.2.0) then 
        a=fb-fa 
      else if(3.0*fb.le.fa.and.fa+fb-sqrt(fa*fb+fb*fb).le.0.5) then 
        bstar=2.0*(fa+fb-sqrt(fa*fb+fb*fb)) 
        a=2.0*(fa-bstar)                    
      else if(fa+fb+temp1.ge.1.5.and.fa+fb/3.0+temp1.le.4.0/3.0) then 
        a=2.0*(2.0*fa+fb-3.0+2.0*temp1) 
      else                         
        tl=1.0/(sqrt(fb/(1.0-fa))+1.0) 
        xastar=1.0-2.0*(1.0-fa)/tl 
        xbstar=1.0+2.0*fb/(1.0-tl) 
        a=1.0/(xastar-xbstar) 
      endif 
      end                              
!c    calculate the transportation from one cell to neighbour cell  
!c    for class I 
      subroutine calflux(fa,fb,u,h,delt,flux) 
      implicit real (a-h,o-z) 
      s=abs(u*delt/h) 
      temp1=sqrt((1.0-fa)*(2.0-fa-fb)) 
      if(3.0*fb.ge.fa.and.3.0*fa-fb.le.2.0) then 
        a=fb-fa 
        bstar=0.5*(3.0*fa-fb) 
        b=bstar*h 
        if(u.ge.0) then 
          flux=s*(a+bstar-a*s/2.0) 
        else 
          flux=s*(a+bstar+a*s/2.0) 
        endif  
      else if(3.0*fb.le.fa.and.fa+fb-sqrt(fa*fb+fb*fb).le.0.5) then 
        bstar=2.0*(fa+fb-sqrt(fa*fb+fb*fb)) 
        a=2.0*(fa-bstar)                    
        b=bstar*h   
        xb=-b/a 
        xbstar=xb/h 
        if(u.ge.0) then 
          flux=s*(a+bstar-a*s/2.0) 
        else 
          if(s+1.0.ge.xbstar) then 
            flux=fb 
          else 
            flux=s*(a+bstar+a*s/2.0) 
          endif 
        endif                 
      else if(fa+fb+temp1.ge.1.5.and.fa+fb/3.0+temp1.le.4.0/3.0) then 
        a=2.0*(2.0*fa+fb-3.0+2.0*temp1) 
        bstar=3.0*(3.0-2.0*fa-2.0/3.0*fb-2.0*temp1) 
        b=bstar*h 
        xa=(h-b)/a 
        xastar=xa/h 
        if(u.ge.0) then 
          if(s.le.1.0-xastar) then 
            flux=s*(a+bstar-a*s/2.0) 
          else 
            flux=s+fa-1.0 
          endif 
        else 
          flux=s*(a+bstar+a*s/2.0) 
        endif 
      else                         
        tl=1.0/(sqrt(fb/(1.0-fa))+1.0) 
        xastar=1.0-2.0*(1.0-fa)/tl 
        xbstar=1.0+2.0*fb/(1.0-tl) 
        a=1.0/(xastar-xbstar) 
        bstar=xbstar/(xbstar-xastar) 
        b=bstar*h 
        if(u.ge.0) then 
          if(s.le.1.0-xastar) then 
            flux=s*(a+bstar-a*s/2.0)   
          else 
            flux=s+fa-1.0      
          endif 
        else 
          if(s+1.0.ge.xbstar) then 
            flux=fb 
          else 
            flux=s*(a+bstar+a*s/2.0) 
          endif 
        endif 
      endif 
      end               
!c    calculate the transportation from one cell to neighbour cell  
!c    for class II       
      subroutine cal2flux(f,h,u,delt,beta,flux1,flux2) 
      implicit real(a-h,o-z) 
       
      if(beta.lt.0) stop 'In subroutine cal2flux' 
      a=-1.0/beta 
      s=abs(u*delt/h) 
      if(f.ge.1.0-1.0/(2.0*beta).and.f.ge.1.0-beta/2.0) then  
        b1=1.0-a-sqrt(-2.0*a*(1.0-f)) 
        xb1=(1.0-b1)/a 
        b=b1*h 
        xb=xb1*h 
        if(s.le.1.0-xb1) then 
          flux1=s*(a+b1-a*s/2.0)    !0.5*a*(1.0-xb1*xb1)+b1*(1.0-xb1) 
        else 
          flux1=f-1.0+s   !0.5*(2.0-a-b1)*(1.0-xb1)+(s-1.0+xb1) 
        endif 
        if(xb1.ge.s) then 
          flux2=s 
        else 
          flux2=xb1+0.5*a*(s*s-xb1*xb1)+b1*(s-xb1) 
        endif     
      else if(f.le.1.0-1.0/(2.0*beta).and.f.ge.1.0/(2.0*beta)) then     
        b1=(f-0.5*a) 
        b=b1*h 
        flux1=s*(a+b1-a*s/2.0)  !0.5*a*(2.0+s)*s+b1*s 
        flux2=s*(0.5*a*s+b1) 
      else if(f.le.1.0-beta/2.0.and.f.ge.beta/2.0) then 
        b1=0.5-a*f 
        b=b1*h 
        xt1=(1.0-b1)/a 
        xt=xt1*h 
        xb1=-b1/a 
        xb=xb1*h 
        if(s.le.1.0-xb1) then 
          flux1=0 
        else if(s.le.1.0-xt1) then 
          flux1=0.5*(xb1-(1.0-s))*(a*(1.0-s)+b1) 
        else  
          flux1=f-1.0+s 
        endif 
        if(s.le.xt1) then 
          flux2=s 
        else if(s.le.xb1) then 
          flux2=f-0.5*(xb1-s)*(a*s+b1) 
        else 
          flux2=f 
        endif     
      else            
        b1=sqrt(-2.0*a*f) 
        b=b1*h 
        xb1=-b1/a 
        xb=xb1*h 
        if(s.le.1.0-xb1) then 
          flux1=0 
        else  
          flux1=0.5*a*(xb1*xb1-(1.0-s)*(1.0-s))+b1*(xb1-(1.0-s)) 
        endif 
        if(s.le.xb1) then 
          flux2=(0.5*a*s+b1)*s 
        else 
          flux2=f 
        endif 
      endif 
      end  
	  
	  
SUBROUTINE VFLAIR(cycles,imax, jmax,wW, wE, wN, wS, delx, dely,delt, &
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
REAL(RP) s,flux,fa,fd,fb,slope1,slope,slope2,flux1,flux2,beta
integer I,J,id,ia,jd,ja

REAL FSN(0:IMAX+1,0:JMAX+1)
REAL(RP), DIMENSION(:,:), ALLOCATABLE ::FS1,FS0
REAL(RP),PARAMETER::EMF=1E-6,EMF1=1.0-EMF
ALLOCATE (FS1(0:IMAX+1,0:JMAX+1))

DO J=0,jmax
 DO I=0,imax
   FS1(I,J)=FS(I,J)   
 END DO
END DO

!c     the fluid transportation in X-direction 
        do i=1,imax 
          do j=1,jmax 
            flux=0 
            s=abs(u(i,j)*delt/delx) 
            if(u(i,j).ge.0) then 
              id=i   
              ia=i+1    
            else 
              id=i+1    
              ia=i 
            endif 
            if(fs(id,j).ge.1.0-emf) then 
              flux=s        
            else if(fs(i,j).le.1.0-emf.and.fs(i,j).ge.emf &
            .and.fs(i+1,j).le.1.0-emf.and.fs(i+1,j).ge.emf)then  
              if(fs(i,j).ge.fs(i+1,j)) then 
                  fa=fs(i,j) 
                  fb=fs(i+1,j)   
                  call calflux(fa,fb,u(i,j),delx,delt,flux) 
              else    
                  fa=fs(i+1,j) 
                  fb=fs(i,j)   
                  call calflux(fa,fb,-u(i,j),delx,delt,flux) 
              endif   
            else if((fs(id,j).le.1.0-emf.and.fs(id,j).ge.emf) & 
             .and.(fs(ia,j).le.emf.or.fs(ia,j).ge.1.0-emf)) then 
              if(fs(id,j+1).le.1.0-emf.and.fs(id,j+1).ge.emf)then  
                if(fs(id,j+1).ge.fs(id,j)) then 
                  fa=max(min(fs(id,j+1),1.0-emf),emf) 
                  fb=max(min(fs(id,j),1.0-emf),emf) 
                  call calslope(fa,fb,slope1) 
                else 
                  fa=max(min(fs(id,j),1.0-emf),emf) 
                  fb=max(min(fs(id,j+1),1.0-emf),emf) 
                  call calslope(fa,fb,slope1) 
                  slope1=-slope1 
                endif  
              endif 
              if(fs(id,j-1).le.1.0-emf.and.fs(id,j-1).ge.emf)then 
                if(fs(id,j).ge.fs(id,j-1)) then 
                  fa=max(min(fs(id,j),1.0-emf),emf) 
                  fb=max(min(fs(id,j-1),1.0-emf),emf) 
                  call calslope(fa,fb,slope2) 
                else 
                  fa=max(min(fs(id,j-1),1.0-emf),emf) 
                  fb=max(min(fs(id,j),1.0-emf),emf) 
                  call calslope(fa,fb,slope2) 
                  slope2=-slope2    
                endif 
              endif 
              if(fs(id,j+1).ge.1.0-emf.and.fs(id,j-1).ge.1.0-emf &
               .or.fs(id,j+1).le.emf.and.fs(id,j-1).le.emf) then 
                slope=0           
              else if(fs(id,j+1).ge.1.0-emf.and.fs(id,j-1).le.emf & 
               .or.fs(id,j+1).le.emf.and.fs(id,j-1).ge.1.0-emf) &
                 then 
                slope=0.5 
              else if(fs(id,j+1).ge.1.0-emf.or.fs(id,j+1).le.emf) &
                      then 
                slope=slope2 
              else if(fs(id,j-1).ge.1.0-emf.or.fs(id,j-1).le.emf) &
                      then 
                slope=slope1 
              else 
                slope=0.5*(slope1+slope2) 
              endif                 
              slope=0.5*(slope1+slope2) 
              if(abs(slope).le.emf) then 
                if(fs(ia,j).le.emf) then 
                  if((1.0-fs(id,j)).ge.s) then 
                    flux=0 
                  else 
                    flux=s-(1.0-fs(id,j)) 
                  endif 
                else 
                  if(fs(id,j).ge.s) then 
                    flux=s 
                  else 
                    flux=fs(id,j) 
                  endif 
                endif 
              else     
                beta=abs(slope)                      
                call cal2flux(fs(id,j),delx,u(i,j),delt,beta,flux1,flux2) 
                if(fs(ia,j).le.emf) then 
                  flux=flux1 
                else 
                  flux=flux2 
                endif     
              endif               
            endif   
            fs1(id,j)=fs1(id,j)-flux 
            fs1(ia,j)=fs1(ia,j)+flux 
          enddo 
        enddo 
!c     the transportation of Flux in Y-direction 
       do i=1,imax 
         do j=1,jmax 
           fs(i,j)=fs1(i,j) 
         enddo 
       enddo     
       do i=1,imax 
         do j=1,jmax    
            flux=0 
            s=abs(v(i,j)*delt/dely) 
            if(v(i,j).ge.0) then 
              jd=j   
              ja=j+1    
            else 
              jd=j+1    
              ja=j 
            endif 
            if(fs(i,jd).ge.1.0-emf) then 
              flux=s 
            else if(fs(i,j).le.1.0-emf.and.fs(i,j).ge.emf &
             .and.fs(i,j+1).le.1.0-emf.and.fs(i,j+1).ge.emf) then 
              if(fs(i,j).ge.fs(i,j+1)) then 
                  fa=fs(i,j) 
                  fb=fs(i,j+1)   
                  call calflux(fa,fb,v(i,j),dely,delt,flux) 
              else    
                  fa=fs(i,j+1) 
                  fb=fs(i,j)   
                  call calflux(fa,fb,-v(i,j),dely,delt,flux) 
              endif 
            else if((fs(i,jd).le.1.0-emf.and.fs(i,jd).ge.emf) & 
             .and.(fs(i,ja).le.emf.or.fs(i,ja).ge.1.0-emf)) then   
              if(fs(i+1,jd).le.1.0-emf.and.fs(i+1,jd).ge.emf)then  
                if(fs(i+1,jd).ge.fs(i,jd)) then 
                  fa=max(min(fs(i+1,jd),1.0-emf),emf) 
                  fb=max(min(fs(i,jd),1.0-emf),emf) 
                  call calslope(fa,fb,slope1) 
                else 
                  fa=max(min(fs(i,jd),1.0-emf),emf) 
                  fb=max(min(fs(i+1,jd),1.0-emf),emf) 
                  call calslope(fa,fb,slope1) 
                  slope1=-slope1 
                endif   
              endif 
              if(fs(i-1,jd).le.1.0-emf.and.fs(i-1,jd).ge.emf)then 
                if(fs(i,jd).ge.fs(i-1,jd)) then 
                  fa=max(min(fs(i,jd),1.0-emf),emf) 
                  fb=max(min(fs(i-1,jd),1.0-emf),emf) 
                  call calslope(fa,fb,slope2) 
                else 
                  fa=max(min(fs(i-1,jd),1.0-emf),emf) 
                  fb=max(min(fs(i,jd),1.0-emf),emf) 
                  call calslope(fa,fb,slope2) 
                  slope2=-slope2 
                endif 
              endif   
              if(fs(i+1,jd).ge.1.0-emf.and.fs(i-1,jd).ge.1.0-emf &
               .or.fs(i+1,jd).le.emf.and.fs(i-1,jd).le.emf) then 
                slope=0    
              else if(fs(i+1,jd).ge.1.0-emf.and.fs(i-1,jd).le.emf &
               .or.fs(i+1,jd).le.emf.and.fs(i-1,jd).ge.1.0-emf) then 
                slope=0.5 
              else if(fs(i+1,jd).ge.1.0-emf.or.fs(i+1,jd).le.emf) then 
                slope=slope2 
              else if(fs(i-1,jd).ge.1.0-emf.or.fs(i-1,jd).le.emf) then 
                slope=slope1 
              else 
                slope=0.5*(slope1+slope2)   
              endif 
              slope=0.5*(slope1+slope2)     
              if(abs(slope).le.emf) then 
                if(fs(i,ja).le.emf) then 
                  if(1.0-fs(i,jd).ge.s) then 
                    flux=0 
                  else 
                    flux=s-(1.0-fs(i,jd)) 
                  endif 
                else 
                  if(fs(i,jd).ge.s) then 
                    flux=s 
                  else 
                    flux=fs(i,jd) 
                  endif 
                endif 
              else 
                beta=abs(slope) 
                call cal2flux(fs(i,jd),dely,v(i,j),delt,beta,flux1,flux2) 
                if(fs(i,ja).le.emf) then 
                  flux=flux1 
                else 
                  flux=flux2 
                endif     
              endif 
            endif 
            fs1(i,jd)=fs1(i,jd)-flux 
            fs1(i,ja)=fs1(i,ja)+flux   
          enddo      
        enddo  
              
        do i=1,imax 
          do j=1,jmax   
            fs(i,j)=min(1.0,max(fs1(i,j),0.0)) 
          enddo 
        enddo     
	
DEALLOCATE (FS1)
END SUBROUTINE VFLAIR                              