c  This file is part of Fbmns3d
c  Copyright (C) 2004 Mourad Ismail
c
c  Fbmns3d is free software; you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation; either version 3 of the Licence, or 
c  (at your option) any later version.
c
c  Fbmns3d is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c
c  You should have received a copy of the GNU General Public License
c  along with Fbmns3d. If not, see <http://www.gnu.org/licenses/>.


c Author : Mourad Ismail   (ismail@ujf-grenoble.fr)
c $Id$
C======================================================================
      SUBROUTINE INITIAL (nx,ny,nz,nmx,nmy,nmz,npt,
     >         R,x,y,z,DnV,DnV0,Uinter,U_n,V_X0,V_Y0,V_Z0,PR0)
C
      IMPLICIT NONE
      INTEGER i,j,k,npt,nmx,nmy,nmz,nx,ny,nz
      DOUBLE PRECISION Uinter(*),DnV(*),DnV0(*),U_n(nmz,nmy,nmx),
     >     x(*),y(*),z(*),R,rho,V_X0(nmz,nmy,nmx),V_Y0(nmz,nmy,nmx),
     >     V_Z0(nmz,nmy,nmx),PR0(nmz,nmy,nmx)
C
      DO i=1,2*npt
         Uinter (i) = 0.D0
      ENDDO
      DO i=1,npt
         Dnv(i)=0.D0
         DnV0(i)=1.D0
      ENDDO
      do k=1,nz
       do j=1,ny
        do i=1,nx
          U_n (k,j,i)=0.d0
          if (k.eq.1) then
             PR0 (k,j,i)=1.d0
          else
             PR0 (k,j,i)=0.d0
          endif
          V_x0(k,j,i)=0.d0
          V_Y0(k,j,i)=0.d0
          V_Z0(k,j,i)=1.d0!*(1.D0-x(i)**2)*(1.D0-y(j)**2)
!               if (k.eq.1) then
!                  V_Z0(k,j,i)=1.d0
!               else
!                  V_Z0(k,j,i)=0.d0
!               endif
        enddo
       enddo
      enddo
c$$$         do j=1,ny
c$$$            do i=1,nx
c$$$          V_x0(nz,j,i)=1.d0
c$$$            EndDo
c$$$         EndDo
C
      RETURN
      END
C
*=======================================================================*
      SUBROUTINE input_2nd_mbre(switch,ll,nx,ny,nz,nmx,nmy,nmz,np,
     >i_temps,dt,nu,x,y,z,ax,ay,az,R,rho,f)
*=======================================================================*
      IMPLICIT NONE
      INTEGER nx,ny,nz,nmx,nmy,nmz,np,i,j,k,l,i_temps,npp,ll
      double precision R,dt,temps,nu
      DOUBLE PRECISION x(*),y(*),z(*)
      double precision f(nmz,nmy,nmx)
      double precision ax(*),ay(*),az(*),rho(*)
      logical switch
C
      temps=i_temps*dt
      if (switch) then
C-PRESSION--PRESSION--PRESSION--PRESSION--PRESSION--PRESSION--PRESSION-
      do k=1,nz  
         do j=1,ny  !A MODIFIER POUR TENIR COMPTE DES NOUVEAUX DDL DE LA PRESSION
            do i=1,nx
               f(k,j,i)=0.d0 
            EndDo
         EndDo
      end do
      else
C-Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse-
      if(ll.eq.1) then
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f(k,j,i)=0.D0
            EndDo
         EndDo
      endDo
      elseif(ll.eq.2) then
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f(k,j,i)=0.D0
            EndDo
         EndDo
      endDo
      elseif(ll.eq.3) then
      do k=1,nz
         do j=1,ny
            do i=1,nx
               f(k,j,i)=0.d0
            EndDo
         EndDo
      endDo
      else
         STOP
      endif
      endif
C
      RETURN
      END
*=======================================================================*
      SUBROUTINE input_analy(switch,l,nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     >      x,y,z,R,u)
*=======================================================================*
      IMPLICIT NONE
      INTEGER nmx,nmy,nmz,i,j,k,l,i_temps,nx,ny,nz
      double precision dt,temps,R,pi
      DOUBLE PRECISION x(*),y(*),z(*),u(nmz,nmy,nmx)
      logical switch
c
      temps=i_temps*dt ; pi=dacos(-1.d0)
c
      if (switch) then
C-PRESSION--PRESSION--PRESSION--PRESSION--PRESSION--PRESSION--PRESSION-
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if(k.eq.1) then 
                  u(k,j,i)=1.d0
!                  u(k,j,i)=0.d0
               else
                  u(k,j,i)=0.d0
               endif
            end do
         end do
      end do
      else
C-Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse--Vitesse-
      if(l.eq.1) then      !--------VX--------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
          u(k,j,i)=0.d0
            EndDo
         EndDo
      endDo
c$$$         do j=1,ny
c$$$            do i=1,nx
c$$$          u(nz,j,i)=1.d0
c$$$            EndDo
c$$$         EndDo
      elseif(l.eq.2) then  !--------VY--------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
          u(k,j,i)=0.d0
            EndDo
         EndDo
      endDo
      elseif(l.eq.3) then  !--------VZ--------------
      do k=1,nz
         do j=1,ny
            do i=1,nx
!               u(k,j,i)=1.d0*(x(i)**2-1.d0)*(y(j)**2-1.d0)
!          u(k,j,i)=16.d0*x(i)*(1.D0-x(i))*y(j)*(1.D0-y(j))
!               if (k.eq.1) then
                  u(k,j,i)=1.d0
 !              else
!                  u(k,j,i)=0.d0
!               endif
            EndDo
         EndDo
      endDo
      else                  !--------ERROR----------
         STOP
      endif
      endif
C
      RETURN
      END
c
c========================================================================c
      subroutine Input_func(d,n,nmx,nmy,nmz,r,x,y,z,ix,iy,iz,a,b,c,
     >                      U,ga)
c========================================================================c
c
      implicit none
      integer i,j,l,n,d,nmx,nmy,nmz,ix(8),iy(8),iz(8)
      double precision r,r2,rho,rhoR2,dxk,dyk,dzk,x(*),y(*),z(*)
      double precision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z,duex_x
      double precision duex_y,duex_z,uex,ucomp,interpol_Q1,ga(5,5,5)
      double precision a(5),b(5),c(5),U(nmz,nmy,nmx)
      double precision pi, beta, beta2
c
      R2=r**2 ; pi=dacos(-1.d0) ; beta=2.d0*pi ; beta2=beta**2
c
      if (d.eq.0) then ! d=0 ----> input exact and computing solutions for 
c                      !            L2 error computing 
      do l=1,n
       do j=1,n
         do i=1,n
           rho=a(i)**2+b(j)**2+c(l)**2
           rhoR2=rho-(R2)         
         if(rho.le.R2)then
       uex  =1.d0
       ucomp=interpol_Q1(nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)
          else  
!       uex=dcos(beta*(rhoR2))
!       ucomp=interpol_Q1(nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)
         endif
           ga(i,j,l)=(uex-ucomp)**2
         enddo
       enddo
      enddo
c
      elseif (d.eq.1) then ! d=1 ----> input exact and computing solutions for 
c                          !            H1 error computing  
      do l=1,n
       do j=1,n
         do i=1,n
           rho=a(i)**2+b(j)**2+c(l)**2
           rhoR2=rho-R2
      dxk=interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)
      dyk=interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)
      dzk=interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)
         if(rho.le.R2)then
           duex_x =  0.d0
           duex_y =  0.d0
           duex_z =  0.d0
          else 
           duex_x = -2*beta*a(i)*dsin(beta*(rhoR2))
           duex_y = -2*beta*b(j)*dsin(beta*(rhoR2))
           duex_z = -2*beta*c(l)*dsin(beta*(rhoR2))
!           duex_x = 0.5d0*a(i)*(2.d0*(a(i)**2)+c(l)**2-2.d0*R2)
!           duex_y =-0.5d0*b(j)*(2.d0*(b(j)**2)+c(l)**2-2.d0*R2)
!           duex_z = 0.5d0*c(l)*(a(i)**2-b(j)**2)
         endif
           ga(i,j,l)=(duex_x-dxk)**2+(duex_y-dyk)**2+(duex_z-dzk)**2
         enddo
       enddo
      enddo
c
      elseif (d.eq.2) then ! d=2 ----> input right hand side
c
      do l=1,n
       do j=1,n
        do i=1,n
          rho=a(i)**2+b(j)**2+c(l)**2
!          if(rho.le.R2) then
           ga(i,j,l)=0.d0
!           else
!           ga(i,j,l)= 6.d0*beta*dsin(beta*(rho-R2)) +
!     >                4.d0*beta2*rho*dcos(beta*(rho-R2))
!           ga(i,j,l)=7.d0*(b(j)**2-a(i)**2)/2.d0
!          endif
        enddo
       enddo
      enddo
c                          
      else
       print*,' '
       print*,'ABORTING.......'
       print*,'Error in d index. d=',d
       print*,'d must be :'
       print*,'     0 for input exact and computing solutions'
       print*,'     1 for input exact and computing solutions gradient'
       print*,'     2 for input right hand side'
       print*,' '
       stop
      endif
c
      return
      end
c
C=====================================================================
      SUBROUTINE DATA_AND_GEOM(ncs,npt,np,nr,r,epsilon,theta,alpha,nu,
     > zeta,dt,itermaxV,itermax,itermax_tp,prec,xxi,yyi,zzi,xl,yl,zl,
     > bcP,bcVX,bcVY,bcVZ,ax,ay,az,switch)
C=====================================================================
      IMPLICIT NONE
      INTEGER i,np,nr,ncs,npt,itermax,itermax_tp,itermaxV
      integer bcP(6),bcVX(6),bcVY(6),bcVZ(6)
      DOUBLE PRECISION r,epsilon,theta,alpha,xxi,yyi,zzi,xl,yl,zl,
     >        prec,nu,dt,zeta
      DOUBLE PRECISION ax(*),ay(*),az(*)
      logical switch
C
      open(unit=29,file='data.data',status='unknown')
      read(29,*)
      read(29,*) xxi , yyi , zzi , xl , yl , zl , switch
      read(29,*)
      read(29,*) (bcVX(i),i=1,6)
      read(29,*) (bcVY(i),i=1,6)
      read(29,*) (bcVZ(i),i=1,6)
      read(29,*) (bcP (i),i=1,6)
      read(29,*)
      read(29,*) itermaxV, itermax, prec
      read(29,*)
      read(29,*) alpha , zeta, theta , epsilon , r , np , nr
      read(29,*)
      read(29,*) nu , dt ,itermax_tp
      read(29,*)      
c
      if (np.ne.0) then
        DO i=1,np
         read(29,*) ax(i),ay(i),az(i)
        ENDDO
        npt=np*ncs
       else
        npt=ncs
      endif
c
      do i=1,6
       if((bcP (i).ne.1).and.(bcP (i).ne.0).or.
     >    (bcVX(i).ne.1).and.(bcVX(i).ne.0).or.
     >    (bcVY(i).ne.1).and.(bcVY(i).ne.0).or.
     >    (bcVZ(i).ne.1).and.(bcVZ(i).ne.0) ) then
        print*,'Error in Boundary Conditions index'
        stop
       endif
      enddo
!      alpha=alpha/(nu*dt)
      CLOSE(29)
C
      RETURN
      END
