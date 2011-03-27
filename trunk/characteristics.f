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
C======================================================================C
      subroutine characteristics (nx,ny,nz,nmx,nmy,nmz,
     > V_x,V_y,V_z,xi,yi,zi,xl,yl,zl,hx,hy,hz,dt,x,y,z,num,U,U_conv)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j,k,i,ind_elt,i_conv,j_conv,k_conv,
     >        nx,ny,nz,ik,ix(8),iy(8),iz(8),num(*)
      double precision x_conv,y_conv,z_conv,dt,hx,hy,hz,
     >  interpol_Q1,x(*),y(*),z(*),u(nmz,nmy,nmx),U_conv(nmz,nmy,nmx),
     >  xi,xl,xf,yi,yl,yf,zi,zl,zf,
     >v_x(nmz,nmy,nmx),v_y(nmz,nmy,nmx),v_z(nmz,nmy,nmx)
C
      xf=xi+xl ; yf=yi+yl ; zf=zi+zl 
      do k=1,nz
        do j=1,ny
          do i=1,nx
             !abscisse du point convecte
             x_conv = x(i)-dt*V_x(k,j,i)
             !ordonnee du point convecte
             y_conv = y(j)-dt*V_y(k,j,i)
             !cote du point convecte
             z_conv = z(k)-dt*V_z(k,j,i)
      ! si la caracteristique sort du domaine, garder l'ancien U_conv (D_BC)       
      IF (x_conv.LE.xi) then
         U_conv(k,j,i)=u(k,j,1)
        elseif(x_conv.GE.xf)then
         U_conv(k,j,i)=u(k,j,nx)
        elseif(y_conv.LE.yi)then
         U_conv(k,j,i)=u(k,1,i)
        elseif(y_conv.GE.yf)then
         U_conv(k,j,i)=u(k,ny,i)
        elseif(z_conv.LE.zi)then
         U_conv(k,j,i)=u(1,j,i)
        elseif(z_conv.GE.zf)then
         U_conv(k,j,i)=u(nz,j,i)
         else
            !sa localisation dans le maillage 
             i_conv=int(dabs(x_conv-xi)/hx)+1
             j_conv=int(dabs(y_conv-yi)/hy)+1
             k_conv=int(dabs(z_conv-zi)/hz)+1
      ind_elt=((k_conv-1)*(ny-1)*(nx-1))+((j_conv-1)*(nx-1))+i_conv
        ik=num(ind_elt)
      call Coord (ik,nx,ny,ix,iy,iz)    
      !U_conv(X)=U(X-dt*U)
      U_conv(k,j,i)=
     >       interpol_Q1(nmx,nmy,nmz,ix,iy,iz,
     >                     x_conv,y_conv,z_conv,x,y,z,U)
      endif
      enddo 
        enddo 
      enddo 
C
      return 
      end
C======================================================================C
      subroutine charact_Vect (nx,ny,nz,nmx,nmy,nmz,xi,yi,zi,xl,yl,zl,
     >  hx,hy,hz,dt,x,y,z,num,V_x,V_y,V_z,UC_x,UC_y,UC_z)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j,k,i,ind_elt,i_conv,j_conv,k_conv,
     >        nx,ny,nz,ik,ix(8),iy(8),iz(8),num(*)
      double precision x_conv,y_conv,z_conv,dt,hx,hy,hz,
     >  interpol_Q1,x(*),y(*),z(*),xi,xl,xf,yi,yl,yf,zi,zl,zf,
     >v_x(nmz,nmy,nmx),v_y(nmz,nmy,nmx),v_z(nmz,nmy,nmx),
     >UC_x(nmz,nmy,nmx),UC_y(nmz,nmy,nmx),UC_z(nmz,nmy,nmx)
C
      print*, ' '
      print*,' Convecting the Velocity by the characteristics...'
C
      xf=xi+xl ; yf=yi+yl ; zf=zi+zl 
      do k=1,nz
        do j=1,ny
          do i=1,nx
             !abscisse du point convecte
             x_conv = x(i)-dt*V_x(k,j,i)
             !ordonnee du point convecte
             y_conv = y(j)-dt*V_y(k,j,i)
             !cote du point convecte
             z_conv = z(k)-dt*V_z(k,j,i)
      ! si la caracteristique sort du domaine, garder l'ancien U_conv (D_BC)       
      IF (x_conv.LE.xi) then
         UC_x(k,j,i)=V_x(k,j,1)
         UC_y(k,j,i)=V_y(k,j,1)
         UC_Z(k,j,i)=V_z(k,j,1)
        elseif(x_conv.GE.xf)then
         UC_x(k,j,i)=V_x(k,j,nx)
         UC_y(k,j,i)=V_y(k,j,nx)
         UC_Z(k,j,i)=V_z(k,j,nx)
        elseif(y_conv.LE.yi)then
         UC_x(k,j,i)=V_x(k,1,i)
         UC_y(k,j,i)=V_y(k,1,i)
         UC_Z(k,j,i)=V_z(k,1,i)
        elseif(y_conv.GE.yf)then
         UC_x(k,j,i)=V_x(k,ny,i)
         UC_y(k,j,i)=V_y(k,ny,i)
         UC_Z(k,j,i)=V_z(k,ny,i)
        elseif(z_conv.LE.zi)then
         UC_x(k,j,i)=V_x(1,j,i)
         UC_y(k,j,i)=V_y(1,j,i)
         UC_Z(k,j,i)=V_z(1,j,i)
        elseif(z_conv.GE.zf)then
         UC_x(k,j,i)=V_x(nz,j,i)
         UC_y(k,j,i)=V_y(nz,j,i)
         UC_Z(k,j,i)=V_z(nz,j,i)
         else
            !sa localisation dans le maillage 
             i_conv=int(dabs(x_conv-xi)/hx)+1
             j_conv=int(dabs(y_conv-yi)/hy)+1
             k_conv=int(dabs(z_conv-zi)/hz)+1
      ind_elt=((k_conv-1)*(ny-1)*(nx-1))+((j_conv-1)*(nx-1))+i_conv
        ik=num(ind_elt)
      call Coord (ik,nx,ny,ix,iy,iz)    
      !U_conv(X)=U(X-dt*U)
      UC_x(k,j,i)=
     >       interpol_Q1(nmx,nmy,nmz,ix,iy,iz,
     >                     x_conv,y_conv,z_conv,x,y,z,V_x)
      UC_y(k,j,i)=
     >       interpol_Q1(nmx,nmy,nmz,ix,iy,iz,
     >                     x_conv,y_conv,z_conv,x,y,z,V_y)
      UC_z(k,j,i)=
     >       interpol_Q1(nmx,nmy,nmz,ix,iy,iz,
     >                     x_conv,y_conv,z_conv,x,y,z,V_z)
      endif
      enddo 
        enddo 
      enddo 
c
      print*, ' '
      print*,' Characteristics.......DONE'
      print*, ' '
C
      return 
      end

