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

C	$Id$	


C--------------------------------------------------------------------C
C-------------------------FUNCTIONS----------------------------------C
C--------------------------------------------------------------------C
C
c*********************************************************************
c
c Function time_cpu
c
c Purpose:
c   returns the used CPU-time in seconds.
c
c Parameters:
c
c Output:
c   time_cpu - used cpu time in seconds.
c
C======================================================================C
      double precision function time_cpu()
C======================================================================C
      implicit none
      real etime, t(2)
c
      time_cpu = dble(etime(t))
      return
      end
C======================================================================C
      FUNCTION prodsca(ns,H,ARH)
C======================================================================C      
      IMPLICIT NONE
      INTEGER ns,i
      DOUBLE PRECISION prodsca,H(*),ARH(*)
C
      prodsca=0.D0
      DO i=1,ns
      prodsca=prodsca+H(i)*ARH(i)
      ENDDO
C
      RETURN
      END
C
C=======================================================================*
      FUNCTION prodsca_3d(nmx,nmy,nmz,nx,ny,nz,f,g)
C=======================================================================*      
      IMPLICIT NONE
      INTEGER nmx,nmy,nmz,nx,ny,nz,i,j,k
      DOUBLE PRECISION prodsca_3d,f(nmz,nmy,nmx),g(nmz,nmy,nmx)
C
      prodsca_3d=0.D0
      DO k=1,nz
       DO j=1,ny
         DO i=1,nx
            prodsca_3d=prodsca_3d+f(k,j,i)*g(k,j,i)
         ENDDO
       ENDDO
      ENDDO
C
      RETURN
      END
C
C======================================================================C
C
C--------------------------------------------------------------------C
C--------------------------ROUTINES----------------------------------C
C--------------------------------------------------------------------C
C
C======================================================================C
      subroutine reorganization (nx_fd,ny_fd,nz_fd,nx,ny,nz,nmx,nmy,nmz,
     > bc,frestore,u,f)
C======================================================================C
C
      implicit none
      integer nmx,nmy,nmz,nx_fd,ny_fd,nz_fd,i,j,k,bc(6)
      integer nx,ny,nz
      double precision f(nmz,nmy,nmx),frestore(nmz,nmy,nmx),
     > u(nmz,nmy,nmx)
C
      do k=1,nz_fd
       do j=1,ny_fd
        do i=1,nx_fd
          frestore(k,j,i)=f(k,j,i)
        enddo
       enddo
      enddo

      if(bc(5).eq.1) then
      do j=1,ny
       do i=1,nx
         f(1,j,i)=u(1,j,i)
       enddo
      enddo
      endif
C
      if(bc(6).eq.1) then
      do j=1,ny
       do i=1,nx
         f(nz,j,i)=u(nz,j,i)
       enddo
      enddo
      endif
C
      if(bc(3).eq.1) then
      do k=1,nz
       do i=1,nx
        f(k,1,i)=u(k,1,i)
       enddo
      enddo      
      endif
C
      if(bc(4).eq.1) then
      do k=1,nz
       do i=1,nx
        f(k,ny,i)=u(k,ny,i)
       enddo
      enddo
      endif
C
      if(bc(1).eq.1) then
      do k=1,nz
       do j=1,ny
        f(k,j,1)=u(k,j,1)
       enddo
      enddo
      endif
C
      if(bc(2).eq.1) then
      do k=1,nz
       do j=1,ny
        f(k,j,nx)=u(k,j,nx)
       enddo
      enddo
      endif
C
      do k=1,nz_fd
       do j=1,ny_fd
        do i=1,nx_fd
         f(k+bc(5),j+bc(3),i+bc(1))=frestore(k,j,i)
        enddo
       enddo
      enddo
C
      return 
      end
C======================================================================C
      subroutine elimination (nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,
     >                        b1,d1,b2,d2,b3,d3,bc,f)
C======================================================================C
C
      implicit none
      integer nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,i,j,k,bc(6)
      doubleprecision f(nmz,nmy,nmx),b1(*),d1(*),b2(*),d2(*),b3(*),d3(*)
C
       do k=1,nz_fd
        do j=1,ny_fd
         do i=1,nx_fd
          f(k,j,i)=f(k+bc(5),j+bc(3),i+bc(1))
         enddo
        enddo
       enddo
C
      if((bc(1).eq.1).and.(bc(2).eq.0)) then
      b1(nx_fd) = b1(nx_fd)/2
      d1(nx_fd) = d1(nx_fd)/2
      endif
      if((bc(3).eq.1).and.(bc(4).eq.0)) then
      b2(ny_fd) = b2(ny_fd)/2
      d2(ny_fd) = d2(ny_fd)/2
      endif
      if((bc(5).eq.1).and.(bc(6).eq.0)) then
      b3(nz_fd) = b3(nz_fd)/2
      d3(nz_fd) = d3(nz_fd)/2
      endif

C
      return 
      end
C======================================================================C
      subroutine Update (nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,a,b,f,U_n)
C======================================================================C
C
      implicit none
      integer nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,i,j,k
      double precision a,b,f(nmz,nmy,nmx),U_n(nmz,nmy,nmx)
C
      do k=1,nz_fd
       do j=1,ny_fd
        do i=1,nx_fd
         u_n(k,j,i)=(a*u_n(k,j,i))+(b*f(k,j,i))
        enddo
       enddo
      enddo
C
      return 
      end
C======================================================================C
      subroutine Degree_Freedom (nx,ny,nz,nx_fd,ny_fd,nz_fd,
     > xl,yl,zl,bc,dbc)
C======================================================================C
C
      implicit none
      integer nx,ny,nz,nx_fd,ny_fd,nz_fd,dbc,bc(6)
      double precision xl,yl,zl
C
       nx_fd = nx - (bc(1)+bc(2))
       ny_fd = ny - (bc(3)+bc(4))
       nz_fd = nz - (bc(5)+bc(6))
      dbc=bc(1)+bc(2)+bc(3)+bc(4)+bc(5)+bc(6)
      print*,'**********************************************'
      print*,' nx_fd=',nx_fd,' ny_fd=',ny_fd,' nz_fd=',nz_fd
      print*,' '
      print*,'It will be better that this values ares in this form:'
      print*,'2^^^^^^^^^^^^^^^q-----------------1'
      print*,' '
      print*,' hx=xl/',nx-1,' hy=yl/',ny-1,' hz=zl/',nz-1
      print*,' hx =',xl/(2.D0*(nx-1))!,' hy =',yl/(ny-1),' hz =',zl/(nz-1)
      print*,'**********************************************'
!      print*,'Are this values ok ?' 
!      pause
C
      return 
      end
c
C====================================================================================C
      subroutine Coord (ik1,nx,ny,ix,iy,iz)
C====================================================================================C
c
      implicit none
      integer nx,ny,i,ik1,npxy,ninter
      integer ix(8),iy(8),iz(8)
c
      npxy=nx*ny 
c
      iz(1) = (ik1/npxy)+1
      if (mod(ik1,npxy).eq.0) iz(1) = iz(1)-1
      ninter = ik1 - ((iz(1)-1)*npxy)
      iy(1) = (ninter/nx)+1
      if (mod( ninter,nx).eq.0)  iy(1) = iy(1)-1
      ix(1) = ninter - (iy(1)-1)*nx
c
      do i=1,3
       iz(i+1)=iz(1) 
       iz(i+4)=iz(1)+1
      enddo
       iz(8)=iz(1)+1
c
      iy(2) = iy(1)   ; iy(5) = iy(1); iy(6) = iy(1)
      iy(3) = iy(1)+1 ; iy(4) = iy(3); iy(7) = iy(3) ; iy(8) = iy(3) 
c
      ix(4) = ix(1)   ; ix(5) = ix(1); ix(8) = ix(1)
      ix(2) = ix(1)+1 ; ix(3) = ix(2); ix(6) = ix(2) ; ix(7) = ix(2) 
c
      return
      end

C======================================================================C
      subroutine n_vect3(OA,OB,u)
! produit vectoriel normalise. u = (OA x OB)/||OA x OB||
C======================================================================C      
      IMPLICIT NONE
      INTEGER i
      DOUBLE PRECISION OA(*),OB(*),u(*),norm
c
      u(1) = OA(2)*OB(3) - OB(2)*OA(3)
      u(2) = OB(1)*OA(3) - OA(1)*OB(3)      
      u(3) = OA(1)*OB(2) - OB(1)*OA(2)
c     
      norm = dsqrt(u(1)**2 + u(2)**2 + u(3)**2)
c
      do i=1,3
       u(i)=u(i)/norm
      enddo
c
      return
      end
C======================================================================C

