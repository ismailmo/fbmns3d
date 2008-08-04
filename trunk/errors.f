C======================================================================C
      Function L_inf_error_fd (nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,bc,f,U)
C======================================================================C
      implicit none
      integer nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,i,j,k,bc(6)
      double precision err,L_inf_error_fd
      double precision f(nmz,nmy,nmx),U(nmz,nmy,nmx)
C
      L_inf_error_fd = 0.d0
      do k=1,nz_fd
         do j=1,ny_fd
            do i=1,nx_fd
      err = dabs(f(k,j,i) - u(k+bc(5),j+bc(3),i+bc(1)))
      L_inf_error_fd = max(err, L_inf_error_fd)
            end do
         end do
      end do
C
      return 
      end
C======================================================================C
      Function Linfinity_error (nx,ny,nz,nmx,nmy,nmz,f,U)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,i,j,k,nx,ny,nz
      double precision err,Linfinity_error
      double precision f(nmz,nmy,nmx),U(nmz,nmy,nmx)
C
      Linfinity_error = 0.d0
      do k=1,nz
         do j=1,ny
            do i=1,nx
      err = dabs(f(k,j,i) - u(k,j,i))
      Linfinity_error = max(err, Linfinity_error)
            end do
         end do
      end do
C
      return 
      end
C======================================================================C
!
!      integ=dx*(w0*f(x0)+w1*(f(x-)+f(x+)))
!
!      integ=dx*(w0*
!                   (dy*(w0*f(x0,y0)+w1*(f(x0,y-)+f(x0,y+))))
!               +w1*(
!                    (dy*(w0*f(x-,y0)+w1*(f(x-,y-)+f(x-,y+)))) +
!                    (dy*(w0*f(x+,y0)+w1*(f(x+,y-)+f(x+,y+))))  ))
C======================================================================C
c   L2 and H1 errors computing
      Subroutine Errors (d,n,nx,ny,nz,nc,nmx,nmy,nmz,r,hx,hy,hz,
     >      num,x,y,z,U,errquad,errquadloc)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,nx,ny,nz,nc,i,j,k,l,n,d
      integer ix(8),iy(8),iz(8),num(*)
      double precision errquad,errquadloc,errquadfatloc,rho,rhoR2
      double precision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z
      double precision pi,r,r2,r3,beta,xmax,xmin,ymin,ymax,zmin,zmax
      double precision dxk,dyk,dzk,duex_x,duex_y,duex_z,hx,hy,hz,r4
      double precision rt,rt2,w0,w1,w2,dx,dy,dz,integ_QN,a(5),b(5),c(5)
      double precision u(nmz,nmy,nmx),ga(5,5,5),x(*),y(*),z(*),integ
c
      if (d.eq.0) then 
       print*,'Computing L2 errors for hx,hy,hz=',hx,hy,hz
       print*,'...'
      elseif(d.eq.1) then
       print*,'Computing H1 errors for hx,hy,hz=',hx,hy,hz
       print*,'...'
      endif
c
      pi=dacos(-1.d0)
      R2=R**2
      R3=(R+0.05D0)**2
      R4=(R+0.1D0)**2      
      beta=2.d0*pi

      errquad=0.D0 ;  errquadloc=0.D0 ; errquadfatloc=0.D0
c
      do k=1,nc
c
       call Coord (num(k),nx,ny,ix,iy,iz)

       xmin=x(ix(1)) ; xmax=x(ix(7))
       ymin=y(iy(1)) ; ymax=y(iy(7))
       zmin=z(iz(1)) ; zmax=z(iz(7))

       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
       call Input_func(d,n,nmx,nmy,nmz,r,x,y,z,ix,iy,iz,a,b,c,
     >                      U,ga)
c
       integ=integ_QN(n,dx,dy,dz,w0,w1,w2,ga)
       errquad=errquad+integ
c
       !square of the distance between (0,0,0) and the cube center
       rho=0.0625d0*( (a(1)+a(2)+a(n-1)+a(n))**2 +
     >                (b(1)+b(2)+b(n-1)+b(n))**2 +
     >                (c(1)+c(2)+c(n-1)+c(n))**2   )
c
       if(rho.gt.R4)then
        errquadfatloc=errquadfatloc+integ
        errquadloc=errquadloc+integ
       elseif(rho.gt.R3)then
         errquadloc=errquadloc+integ
       endif
c
      enddo
c
      return
      end
