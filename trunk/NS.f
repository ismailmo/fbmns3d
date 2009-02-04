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
c
c Author : Mourad Ismail   (ismail@ujf-grenoble.fr)
c
      program FBM_Navier_Stokes
C{
      implicit none
      integer nmx,nmy,nmz,nlmax,ldw,liw,ldwP,liwP
      integer nn,nxyz,px,py,pz,nx,ny,nz,ns,nc,pc
      integer npmax,nptmax,nss,ncs,nrmax,nr2max
c
c Parameters and data which user can define
c
      parameter (nmx=160,nmy=160,nmz=400)
c
      parameter (nxyz=19)   ! 2**q+1 ------> h=1/2**q
!      parameter (px = 129, py = 129, pz = 129)
      parameter (px = 19, py = 19, pz = 74)
!      parameter (px = 14, py = 17, pz = 19)
      parameter (nx = 2*px-1, ny = 2*py-1, nz = 2*pz-1, nn=30)
!      parameter (nx = px, ny = py, nz = pz, nn=30)
      parameter (npmax=56,nrmax=3) 
C     
      !parameter (ldw = 6*nlmax*(nx+ny) + max(9*nx,7*ny*nz) + max(9*ny,10*nz)) 
      parameter (nlmax=6 ,ldw = 6*nlmax*(nx+ny) + (7*ny*nz) +(10*nz))
      parameter (liw = 2*(4**nlmax-1)/3 + 4*nlmax + 5*nx + 5*ny + 7)
      parameter (         ldwp= 6*nlmax*(px+py) + (7*py*pz) +(10*pz))
      parameter (liwp= 2*(4**nlmax-1)/3 + 4*nlmax + 5*px + 5*py + 7)
c 
      parameter (ns=nx*ny*nz, nc=(nx-1)*(ny-1)*(nz-1))
      parameter (pc=(px-1)*(py-1)*(pz-1))
      parameter (nss=6*(nn-1)*(nn-1)+2, ncs=nss-2)
      parameter (nptmax=npmax*ncs,nr2max=2**nrmax)
c 
      integer iw(liw), ierr , iwP(liwp)
      double precision  dw(ldw),dwP(ldwp)
C
      integer num(nc),NumP(pc),bcVX(6),bcVY(6),bcVZ(6),bcP(6)
      double precision a1(nx), b1(nx), c1(nx), d1(nx)
      double precision a2(ny), b2(ny), c2(ny), d2(ny)
      double precision a3(nz), b3(nz), c3(nz), d3(nz)
      double precision x (nx), y (ny), z (nz)
      double precision xp(px), yp(py), zp(pz)
      double precision ax(npmax),ay(npmax),az(npmax),rho(npmax) 
c
      integer nums(8,ncs),NREP(2*ncs),ncb(npmax),numb(nc,npmax)
      integer nrepS(2*nss)
      double precision SS(3,2*nss),Area(ncs),Uinter(2*nptmax),
     > Uinter0(2*nptmax)
      double precision DnVX(nptmax),DnVX0(nptmax),DnVdiff(nptmax)
      double precision DnVY(nptmax),DnVY0(nptmax)
      double precision DnVZ(nptmax),DnVZ0(nptmax)
c
      double precision f(nmz,nmy,nmx), u(nmz,nmy,nmx),g(nmz,nmy,nmx)
      double precision Uo(nmz,nmy,nmx),U1(nmz,nmy,nmx),f0(nmz,nmy,nmx)
c
      double precision Linfinity_error,prodsca,prodsca_3d, Integ_G
c 
      double precision cpu , time_cpu
      double precision xi, yi, zi, xl, yl, zl, R, hx, hy, hz, hv
      double precision hpx, hpy, hpz, hpv
      double precision alpha, zeta, theta, epsilon, prec,  nu, dt, dtp
      integer np,npt,nr,itermax,itermaxP,itermax_tp,dbcP,dbcX,dbcY,dbcZ
      integer nx_fdX, ny_fdX, nz_fdX,nx_fdP, ny_fdP, nz_fdP
      integer nx_fdY, ny_fdY, nz_fdY,nx_fdZ, ny_fdZ, nz_fdZ
c     
      logical switch
      character*4 xfn*11
      character*4 yfn*11
      character*4 zfn*11
      character*4 pfn*11
      character*4 vfn*11
      character*4 cfn*11
      integer i_temps,i,k,j,l,ii,kk,p,nshowv,nshow,nshow2,ic,jc,kc
      double precision normednv,beta,pi
      double precision V_X0(nmz,nmy,nmx),V_Y0(nmz,nmy,nmx)
      double precision V_Z0(nmz,nmy,nmx),UC_x(nmz,nmy,nmx)
      double precision UC_y(nmz,nmy,nmx),UC_z(nmz,nmy,nmx)
      double precision PR0 (nmz,nmy,nmx),PR1 (nmz,nmy,nmx)
      double precision PRF(nmz,nmy,nmx)
      double precision bx(-1:100),by(-1:100),bz(-1:100)
      doubleprecision  Integ_GPX,Integ_GPY,Integ_GPZ,Pinter(2*nptmax)
      doubleprecision PinterS(2*nss),debitZF
      integer ptsForce(ncs)
      data vfn,cfn,xfn/'Vof_Vel.000','Vof_Coo.000','DX__V_X.000'/
      data yfn,zfn/'DX__V_Y.000','DX__V_Z.000'/
      data     pfn/'DX__Pre.000'/
C     } 
C------------------------------Geometrie--------------------------------C
C     Domaine = ]xxi,xxi+xl[x]yyi,yyi+yl[x]zzi,zzi+zl[
C     gamma = np spheres S(a,r) , gammaprime = np spheres S(a,r+epsilon)
C     a=(ax(np),ay(np),az(np))
C-----------------------------------------------------------------------C
C     
      pi=dacos(-1.d0) ; i_temps=0 ; nshowv=(nx/2) ;  nshow=(ncs/2)
      nshow2=2*ncs+(ncs/2)
c     
c     cpu = Time_Cpu() 
c     
      CALL DATA_AND_GEOM(ncs,npt,np,nr,R,epsilon,theta,alpha,nu,
     >     zeta,dt,itermaxP,itermax,itermax_tp,prec,xi,yi,zi,xl,yl,zl,
     >     bcP,bcVX,bcVY,bcVZ,ax,ay,az,switch)
      call Degree_Freedom (px,py,pz,nx_fdP,ny_fdP,nz_fdP,
     >     xl,yl,zl,bcP,dbcP)
      call Degree_Freedom (nx,ny,nz,nx_fdX ,ny_fdX ,nz_fdX ,
     >     xl,yl,zl,bcVX,dbcX )
      call Degree_Freedom (nx,ny,nz,nx_fdY ,ny_fdY ,nz_fdY ,
     >     xl,yl,zl,bcVY,dbcY )
      call Degree_Freedom (nx,ny,nz,nx_fdZ ,ny_fdZ ,nz_fdZ ,
     >     xl,yl,zl,bcVZ,dbcZ )
c     Maillages de OMEGA, gamma et omega
      CALL MSH_OMG (nx,ny,nz,hx,hy,hz,hv,xi,yi,zi,xl,yl,zl,x,y,z,NUM)  
      CALL MSH_OMG (px,py,pz,hpx,hpy,hpz,hpv,xi,yi,zi,xl,yl,zl,
     >     xp,yp,zp,NumP)  
c     
      if (np.ne.0) then
         call MSH_sphere (nn,nss,ncs,r,SS,NUMS)      
         CALL MSH_omega  (nss,ncs,epsilon,r,nums,SS)
c     Calcul des aires des elements de gamma
         call  area_spherical_quads (ncs,nums,SS,R,area)
c     CALL Surface_Gamma (ncs,nss,NUMS,SS,Area)
c     Localisation des sommets de gamma dans le maillage global
         CALL Locate_G_Gp(nx,ny,nz,hx,hy,hz,xi,yi,zi,ncs,nss,
     >        nums,R,SS,nrep)
      call LocateSG(nx,ny,nz,hx,hy,hz,xi,yi,zi,nss,
     >                    SS,nreps)
!         call ptsForceLocate(nx,ny,nz,hx,hy,hz,xi,yi,zi,ncs,nss,
!     >                    nums,R,SS,ptsForce)
      endif
c     
      Call Initial  (nx,ny,nz,nmx,nmy,nmz,npt,
     >     R,x,y,z,DnVX,DnVX0,Uinter,f0,V_X0,V_Y0,V_Z0,PR0)
      do i=1,npt
         DnVX(i)=0.d0 ; DnVX0(i)=1.d0
         DnVY(i)=0.d0 ; DnVY0(i)=1.d0
         DnVZ(i)=0.d0 ; DnVZ0(i)=1.d0
      enddo
      do k=1,nz
         do j=1,ny
            do i=1,nx
               PR0 (k,j,i)=0.d0
               PR1 (k,j,i)=0.d0
            enddo
         enddo
      enddo
c------------------------Iterations en temps------------------------------C
      print*,'================================'
      print*,' starting time  iterations....  '
      print*,'================================'
c     
c     bx(-1)=0.;by(-1)=0.;bz(-1)=0.
c     bx(0)=ax(1);by(0)=ay(1);bz(0)=az(1)
c     do i=1,itermax_tp
c     bx(i)=ax(1);by(i)=ay(1);bz(i)=az(1)
c     enddo
      dtp=dt
      DO i_temps=1,itermax_tp
c     
c     if(i_temps.ge.100) dt=dtp/2
c     ay(1)=0.25d0*dcos(0.5*i_temps*dt)
c     az(1)=0.5d0*dsin(0.5*i_temps*dt)
c     C-------------------------------------------------------------------------C
c     C------------------------Computing the Velocity---------------------------C
c     C-------------------------------------------------------------------------C
c     
c     update U_n (u_n=f qui vient d'etre calculer a l'iteration precedente)
c     terme de convection
c     U_conv=U_n(X-dt*velocity)
         call charact_Vect (nx,ny,nz,nmx,nmy,nmz,xi,yi,zi,xl,yl,zl,
     >        hx,hy,hz,dt,x,y,z,num,V_x0,V_y0,V_z0,UC_x,UC_y,UC_z)
c     
         print*,'-----------------------------------------------'
         print*,' Computing the V_XXXXXX at time', i_temps,'x',dt,' ...'
         print*,'-----------------------------------------------'
C     
         beta=dt*nu
C     
c     calcul des matrices beta*A et M
         call mkdmt3d(nx,ny,nz,hx,hy,hz,bcVX,beta,
     >        a1 ,b1 ,c1 ,d1 ,a2 ,b2 ,c2 ,d2 ,a3 ,b3 ,c3 ,d3 )
         call right_side(.false.,1,nx,ny,nz,nmx,nmy,nmz,np,dt,i_temps,
     >        nu,R,x,y,z,ax,ay,az,rho,c1,d1,c2,d2,c3,d3,uo,g)
C     
c     u_n= la matrice de masse X U_conv
         Call Prod_MU (nx,ny,nz,c1,d1,c2,d2,c3,d3,
     >        uo,nmy,nmz,UC_x,nmy,nmz) 
c     2nd mbre sans C.L
         call Update (nx,ny,nz,nmx,nmy,nmz,dt,1.d0,Uo,g)
         call AssembleDXPv (2,nx,ny,nz,nc,nmx,nmy,nmz,-dt,            
     >        x,y,z,num,bcVX,PR1,g)
C     
C------------------Debut de l'algo du pt fixe-------------------------------C
C     
         print*,' '
         print*, 'starting the  fix point iterations for V_XXXXXX....'
         print*,' '
c     
         if (np.ne.0) then
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     do l=1,np
                        rho(l)=(x(i)-ax(l))**2+(y(j)-ay(l))**2+
     >                       (z(k)-az(l))**2
                        if(rho(l).le.r**2) g(k,j,i)=0.D0
                     enddo
                  enddo
               enddo
            enddo
         endif
C     
         if (dbcX.ne.0) then
c     Dirichlet Boundary conditions
            call Dirich_B_C(.false.,1,nx_fdX,ny_fdX,nz_fdX,nx,ny,nz,
     >           nmx,nmy,nmz,1.d0,i_temps,dt,R,bcVX,a1,b1,c1,d1,a2,b2,
     >           c2,d2,a3,b3,c3,d3,x,y,z,uo,UC_x,g)
         endif
c     
         DO ii=1,itermax        ! point fix iteration on V_XXX
c     
            if (np.ne.0) then
               CALL Normal_Derivative (ncs,np,r,theta,epsilon,beta,
     >              nums,ax,ay,az,Uinter,SS,DnVX)
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
               CALL S_SNDMBRE (np,nss,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,
     >              R,nrep,num,nums,bcVX,ax,ay,az,SS,x,y,z,Area,DnVX,g,
     >              F0)
            else
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
            endif
c     
            call dcq3d(nx_fdX,ny_fdX,nz_fdX,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.true.,ierr)
c     Solve the given problem with the subroutine dcq3d.
            call dcq3d(nx_fdX,ny_fdX,nz_fdX,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.false.,
     &           ierr) 
            if (ierr.ne.0) then
               print *, 'Error no ', ierr, ' in solution'
               stop
            end if
C     
c     maximum pointwise error.
            CALL input_analy(.false.,1,nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     >           x,y,z,R,u)
c     
c     reorganisation de la solution en tenant compte des D_BC
c     
            if (dbcX.ne.0) then
               call reorganization (nx_fdX,ny_fdX,nz_fdX,nx,ny,nz,nmx,
     >              nmy,nmz,bcVX,uo,u,f0)
            endif
c     
            if(np.eq.0) goto 234
            CALL INTERPOL_C(np,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >           nrep,num,nums,ax,ay,az,x,y,z,SS,f0,uinter)
            DO i=1,npt
               DnVdiff(i)=DnvX(i)-DnvX0(i)
               DnVX0(i)   =DnVX(i)
            ENDDO
            normeDnV=dsqrt(prodsca(npt,DnVdiff,DnVdiff)) !/theta
            print*,'*  Norme DiffDnv =',normeDnV
            if(mod(ii,10).eq.0) then
               print*, '  '
               print*, '***********************************************'
               print*,'*',ii,' fix point iteration for XXX-Velocity'
               print*,'* DnV(ncs/2) on sphere number 1 =',DnVX(nshow)
               print*,'* Uinter(ncs/2) on sphere number 1 =',
     >              Uinter(nshow)
               print*,'* Uinter(2ncs+ncs/2) on sphere number 2 =',
     >              Uinter(nshow2)
               print*,'* Norme DiffDnv =',normeDnV
               print*, '***********************************************'
               print*, '  '
            endif
            if(normeDnV.le.prec) goto 234
C     
         ENDDO                  ! point fix iteration on V_XXX
C     
         print*, ' '
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, '!!*           VELOCITY_XXX               *!!'
         PRINT*, '!!* MAXIMAL ITERATION NUMBER IS REACHED  *!!'
         print*, '!!*           VELOCITY_XXX               *!!'
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, ' '
C     
C------------------Fin de l'algo du pt fixe-------------------------------C
C     
 234     print*, '###################################################' 
         print*, '# V_XXX COMPUTING at time',i_temps,'x',dt,'...done' 
         print*, '# EXIT AFTER',ii,' ITERATIONS OF FIX POINT ALGORITHM'
         print*, '# NORM =',normeDnV
         print*, '#Uinter(ncs/2) on sphere number 1 =',Uinter(nshow)
         print*, '#Uinter(2ncs+ncs/2) on sphere number 2 =',
     >        Uinter(nshow2)
         print*, '###################################################'
c     
         call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,f0,V_X0) 
C     
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     DO i=1,npt                  !%%
!     DnVdiff(i)= 0.d0            !%%
!     DnV0(i)   = 1.d0            !%%
!     DnV(i)    = 0.d0            !%%
!     ENDDO                       !%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     
!--------------------VYYYYYYYYYYYYYYYYYYYYYYY-------------------------------
C     
c     calcul des matrices beta*A et M
         call mkdmt3d(nx,ny,nz,hx,hy,hz,bcVY,beta,
     >        a1 ,b1 ,c1 ,d1 ,a2 ,b2 ,c2 ,d2 ,a3 ,b3 ,c3 ,d3 )
         call right_side(.false.,2,nx,ny,nz,nmx,nmy,nmz,np,dt,i_temps,
     >        nu,R,x,y,z,ax,ay,az,rho,c1,d1,c2,d2,c3,d3,uo,g)
C     
c     update U_n (u_n=f qui vient d'etre calculer a l'iteration precedente)
c     terme de convection
c     U_conv=U_n(X-dt*velocity)
c     u_n= la matrice de masse X U_conv
         Call Prod_MU (nx,ny,nz,c1,d1,c2,d2,c3,d3,
     >        uo,nmy,nmz,UC_y,nmy,nmz) 
c     2nd mbre sans C.L
         call Update (nx,ny,nz,nmx,nmy,nmz,dt,1.d0,Uo,g)
         call AssembleDYPv (2,nx,ny,nz,nc,nmx,nmy,nmz,-dt,
     >        x,y,z,num,bcVY,PR1,g)
C     
C------------------Debut de l'algo du pt fixe-------------------------------C
C     
         print*,' '
         print*, 'starting the  fix point iterations for V_YYYYYY....'
         print*,' '
c     
         if (np.ne.0) then
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     do l=1,np
                        rho(l)=(x(i)-ax(l))**2+(y(j)-ay(l))**2+
     >                       (z(k)-az(l))**2
                        if(rho(l).le.r**2) g(k,j,i)=0.D0
                     enddo
                  enddo
               enddo
            enddo
         endif
c     
         if (dbcY.ne.0) then
c     Dirichlet Boundary conditions
            call Dirich_B_C(.false.,2,nx_fdY,ny_fdY,nz_fdY,nx,ny,nz,
     >           nmx,nmy,nmz,1.d0,i_temps,dt,R,bcVY,a1,b1,c1,d1,a2,b2,
     >           c2,d2,a3,b3,c3,d3,x,y,z,uo,UC_y,g)
         endif
c     
         DO ii=1,itermax        ! point fix iteration on V_YYY
c     
            if (np.ne.0) then
               CALL Normal_Derivative (ncs,np,r,theta,epsilon,beta,
     >              nums,ax,ay,az,Uinter,SS,DnVY)
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
               CALL S_SNDMBRE (np,nss,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,
     >              R,nrep,num,nums,bcVY,ax,ay,az,SS,x,y,z,Area,DnVY,g,
     >              F0)
            else
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
            endif
c     
            call dcq3d(nx_fdY,ny_fdY,nz_fdY,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.true.,ierr)
c     Solve the given problem with the subroutine dcq3d.
            call dcq3d(nx_fdY,ny_fdY,nz_fdY,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.false.,
     &           ierr)
            if (ierr.ne.0) then
               print *, 'Error no ', ierr, ' in solution'
               stop
            end if
C     
c     maximum pointwise error.
            CALL input_analy(.false.,2,nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     >           x,y,z,R,u)
c     
c     reorganisation de la solution en tenant compte des D_BC
c     
            if (dbcY.ne.0) then
               call reorganization (nx_fdY,ny_fdY,nz_fdY,nx,ny,nz,nmx,
     >              nmy,nmz,bcVY,uo,u,f0)
            endif
c     
            if(np.eq.0) goto 235
            CALL INTERPOL_C(np,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >           nrep,num,nums,ax,ay,az,x,y,z,SS,f0,uinter)
            DO i=1,npt
               DnVdiff(i)=DnvY(i)-DnvY0(i)
               DnVY0(i)   =DnVY(i)
            ENDDO
            normeDnV=dsqrt(prodsca(npt,DnVdiff,DnVdiff)) !/theta
            print*,'*  Norme DiffDnv =',normeDnV
            if(mod(ii,10).eq.0) then
               print*, '  '
               print*, '***********************************************'
               print*,'*',ii,' fix point iteration for YYY-Velocity'
               print*,'* DnV(ncs/2) on sphere number 1 =',DnVY(nshow)
               print*,'* Uinter(ncs/2) on sphere number 1 =',
     >              Uinter(nshow)
               print*,'* Uinter(2ncs+ncs/2) on sphere number 2 =',
     >              Uinter(nshow2)
               print*,'* Norme DiffDnv =',normeDnV
               print*, '***********************************************'
               print*, '  '
            endif
            if(normeDnV.le.prec) goto 235
c     
         ENDDO                  ! point fix iteration on V_YYY
c     
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, '!!*           VELOCITY_YYY               *!!'
         PRINT*, '!!* MAXIMAL ITERATION NUMBER IS REACHED  *!!'
         print*, '!!*           VELOCITY_YYY               *!!'
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, ' '
C     
C------------------Fin de l'algo du pt fixe-------------------------------C
C     
 235     print*, '###################################################' 
         print*, '# V_YYY COMPUTING at time',i_temps,'x',dt,'...done' 
         print*, '# EXIT AFTER',ii,' ITERATIONS OF FIX POINT ALGORITHM'
         print*, '# NORM =',normeDnV
         print*, '#Uinter(ncs/2) on sphere number 1 =',Uinter(nshow)
         print*, '#Uinter(2ncs+ncs/2) on sphere number 2 =',
     >        Uinter(nshow2)
         print*, '###################################################'
c     
         call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,f0,V_Y0) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     DO i=1,npt                  !%%
!     DnVdiff(i)= 0.d0            !%%
!     DnV0(i)   = 1.d0            !%%
!     DnV(i)    = 0.d0            !%%
!     ENDDO                       !%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
!-----------------------VZZZZZZZZZZZZZZZZZZZZ-------------------------------
C     
c     calcul des matrices beta*A et M
         call mkdmt3d(nx,ny,nz,hx,hy,hz,bcVZ,beta,
     >        a1 ,b1 ,c1 ,d1 ,a2 ,b2 ,c2 ,d2 ,a3 ,b3 ,c3 ,d3 )
         call right_side(.false.,3,nx,ny,nz,nmx,nmy,nmz,np,dt,i_temps,
     >        nu,R,x,y,z,ax,ay,az,rho,c1,d1,c2,d2,c3,d3,uo,g)
C     
c     update U_n (u_n=f qui vient d'etre calculer a l'iteration precedente)
c     terme de convection
c     U_conv=U_n(X-dt*velocity)
c     u_n= la matrice de masse X U_conv
         Call Prod_MU (nx,ny,nz,c1,d1,c2,d2,c3,d3,
     >        uo,nmy,nmz,UC_z,nmy,nmz) 
c     2nd mbre sans C.L
         call Update (nx,ny,nz,nmx,nmy,nmz,dt,1.d0,Uo,g)
         call AssembleDZPv (2,nx,ny,nz,nc,nmx,nmy,nmz,-dt,
     >        x,y,z,num,bcVZ,PR1,g)
C     
         print*,' '
         print*, 'starting the  fix point iterations for V_ZZZZZZ....'
         print*,' '
c     
         if (np.ne.0) then
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     do l=1,np
                        rho(l)=(x(i)-ax(l))**2+(y(j)-ay(l))**2+
     >                       (z(k)-az(l))**2
                        if(rho(l).le.r**2) g(k,j,i)=0.D0
                     enddo
                  enddo
               enddo
            enddo
         endif
C     
         if (dbcZ.ne.0) then
c     Dirichlet Boundary conditions
            call Dirich_B_C(.false.,3,nx_fdZ,ny_fdZ,nz_fdZ,nx,ny,nz,
     >           nmx,nmy,nmz,1.d0,i_temps,dt,R,bcVZ,a1,b1,c1,d1,a2,b2,
     >           c2,d2,a3,b3,c3,d3,x,y,z,uo,UC_z,g)
         endif
C     
C------------------Debut de l'algo du pt fixe-------------------------------C
C     
         
c     
         DO ii=1,itermax        ! point fix iteration on V_ZZZ
C     
            if (np.ne.0) then
               CALL Normal_Derivative (ncs,np,r,theta,epsilon,beta,
     >              nums,ax,ay,az,Uinter,SS,DnVZ)
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
               CALL S_SNDMBRE (np,nss,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,
     >              R,nrep,num,nums,bcVZ,ax,ay,az,SS,x,y,z,Area,DnVZ,g,
     >              F0)
            else
               call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,g,f0)
            endif
c     
            call dcq3d(nx_fdZ,ny_fdZ,nz_fdZ,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.true.,ierr)
c     Solve the given problem with the subroutine dcq3d.
            call dcq3d(nx_fdZ,ny_fdZ,nz_fdZ,f0,nmy,nmz,a1,b1,c1,d1,
     &           a2,b2,c2,d2,a3,b3,c3,d3,1.d0,dw,ldw,iw,liw,.false.,
     &           ierr)
            if (ierr.ne.0) then
               print *, 'Error no ', ierr, ' in solution'
               stop
            end if
C     
c     maximum pointwise error.
            CALL input_analy(.false.,3,nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     >           x,y,z,R,u)
c     
c     reorganisation de la solution en tenant compte des D_BC
c     
            if (dbcZ.ne.0) then
               call reorganization (nx_fdZ,ny_fdZ,nz_fdZ,nx,ny,nz,nmx,
     >              nmy,nmz,bcVZ,uo,u,f0)
            endif
c     
            if(np.eq.0) goto 236
            CALL INTERPOL_C(np,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >           nrep,num,nums,ax,ay,az,x,y,z,SS,f0,uinter)
            DO i=1,npt
               DnVdiff(i)=DnvZ(i)-DnvZ0(i)
               DnVZ0(i)   =DnVZ(i)
            ENDDO
            normeDnV=dsqrt(prodsca(npt,DnVdiff,DnVdiff)) !/theta
            print*,'*  Norme DiffDnv =',normeDnV
            if(mod(ii,10).eq.0) then
               print*, '  '
               print*, '***********************************************'
               print*,'*',ii,' fix point iteration for ZZZ-Velocity'
               print*,'* DnV(ncs/2) on sphere number 1 =',DnVZ(nshow)
               print*,'* Uinter(ncs/2) on sphere number 1 =',
     >              Uinter(nshow)
               print*,'* Uinter(2ncs+ncs/2) on sphere number 2 =',
     >              Uinter(nshow2)
               print*,'* Norme DiffDnv =',normeDnV
               print*, '***********************************************'
               print*, '  '       
            endif
            if(normeDnV.le.prec) goto 236
C     
         ENDDO                  ! point fix iteration on V_ZZZ
C     
         print*, ' '
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, '!!*           VELOCITY_ZZZ               *!!'
         PRINT*, '!!* MAXIMAL ITERATION NUMBER IS REACHED  *!!'
         print*, '!!*           VELOCITY_ZZZ               *!!'
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, ' '
C     
C------------------Fin de l'algo du pt fixe-------------------------------C
C     
 236     print*, '###################################################' 
         print*, '# V_ZZZ COMPUTING at time',i_temps,'x',dt,'...done' 
         print*, '# EXIT AFTER',ii,' ITERATIONS OF FIX POINT ALGORITHM'
         print*, '# NORM =',normeDnV
         print*, '#Uinter(ncs/2) on sphere number 1 =',Uinter(nshow)
         print*, '#Uinter(2ncs+ncs/2) on sphere number 2 =',
     >        Uinter(nshow2)
         print*, '###################################################'
C     
         call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,f0,V_Z0) 

c     calcul debit
         ! toto debit ??????????????
         call DebitZ (2,nx,ny,nz,nc,nmx,nmy,nmz,hx,hy,hz,
     >      num,x,y,z,V_Z0,debitZF)
         write(33,*) debitZF
         write(*,*) '---->Debit ', i_temps, ' ', debitZF

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     DO i=1,npt                  !%%
c     DnVdiff(i)= 0.d0            !%%
c     DnV0(i)   = 1.d0            !%%
c     DnV(i)    = 0.d0            !%%
c     ENDDO                       !%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     
c     C-------------------------------------------------------------------------C
c     C------------------------Computing the PRESSURE---------------------------C
c     C-------------------------------------------------------------------------C
C     
         print*,'----------------------------------------------'
         print*,' Computing the Pressure at time', i_temps,'x',dt,' ...'
         print*,'----------------------------------------------'
c     calcul des matrices 1*A et M
         call mkdmt3d(px,py,pz,hpx,hpy,hpz,bcP,1.D0,
     >        a1 ,b1 ,c1 ,d1 ,a2 ,b2 ,c2 ,d2 ,a3 ,b3 ,c3 ,d3 )
         call right_side(.true.,1,px,py,pz,nmx,nmy,nmz,np,dt,i_temps,nu,
     >        R,xp,yp,zp,ax,ay,az,rho,c1,d1,c2,d2,c3,d3,uo,g)
c     
c     do k=1,nz
c     do j=1,ny
c     do i=1,nx
c     V_X0(k,j,i)=0.d0
c     V_Y0(k,j,i)=0.d0
c     V_Z0(k,j,i)=0.25d0*(x(i)**2-1.d0)*(y(j)**2-1.d0)
c     enddo
c     enddo
c     enddo
         call  AssembleDivUq (switch,2,px,py,pz,pc,nmx,nmy,nmz,-1.d0/dt,
     >        xp,yp,zp,x,y,z,numP,bcP,V_X0,V_Y0,V_Z0,g) 
c     call test_quad (2,px,py,pz,pc,nmx,nmy,nmz,
c     >                ax,ay,az,x,y,z,R,num,bcP,PR1,g)
C     
         if(np.ne.0) then
            CALL Cubes_Num_B (px,py,pz,pc,np,NumP,xi,yi,zi,hpx,hpy,hpz,
     >           R,ax,ay,az,ncb,numb)
            do k=1,pz
               do j=1,py
                  do i=1,px
                     do l=1,np
                        rho(l)=(xp(i)-ax(l))**2+(yp(j)-ay(l))**2+
     >                       (zp(k)-az(l))**2
                        if(rho(l).le.r**2) g(k,j,i)=0.D0
                     enddo
                  enddo
               enddo
            enddo
         endif
c     
         if (dbcP.ne.0) then
c     Dirichlet Boundary conditions
            call Dirich_B_C(.true.,1,nx_fdP,ny_fdP,nz_fdP,px,py,pz,
     >           nmx,nmy,nmz,zeta,i_temps,dt,R,bcP,a1,b1,c1,d1,a2,b2,c2,
     >           d2,a3,b3,c3,d3,x,y,z,uo,u1,g)
         endif
c     
         DO ii=1,itermaxP       ! point fix iteration on pressure
C     
            call Update (px,py,pz,nmx,nmy,nmz,0.d0,1.d0,g,f)
c     
            if(np.ne.0)then
               CALL SND_B_Neu (px,py,pz,pc,nmx,nmy,nmz,nr,np,ncb,R,
     >              hpx,hpy,hpz,ax,ay,az,xp,yp,zp,numb,bcP,uo,f)
            endif
c     
c     Solve the given problem with the subroutine dc3d.
            call dcq3d(nx_fdP,ny_fdP,nz_fdP,f,nmy,nmz,a1,b1,c1,d1,a2,b2,
     >           c2,d2,a3,b3,c3,d3,zeta,dwP,ldwP,iwP,liwP,.true.,ierr)
            call dcq3d(nx_fdP,ny_fdP,nz_fdP,f,nmy,nmz,a1,b1,c1,d1,a2,b2,
     >           c2,d2,a3,b3,c3,d3,zeta,dwP,ldwP,iwP,liwP,.false.,ierr)
            if (ierr.ne.0) then
               print *, 'Error no ', ierr, ' in solution'
               stop
            end if
c     
c     maximum pointwise error.
            CALL input_analy(.true.,1,px,py,pz,nmx,nmy,nmz,i_temps,dt,
     >           xp,yp,zp,R,u)
c     
c     reorganisation de la solution en tenant compte des D_BC
c     
            if (dbcP.ne.0) then
               call reorganization (nx_fdP,ny_fdP,nz_fdP,px,py,pz,nmx,
     >              nmy,nmz,bcP,u1,u,f)
            endif
c     
            if (np.eq.0) goto 244
c     
            do k=1,pz
               do j=1,py
                  dO i=1,px
                     u1(k,j,i)=dabs(f(k,j,i)-uo(k,j,i))
                  enddo
               enddo
            enddo
            normeDnV=dsqrt(prodsca_3d(nmx,nmy,nmz,px,py,pz,u1,u1))
c     
            print*,'*  Norm Diff =',normeDnV
c     
            call Update (px,py,pz,nmx,nmy,nmz,0.d0,1.d0,f,uo)
c     
            if(mod(ii,10).eq.0) then
               print*, '  '
               print*, '***********************************************'
               print*,'*',ii,' fix point iteration for --PRESSURE--'
               print*,'*  u(n/2)  =',f(nshowv,nshowv,nshowv)
               print*,'*  Norm Diff =',normeDnV
               print*, '***********************************************'
               print*, '  '
            endif
            if(normeDnV.le.prec) goto 244
c     
         ENDDO                  ! point fix iteration on pressure
c     
         print*, ' '
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, '!!*            PRESSURE                  *!!'
         PRINT*, '!!* MAXIMAL ITERATION NUMBER IS REACHED  *!!'
         print*, '!!*            PRESSURE                  *!!'
         print*, '!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!*!!'
         print*, ' '
C----------------FIN des Iterations de l'algo du point fixe---------C
C     
 244     print*, '###################################################' 
         print*, '# PRESSURE COMPUTING at time',i_temps,'x',dt,
     >        ' ...done' 
c         '
         print*, '# EXIT AFTER',ii,' ITERATIONS OF FIX POINT ALGORITHM'
         print*, '# NORM =',normeDnV
         print*, '###################################################'
c     
         if(switch) then
            call PROJECT_GF (nx,ny,nz,px,py,pz,nmx,nmy,nmz,numP,
     >           x,y,z,xp,yp,zp,f,PRF)
         else
            call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,f,PRF)
         endif
c     
         call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,PR1,PR0)
         call Update (nx,ny,nz,nmx,nmy,nmz,0.d0,1.d0,PRF,PR1)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  if(i.ne.1.and.i.ne.nx) then
                     UC_x(k,j,i)=V_X0(k,j,i)-dt*
     >                    ((PR1(k,j,i+1)-PR1(k,j,i))/hx)
                  else
                     UC_x(k,j,i)=V_X0(k,j,i)
                  endif
                  if(j.ne.1.and.j.ne.ny) then
                     UC_y(k,j,i)=V_Y0(k,j,i)-dt*
     >                    ((PR1(k,j+1,i)-PR1(k,j,i))/hy)
                  else
                     UC_y(k,j,i)=V_Y0(k,j,i)
                  endif
                  if(k.ne.1.and.k.ne.nz) then
                     UC_z(k,j,i)=V_Z0(k,j,i)-dt*
     >                    ((PR1(k+1,j,i)-PR1(k,j,i))/hz)
                  else
                     UC_z(k,j,i)=V_Z0(k,j,i)
                  endif
                  g(k,j,i)=dsqrt(UC_X(k,j,i)**2+UC_Y(k,j,i)**2+
     >                 UC_Z(k,j,i)**2)
               end do
            end do
         end do
c%%%%%%%%%%%%
         WRITE(35,*)'$DATA = CONTOUR'
         WRITE(35,*)'%meshplot = off'
         WRITE(35,*)'%contstyle=2'
         WRITE(35,*)'%xmin =', xi,' xmax=', xi+xl
         WRITE(35,*)'%ymin =', zi,' ymax=', zi+zl
         WRITE(35,*)'%nx =  ',nx
         WRITE(35,*)'%ny =  ',nz
         WRITE(35,*)'%nsteps =',50
         
         do k=1,nz
            do j=ny/2,ny/2
               do i=1,nx
                  write(35,*) PR1(k,j,i)
               enddo
            enddo
         enddo
c%%%%%%%%%%%%
         do kc=1,nz
            do jc=1,ny
               do ic=1,nx
                  do l=1,np
                     rho(l)=(x(ic)-ax(l))**2+(y(jc)-ay(l))**2+
     >                    (z(kc)-az(l))**2
                     if(rho(l).le.(r+hx)**2) then 
                        UC_x(kc,jc,ic)=0.d0
                        UC_y(kc,jc,ic)=0.d0
                        UC_z(kc,jc,ic)=0.d0
                        g(kc,jc,ic)=0.d0
                     endif
                  enddo
               enddo
            enddo
         enddo
c%%%%%%%%%%%%
         call Visu_mtv (nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     >        xi,yi,zi,xl,yl,zl,g,U)
c     call View_gradient (nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
c     > xi,yi,zi,xl,yl,zl,x,y,z,UC_x,UC_y,UC_z)
         if (i_temps.eq.1.or.mod(i_temps,2).eq.0) then
!            call Vof (nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,
!     >           vfn,cfn,x,y,z,UC_x,UC_Z)
c      call Visu_DX (0,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,vfn,
c     >            x,y,z,g)
!            call Visu_DX (1,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,xfn,
!     >           x,y,z,UC_X)
!            call Visu_DX (2,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,yfn,
!     >           x,y,z,UC_Y)
!            call Visu_DX (3,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,zfn,
!     >           x,y,z,UC_Z)
!            call Visu_DX (4,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,pfn,
!     >           x,y,z,PR1)
         endif  
c      call Visu_DX (0,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,pfn,
c     >            x,y,z,PR1)
C     
c     write(38,*)i_temps, Integ_G(ncs, area, DnVZ)
c     ax(1) = ax(1)
c     ay(1)=0.25d0*dcos(0.5*i_temps*dt)
c     bx(i_temps)=nu*(2.*bx(i_temps-1)-bx(i_temps-2)-
c     >             0.5*dt*Integ_G(ncs, area, DnVX) )
c     by(i_temps)=nu*(2.*by(i_temps-1)-by(i_temps-2)-
c     >             0.5*dt*Integ_G(ncs, area, DnVY) )
c     bz(i_temps)=nu*(2.*bz(i_temps-1)-bz(i_temps-2)-
c     >             0.5*dt*Integ_G(ncs, area, DnVZ) )
c     ax(1)=bx(i_temps);ay(1)=by(i_temps);az(1)=bz(i_temps)
c     
c     write(36,*)i_temps,ax(1)
c     write(37,*)i_temps,ay(1)
c     write(38,*)i_temps,az(1)
         if (np.ne.0) then
            CALL INTERPOL_C(np,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >           nrep,num,nums,ax,ay,az,x,y,z,SS,PR1,Pinter)
            call Integ_GP(ncs, area, nums, SS, R, Pinter,
     >           Integ_GPX,Integ_GPY,Integ_GPZ)
            
            write(36,*)i_temps,
     >           -Integ_GPX , nu*Integ_G(ncs, area, DnVX)
            write(37,*)i_temps,
     >           -Integ_GPY , nu*Integ_G(ncs, area, DnVY)
            write(38,*)i_temps,
     >           -Integ_GPZ , nu*Integ_G(ncs, area, DnVZ)
            write(39,*)i_temps,dsqrt(
     >        (-Integ_GPX + nu*Integ_G(ncs, area, DnVX))**2 +
     >        (-Integ_GPY + nu*Integ_G(ncs, area, DnVY))**2 +    
     >        (-Integ_GPZ + nu*Integ_G(ncs, area, DnVZ))**2
     >                              )
          endif
      END DO                    ! Time iteration
!%%%%%%%%
      print*, 2.d0*r/nu,dsqrt(
     >        (-Integ_GPX + nu*Integ_G(ncs, area, DnVX))**2 +
     >        (-Integ_GPY + nu*Integ_G(ncs, area, DnVY))**2 +    
     >        (-Integ_GPZ + nu*Integ_G(ncs, area, DnVZ))**2
     >                              )
      print*, Integ_GPX,Integ_GPY,Integ_GPZ
!%%%%%%%%
      cpu = Time_Cpu() - cpu
c     
c     
      print*, ' '
      print*, '#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#'
      print*, '#*#*#*#*#*#*#*#*#*#*#TOTAL COMPUTING...DONE#*#*#*#*#*#*#'
      print*, '#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#**#*#*#*#*#*#*#'
      print*, ' '
      print*, '##########################################'
      print*, '# Total cpu time =',cpu,'CPU seconds ' 
      print*, '##########################################'
      print*, ' '
      print*,'Prepare visualization''s files....'
      print*, ' '
c     
c      if (np.ne.0) call view_sphere (np,ncs,SS,ax,ay,az,nums)
!??????????
            call INTERPOL_S(np,ncs,nss,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >             nrepS,num,nums,ax,ay,az,x,y,z,SS,pr1,PinterS)
            call sphere4D (nss,ncs,SS,ax,ay,az,nums,PinterS)
!??????????
C-------------------fin des iterations en temps------------------------------C
C     
      close(36)
      close(37)
      close(38)
      stop
      end
      
