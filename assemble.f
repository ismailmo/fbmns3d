
c Author : Mourad Ismail   (ismail@ujf-grenoble.fr)

C--------------------------------------------------------------------C
C--------------------------ROUTINES----------------------------------C
C--------------------------------------------------------------------C
C
C========================================================================C
      subroutine mkdmt3d(n1,n2,n3,h1,h2,h3,bc,beta,
     >                   a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3)
C========================================================================C
      integer n1, n2, n3
      double precision a1(n1), b1(n1), c1(n1), d1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision h1, h2, h3, beta
c
      integer i , bc(6)
c
      if (beta.lt.0.d0) then
       print*, 'Error in Matrix assemble. Beta must be positive...ABORT'
       stop
      endif
c
      b1(1) = (1.d0+bc(1))*(1.d0/h1)*beta
      a1(1) =  0.d0
      d1(1) = (1.d0+bc(1))*(h1/3.d0)!*beta
      c1(1) =  0.d0
      do i=2,n1-1
         b1(i) =  (2.d0/h1)*beta
         a1(i) = (-1.d0/h1)*beta
         d1(i) = ((2.d0*h1)/3.d0)!*beta
         c1(i) =  (h1/6.d0)!*beta
      end do
      b1(n1) = (1.d0+bc(2))*(1.d0/h1)*beta
      a1(n1) = (-1.d0/h1)*beta
      d1(n1) = (1.d0+bc(2))*(h1/3.d0)!*beta
      c1(n1) =  (h1/6.d0)!*beta
c
      b2(1) =  (1.d0+bc(3))*(1.d0/h2)*beta
      a2(1) =  0.d0
      d2(1) =  (1.d0+bc(3))*(h2/3.d0)!*beta
      c2(1) =  0.d0
      do i=2,n2-1
         b2(i) =  (2.d0/h2)*beta
         a2(i) = (-1.d0/h2)*beta
         d2(i) = ((2.d0*h2)/3.d0)!*beta
         c2(i) =  (h2/6.d0)!*beta
      end do
      b2(n2) =  (1.d0+bc(4))*(1.d0/h2)*beta
      a2(n2) =  (-1.d0/h2)*beta
      d2(n2) =  (1.d0+bc(4))*(h2/3.d0)!*beta
      c2(n2) =  ( h2/6.d0)!*beta

c
      b3(1) =  (1.d0+bc(5))*(1.d0/h3)*beta
      a3(1) =   0.d0
      d3(1) =  (1.d0+bc(5))*(h3/3.d0)!*beta
      c3(1) =   0.d0
      do i=2,n3-1
         b3(i) =  (2.d0/h3)*beta
         a3(i) = (-1.d0/h3)*beta
         d3(i) = ((2.d0*h3)/3.d0)!*beta
         c3(i) =  (h3/6.d0)!*beta
      end do
      b3(n3) =  (1.d0+bc(6))*(1.d0/h3)*beta
      a3(n3) =  (-1.d0/h3)*beta
      d3(n3) =  (1.d0+bc(6))*(h3/3.d0)!*beta
      c3(n3) =  ( h3/6.d0)!*beta
c
      return
      end
c
C=======================================================================C
      subroutine right_side(switch,l,nx,ny,nz,nmx,nmy,nmz,np,dt,i_temps,
     >      nu,R,x,y,z,ax,ay,az,rho,c1,d1,c2,d2,c3,d3,u,f)
C=======================================================================C
      integer nmx,nmy,nmz,nx,ny,nz,np,i_temps,l
      double precision dt,nu,R,x(nx),y(ny),z(nz)
      double precision ax(*),ay(*),az(*),rho(*)
      double precision c1(nx),d1(nx),c2(ny),d2(ny),c3(nz),d3(nz)
      double precision f(nmz,nmy,nmx),u(nmz,nmy,nmx)
      logical switch
c
      print*,'Assemble the right hand side....'
      !second membre (brut !!!)
      call input_2nd_mbre(switch,l,nx,ny,nz,nmx,nmy,nmz,np,i_temps,dt,
     > nu,x,y,z,ax,ay,az,R,rho,u) 
      !Assemblage du second membre
      !Make the sample right-hand side vector (X par la matrice de masse)
      !f= la matrice de masse X U
      Call Prod_MU (nx,ny,nz,c1,d1,c2,d2,c3,d3,
     >                f,nmy,nmz,u,nmy,nmz) 
      print*, 'done'
      return
      end
c
C======================================================================C
      Subroutine right_side_quad (n,nx,ny,nz,nc,nmx,nmy,nmz,r,hx,hy,hz,
     >      num,x,y,z,F)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,nx,ny,nz,nc,i,j,k,l,ks,n
      integer ix(8),iy(8),iz(8),num(*)
      double precision evalqi,rho,beta2,hx,hy,hz,integ_QN
      double precision pi,r,r2,beta,xmax,xmin,ymin,ymax,zmin,zmax
      double precision rt,rt2,w0,w1,w2,dx,dy,dz,a(5),b(5),c(5)
      double precision f(nmz,nmy,nmx),ga(5,5,5),x(*),y(*),z(*),g(5,5,5)
c
      print*,'Assembling the right hand side with Gauss quad formula'
      print*,'...'
c
      pi=dacos(-1.d0) ; R2=R**2
      beta=2.d0*pi ; beta2=beta**2
c
      do k=1,nz
       do j=1,ny
        do i=1,nx
         F(k,j,i)=0.d0
        enddo
       enddo
      enddo
c
      do k=1,nc
c
       call Coord (num(k),nx,ny,ix,iy,iz)
c
      xmin=x(ix(1)) ; xmax=x(ix(7))
      ymin=y(iy(1)) ; ymax=y(iy(7))
      zmin=z(iz(1)) ; zmax=z(iz(7))
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
       call Input_func(2,n,nmx,nmy,nmz,r,x,y,z,ix,iy,iz,a,b,c,
     >                      f,g)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            ga(i,j,l)=g(i,j,l)*evalqi(ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)
          enddo
        enddo
       enddo
c
      F(iz(ks),iy(ks),ix(ks))=F(iz(ks),iy(ks),ix(ks))+
     >                         integ_QN(n,dx,dy,dz,w0,w1,w2,ga)
      enddo
      enddo
c
      print*,' '
      print*,'done'
      print*,' '
c
      return
      end
!=======================================================================!
      subroutine Prod_MU(n1,n2,n3,c1,d1,c2,d2,c3,d3,
     &                f,ldf2,ldf3,u,ldu2,ldu3)
!=======================================================================!
      integer n1, n2, n3, ldf2, ldf3, ldu2, ldu3
      double precision c1(n1), d1(n1)
      double precision c2(n2), d2(n2)
      double precision c3(n3), d3(n3)
      double precision f(ldf3,ldf2,n1), u(ldu3,ldu2,n1)
c
      integer i, j, k
c
      do i=1,n1
         do j=1,n2
            do k=1,n3
               f(k,j,i) = 0.d0
            end do
         end do
      end do
c
c i=1
c
      if (n1.gt.0) then
         call tmatxv(n2,n3,c2,d2,c3,d3,d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
      end if
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
      end if
c
      do i=2,n1-1
         call tmatxv(n2,n3,c2,d2,c3,d3,c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
      end do
c
c i=n1
c
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
      end if
c
      return
      end
C=======================================================================C
      subroutine Dirich_B_C(switch,l,nx_fd,ny_fd,nz_fd,nx,ny,nz,
     > nmx,nmy,nmz,alpha,i_temps,
     > dt,R,bc,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,x,y,z,uo,u1,f)
C=======================================================================C
C
      implicit none
      integer nx_fd,ny_fd,nz_fd,nx,ny,nz,i_temps,i,j,k,l,bc(6)
      integer nmx,nmy,nmz
      double precision alpha,dt,R
      double precision a1(*),b1(*),c1(*),d1(*),a2(*),b2(*),c2(*),d2(*)
      double precision a3(*),b3(*),c3(*),d3(*),x(*),y(*),z(*)
      double precision f(nmz,nmy,nmx),uo(nmz,nmy,nmx),u1(nmz,nmy,nmx)
      logical switch
C
      print*,'Modify the right side for Dirichlet Boundary condit...'
      CALL input_analy(switch,l,nx,ny,nz,nmx,nmy,nmz,i_temps,
     > dt,x,y,z,R,u1)
C
      do k=1,nz
       do j=1,ny
         do i=1,nx
          uo(k,j,i)=0.d0
         enddo
       enddo
      enddo
C
      if (bc(1).eq.1) then
       do k=1,nz
        do j=1,ny
          uo(k,j,1)=u1(k,j,1)
        enddo
       enddo
      endif
C
      if (bc(2).eq.1) then
       do k=1,nz
        do j=1,ny
          uo(k,j,nx)=u1(k,j,nx)
        enddo
       enddo
      endif
C
      if (bc(3).eq.1) then
       do k=1,nz
        do i=1,nx
          uo(k,1,i)= u1(k,1,i)
        enddo
       enddo
      endif
C
      if (bc(4).eq.1) then
       do k=1,nz
        do i=1,nx
          uo(k,ny,i)=u1(k,ny,i)
        enddo
       enddo
      endif
C
      if (bc(5).eq.1) then
       do j=1,ny
        do i=1,nx
          uo(1,j,i)=u1(1,j,i)
        enddo
       enddo
      endif
C
      if (bc(6).eq.1) then
       do j=1,ny
        do i=1,nx
          uo(nz,j,i)=u1(nz,j,i)
        enddo
       enddo
      endif
C
      !modification du 2nd mbre f en tenant compte de u (non hom D_BC)
      !produit de uo par A+alpha*M:!!! matrices Neumann----------------
      call smxv(nx,ny,nz,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,
     >                alpha,u1,nmy,nmz,uo,nmy,nmz)
c
      do k=1,nz
         do j=1,ny
            do i=1,nx
             f(k,j,i)=f(k,j,i)-u1(k,j,i)
            enddo
         enddo
      enddo
      !elimination des points correpondants a D_BC
       Call Elimination (nx_fd,ny_fd,nz_fd,nmx,nmy,nmz,
     >                   b1,d1,b2,d2,b3,d3,bc,f)
C
      print*, 'done'
      return
      end
*=======================================================================*      
      SUBROUTINE Normal_Derivative (ncs,np,r,theta,epsilon,beta,
     >                          nums,ax,ay,az,Uinter,SS,DnV)
c=========================================================================
C
      IMPLICIT NONE
      INTEGER k,j,kp,np,p,ncs,nums(8,ncs)
      DOUBLE PRECISION r,epsilon,theta,rho,R2,beta
      DOUBLE PRECISION pi,rho2,rhoR2
      DOUBLE PRECISION Uinter(*),DnV(*),SS(3,*),ax(*),ay(*),az(*)
      R2=r**2 ; pi=dacos(-1.d0) 
C
!Le tableau DnV contient les valeurs de la derivee normale de v calculees aux centres
!des elements de gamma de chaque particule. Il est de taille npt=np*ncs.
!pour une particule p, ces valeurs sont sochees dans la partie correspondant a 
! i=(p-1)ncs+1,......,p*ncs
C
      do p=1,np
       DO k=1,ncs
        kp=(  p-1)*ncs+k
         j=(2*p-1)*ncs+k
        DnV(kp)=-((1.D0-theta)/(epsilon))*(Uinter(j))*beta*
     >      (1.D0+epsilon/r)+theta*DnV(kp)
!        Dnv(kp)=-2.d0*beta*r
       enddo
      ENDDO 
      RETURN
      END
c=========================================================================
	SUBROUTINE S_SNDMBRE (np,nsse,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >     nrep,num,nums,bc,ax,ay,az,SS,x,y,z,AireS,DnV,g,F)
*=======================================================================*	
      IMPLICIT NONE
      INTEGER nsse,ncs,i,k,kj,i2,i3,i4,nx,ny,p,trans1,np,trans,
     >        ix(8),iy(8),iz(8),ixfd,iyfd,izfd,nz,nmx,nmy,nmz
      INTEGER num(*),nums(8,*),nrep(*),knrep,bc(6)
      DOUBLEPRECISION evalqi,hx,hy,hz,x1,y1,z1,R,
     > x3,y3,z3,xc,yc,zc,pi,rho
      DOUBLEPRECISION F(nmz,nmy,nmx),DnV(*),SS(3,*),x(*),y(*),z(*),
     >                Aires(*),ax(*),ay(*),az(*),g(nmz,nmy,nmx)
c
      pi=dacos(-1.d0)
c
      DO p=1,np
       trans =(p-1)*ncs
       trans1=int(ax(p)/hx)+int(ay(p)/hy)*(nx-1)+
     >        int(az(p)/hz)*(nx-1)*(ny-1)
	DO k=1,ncs
	 i   =nums(1,k) 
	 i2  =nums(2,k) 
	 i3  =nums(3,k) 
	 i4  =nums(4,k) 
c
	 knrep=num(nrep(k)+trans1)
c
	 x1=SS(1,i ) ; y1=SS(2,i ) ; z1=SS(3,i )
	 x3=SS(1,i3) ; y3=SS(2,i3) ; z3=SS(3,i3)
	 xc=0.5D0*(x1+x3)  ; yc=0.5D0*(y1+y3)  ; zc=0.5D0*(z1+z3)
         rho=dsqrt(xc**2+yc**2+zc**2)
         xc=ax(p)+xc*R/rho ; yc=ay(p)+yc*R/rho ; zc=az(p)+zc*R/rho
c
        call Coord (knrep,nx,ny,ix,iy,iz)
c
	DO kj=1,8
         ixfd=ix(kj)-bc(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(kj)-bc(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(kj)-bc(5) !decalage du a l'elimination des sommets Dirichlet 
	F(izfd,iyfd,ixfd)=
     >        F(izfd,iyfd,ixfd) + (
     >                evalqi(kj,ix,iy,iz,xc,yc,zc,x,y,z)*
     >                DnV(k+trans))*Aires(k)
	ENDDO
	ENDDO
      ENDDO
C
	RETURN
	END
!=======================================================================!
      subroutine AUo(n1,n2,n3,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,
     &                ch,f,ldf2,ldf3,u,ldu2,ldu3)
!=======================================================================!
      integer n1, n2, n3, ldf2, ldf3, ldu2, ldu3
      double precision a1(n1), b1(n1), c1(n1), d1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision f(ldf3,ldf2,n1), u(ldu3,ldu2,n1), ch
c
      integer i, j, k
c
      do i=1,n1
         do j=1,n2
            do k=1,n3
               f(k,j,i) = 0.d0
            end do
         end do
      end do
c
c i=1
c
      if (n1.gt.0) then
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
      end if
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
      end if
c
      do i=2,n1-1
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
      end do
c
c i=n1
c
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
      end if
c
      return
      end
C=======================================================================C      
	SUBROUTINE SND_B_Neu (nx,ny,nz,nc,nmx,nmy,nmz,nr,np,ncb,R,
     >          hx,hy,hz,ax,ay,az,x,y,z,numb,bcV,u,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER ncb(*),j,k,ik1,ixfd,iyfd,izfd,ix(8),iy(8),iz(8),nc,
     >   np,p,numb(nc,np),nx,ny,nz,nmx,nmy,nmz,bcV(6)
	DOUBLEPRECISION hx,hy,hz,R,F(nmz,nmy,nmx),u(nmz,nmy,nmx),
     >           x(*),y(*),z(*),R2,ax(*),ay(*),az(*)
        integer nr,nr2,nr23,nr21,ir,jr,kr
        doubleprecision hxr,hyr,hzr,vol,R2D,R2D1,
     >  xc,yc,zc,interpol_Q1_X,
     >  interpol_Q1_Y,interpol_Q1_Z,
     >  eval_GX_qi,eval_GY_qi,eval_GZ_qi,rho_raf
C
      nr2=2**nr   ; nr21=2*nr2  
      hxr=hx/nr21 ; hyr=hy/nr21 ; hzr=hz/nr21 
      nr23=nr2**3 ; vol=(hx*hy*hz)/nr23
      R2=R**2
c      
      DO p=1,np
	DO k=1,ncb(p)
	ik1=numb(k,p)
        call  Coord (ik1,nx,ny,ix,iy,iz)
c
            do kr=1,nr2
             do jr=1,nr2
              do ir=1,nr2
c
            xc=x(ix(1))+(2*ir-1)*hxr
            yc=y(iy(1))+(2*jr-1)*hyr
            zc=z(iz(1))+(2*kr-1)*hzr
            vol=(hx*hy*hz)/nr23
          rho_raf=(xc-ax(p))**2+(yc-ay(p))**2+(zc-az(p))**2
c***********************************
!      print*, '--------------ATEENTION------------------'
!      print*, 'Note that SND_B_Neu uses the BAD option but it works...'
!      print*, '-----------------------------------------'
!            R2D =R2+0.5d0*(hx**2+hy**2+hz**2)
!            R2D1=R2-0.5d0*(hx**2+hy**2+hz**2)
!            if(rho_raf.le.R2D.and.rho_raf.gt.R2) then
!             vol=0.25d0*vol
!               elseif(rho_raf.le.R2.and.rho_raf.gt.R2D1)then
!             vol=0.5d0*vol  
!            endif
!          if (rho_raf.LE.R2D) then  
c************************************
c
          if (rho_raf.LE.R2) then  
c
         do j=1,8
          ixfd=ix(j)-bcV(1) !decalage du a l'elimination des sommets Dirichlet
          iyfd=iy(j)-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
          izfd=iz(j)-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
            F(izfd,iyfd,ixfd) = F(izfd,iyfd,ixfd)+vol*(
     > (interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,xc,yc,zc,x,y,z,U)*
     >  eval_GX_qi (j,ix,iy,iz,xc,yc,zc,x,y,z))     +
     > (interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,xc,yc,zc,x,y,z,U)*
     >  eval_GY_qi (j,ix,iy,iz,xc,yc,zc,x,y,z))     +
     > (interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,xc,yc,zc,x,y,z,U)*
     >  eval_GZ_qi (j,ix,iy,iz,xc,yc,zc,x,y,z)) 
     >                                                      )
         enddo
        endif
c
              enddo
             enddo
            enddo
         ENDDO
      ENDDO
C
         RETURN
         END
C=======================================================================C      
	SUBROUTINE SND_B_Neu_quad (n,nx,ny,nz,nc,nmx,nmy,nmz,np,ncb,R,
     >                ax,ay,az,x,y,z,numb,bcV,u,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER ncb(*),i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcV(6),
     >   ixfd,iyfd,izfd,np,p,numb(nc,np),nx,ny,nz,nmx,nmy,nmz
	DOUBLEPRECISION R,F(nmz,nmy,nmx),u(nmz,nmy,nmx),
     >           x(*),y(*),z(*),R2,ax(*),ay(*),az(*)
        doubleprecision xc,yc,zc,interpol_Q1_X,interpol_Q1_Y,
     >  interpol_Q1_Z,eval_GX_qi,eval_GY_qi,eval_GZ_qi,rho_raf
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
C
      R2=R**2
c
      DO p=1,np
	DO k=1,ncb(p)
c
        call  Coord (numb(k,p),nx,ny,ix,iy,iz)
c

      xmin=x(ix(1)) ; xmax=x(ix(7)) ; xc=0.5d0*(xmin+xmax)
      ymin=y(iy(1)) ; ymax=y(iy(7)) ; yc=0.5d0*(ymin+ymax)
      zmin=z(iz(1)) ; zmax=z(iz(7)) ; zc=0.5d0*(zmin+zmax)
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            rho_raf=(a(i)-ax(p))**2+(b(j)-ay(p))**2+(c(l)-az(p))**2
       if (rho_raf.LE.R2) then  
            ga(i,j,l)=(
     > (interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >             eval_GX_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z))     +
     > (interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >             eval_GY_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z))     +
     > (interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >  eval_GZ_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)) 
     >                                                      )
       else
            ga(i,j,l)=0.d0
       endif
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)-bcV(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+
     >                    integ_QN(n,dx,dy,dz,w0,w1,w2,ga)
c
      enddo
      enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE AssembleDXPv (n,nx,ny,nz,nc,nmx,nmy,nmz,coef,
     >                                       x,y,z,num,bcV,PR,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcV(6),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz
	DOUBLEPRECISION F(nmz,nmy,nmx),PR(nmz,nmy,nmx),x(*),y(*),z(*)
        doubleprecision interpol_Q1_X,evalqi,coef
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c

      xmin=x(ix(1)) ; xmax=x(ix(7)) 
      ymin=y(iy(1)) ; ymax=y(iy(7)) 
      zmin=z(iz(1)) ; zmax=z(iz(7)) 
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            ga(i,j,l)=evalqi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)*
     > interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,PR)
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcV(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+ (
     >                  coef*integ_QN(n,dx,dy,dz,w0,w1,w2,ga) )
c
       enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE AssembleDYPv (n,nx,ny,nz,nc,nmx,nmy,nmz,coef,
     >                                       x,y,z,num,bcV,PR,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcV(6),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz
	DOUBLEPRECISION F(nmz,nmy,nmx),PR(nmz,nmy,nmx),x(*),y(*),z(*)
        doubleprecision interpol_Q1_Y,evalqi,coef
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c

      xmin=x(ix(1)) ; xmax=x(ix(7)) 
      ymin=y(iy(1)) ; ymax=y(iy(7)) 
      zmin=z(iz(1)) ; zmax=z(iz(7)) 
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            ga(i,j,l)=evalqi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)*
     > interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,PR)
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcV(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+ (
     >                  coef*integ_QN(n,dx,dy,dz,w0,w1,w2,ga) )
c
       enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE AssembleDZPv (n,nx,ny,nz,nc,nmx,nmy,nmz,coef,
     >                                       x,y,z,num,bcV,PR,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcV(6),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz
	DOUBLEPRECISION F(nmz,nmy,nmx),PR(nmz,nmy,nmx),x(*),y(*),z(*)
        doubleprecision interpol_Q1_Z,evalqi,coef
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c

      xmin=x(ix(1)) ; xmax=x(ix(7)) 
      ymin=y(iy(1)) ; ymax=y(iy(7)) 
      zmin=z(iz(1)) ; zmax=z(iz(7)) 
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            ga(i,j,l)=evalqi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)*
     > interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,PR)
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcV(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+ (
     >                  coef*integ_QN(n,dx,dy,dz,w0,w1,w2,ga) )
c
       enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE AssembleDivUq (switch,n,nx,ny,nz,nc,nmx,nmy,nmz,coef,
     >                              xp,yp,zp,x,y,z,num,bcP,Ux,Uy,Uz,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcP(6),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz,
     >   ixf(8),iyf(8),izf(8)
        logical switch
	DOUBLEPRECISION F(nmz,nmy,nmx),Ux(nmz,nmy,nmx),Uy(nmz,nmy,nmx),
     >    Uz(nmz,nmy,nmx),x(*),y(*),z(*),xp(*),yp(*),zp(*)
        doubleprecision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z,
     >                  evalqi,coef
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c
      xmin=xp(ix(1)) ; xmax=xp(ix(7)) 
      ymin=yp(iy(1)) ; ymax=yp(iy(7)) 
      zmin=zp(iz(1)) ; zmax=zp(iz(7)) 
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
       if (switch) then
        do ks=1,8
          ixf(ks)=2*ix(ks)-1 
          iyf(ks)=2*iy(ks)-1  
          izf(ks)=2*iz(ks)-1
        enddo
       else
        do ks=1,8
          ixf(ks)=ix(ks)
          iyf(ks)=iy(ks)
          izf(ks)=iz(ks)
        enddo
       endif
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            ga(i,j,l)=evalqi (ks,ix,iy,iz,a(i),b(j),c(l),xp,yp,zp)*(
     > interpol_Q1_X (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Ux)+
     > interpol_Q1_Y (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Uy)+
     > interpol_Q1_Z (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Uz) ) 
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcP(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcP(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcP(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+ (
     >                  coef*integ_QN(n,dx,dy,dz,w0,w1,w2,ga) )
c
       enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE AssembleUGradq (switch,n,nx,ny,nz,nc,nmx,nmy,nmz,coef,
     >                              xp,yp,zp,x,y,z,num,bcP,Ux,Uy,Uz,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcP(6),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz,
     >   ixf(8),iyf(8),izf(8)
        logical switch
	DOUBLEPRECISION F(nmz,nmy,nmx),Ux(nmz,nmy,nmx),Uy(nmz,nmy,nmx),
     >    Uz(nmz,nmy,nmx),x(*),y(*),z(*),xp(*),yp(*),zp(*)
        doubleprecision interpol_Q1,coef
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn,GradQi(3)
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c
      xmin=xp(ix(1)) ; xmax=xp(ix(7)) 
      ymin=yp(iy(1)) ; ymax=yp(iy(7)) 
      zmin=zp(iz(1)) ; zmax=zp(iz(7)) 
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
       if (switch) then
        do ks=1,8
          ixf(ks)=2*ix(ks)-1 
          iyf(ks)=2*iy(ks)-1  
          izf(ks)=2*iz(ks)-1
        enddo
       else
        do ks=1,8
          ixf(ks)=ix(ks)
          iyf(ks)=iy(ks)
          izf(ks)=iz(ks)
        enddo
       endif
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
           call evalGradQi (ks,ix,iy,iz,a(i),b(j),c(l),xp,yp,zp,GradQi)
           ga(i,j,l)=
     >   interpol_Q1 (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Ux)*
     >      GradQi(1)   + 
     >   interpol_Q1 (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Uy)*
     >      GradQi(2)   + 
     >   interpol_Q1 (nmx,nmy,nmz,ixf,iyf,izf,a(i),b(j),c(l),x,y,z,Uz)*
     >      GradQi(3)

          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcP(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcP(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcP(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+ (
     >                  coef*integ_QN(n,dx,dy,dz,w0,w1,w2,ga) )
c
       enddo
      enddo
c
      RETURN
      END
C=======================================================================C
      subroutine smxv(n1,n2,n3,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,
     &                ch,f,ldf2,ldf3,u,ldu2,ldu3)
C=======================================================================C
      integer n1, n2, n3, ldf2, ldf3, ldu2, ldu3
      double precision a1(n1), b1(n1), c1(n1), d1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision f(ldf3,ldf2,n1), u(ldu3,ldu2,n1), ch
c
      integer i, j, k
c
      do i=1,n1
         do j=1,n2
            do k=1,n3
               f(k,j,i) = 0.d0
            end do
         end do
      end do
c
c i=1
c
      if (n1.gt.0) then
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(1)+ch*d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(1),f(1,1,1),ldf3,
     &        u(1,1,1),ldu3)
      end if
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(2)+ch*c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(2),f(1,1,1),ldf3,
     &        u(1,1,2),ldu3)
      end if
c
      do i=2,n1-1
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(i)+ch*c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(i),f(1,1,i),ldf3,
     &        u(1,1,i-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(i)+ch*d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(i),f(1,1,i),ldf3,
     &        u(1,1,i),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(i+1)+ch*c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(i+1),f(1,1,i),ldf3,
     &        u(1,1,i+1),ldu3)
      end do
c
c i=n1
c
      if (n1.gt.1) then
         call tmatxv(n2,n3,c2,d2,c3,d3,a1(n1)+ch*c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,c1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1-1),ldu3)
         call tmatxv(n2,n3,c2,d2,c3,d3,b1(n1)+ch*d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
         call tmatxv(n2,n3,a2,b2,c3,d3,d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
         call tmatxv(n2,n3,c2,d2,a3,b3,d1(n1),f(1,1,n1),ldf3,
     &        u(1,1,n1),ldu3)
      end if
c
      return
      end
C========================================================================
      FUNCTION eval_GX_qi (i,ix,iy,iz,x,y,z,a,b,c)
C
      IMPLICIT NONE
      INTEGER i,ix(8),iy(8),iz(8)
      DOUBLE PRECISION H1,ax,ay,az,bx,by,bz,x,y,z,
     >                 a(*),b(*),c(*),eval_GX_qi
C
      eval_GX_qi=0.D0
c
      ax=a(ix(1))
      ay=b(iy(1))     
      az=c(iz(1))
      bx=a(ix(7))
      by=b(iy(7))     
      bz=c(iz(7))
c
      if ((bx.eq.ax).or.(by.eq.ay).or.(bz.eq.az)) then
         print*, 'error in eval_GX_qi... nodes a and b coincide'
         stop
      endif
c
      H1=1.D0/((bx-ax)*(by-ay)*(bz-az))
      IF(i.EQ.1) eval_GX_qi=-H1*(by-y)*(bz-z)
      IF(i.EQ.2) eval_GX_qi= H1*(by-y)*(bz-z)
      IF(i.EQ.3) eval_GX_qi= H1*(y-ay)*(bz-z)
      IF(i.EQ.4) eval_GX_qi=-H1*(y-ay)*(bz-z)
      IF(i.EQ.5) eval_GX_qi=-H1*(by-y)*(z-az)
      IF(i.EQ.6) eval_GX_qi= H1*(by-y)*(z-az)
      IF(i.EQ.7) eval_GX_qi= H1*(y-ay)*(z-az)
      IF(i.EQ.8) eval_GX_qi=-H1*(y-ay)*(z-az)
C
      RETURN
      END
C========================================================================
      FUNCTION eval_GY_qi (i,ix,iy,iz,x,y,z,a,b,c)
C
      IMPLICIT NONE
      INTEGER i,ix(8),iy(8),iz(8)
      DOUBLE PRECISION H1,ax,ay,az,bx,by,bz,x,y,z,
     >                 a(*),b(*),c(*),eval_GY_qi
C
      eval_GY_qi=0.D0
c
      ax=a(ix(1))
      ay=b(iy(1))     
      az=c(iz(1))
      bx=a(ix(7))
      by=b(iy(7))     
      bz=c(iz(7))
c
      if ((bx.eq.ax).or.(by.eq.ay).or.(bz.eq.az)) then
         print*, 'error in eval_GY_qi... nodes a and b coincide'
         stop
      endif
c
      H1=1.D0/((bx-ax)*(by-ay)*(bz-az))
      IF(i.EQ.1) eval_GY_qi=-H1*(bx-x)*(bz-z)
      IF(i.EQ.2) eval_GY_qi=-H1*(x-ax)*(bz-z)
      IF(i.EQ.3) eval_GY_qi= H1*(x-ax)*(bz-z)
      IF(i.EQ.4) eval_GY_qi= H1*(bx-x)*(bz-z)
      IF(i.EQ.5) eval_GY_qi=-H1*(bx-x)*(z-az)
      IF(i.EQ.6) eval_GY_qi=-H1*(x-ax)*(z-az)
      IF(i.EQ.7) eval_GY_qi= H1*(x-ax)*(z-az)
      IF(i.EQ.8) eval_GY_qi= H1*(bx-x)*(z-az)
C
      RETURN
      END
C========================================================================
      FUNCTION eval_GZ_qi (i,ix,iy,iz,x,y,z,a,b,c)
C
      IMPLICIT NONE
      INTEGER i,ix(8),iy(8),iz(8)
      DOUBLE PRECISION H1,ax,ay,az,bx,by,bz,x,y,z,
     >                 a(*),b(*),c(*),eval_GZ_qi
C
      eval_GZ_qi=0.D0
c
      ax=a(ix(1))
      ay=b(iy(1))     
      az=c(iz(1))
      bx=a(ix(7))
      by=b(iy(7))     
      bz=c(iz(7))
c
      if ((bx.eq.ax).or.(by.eq.ay).or.(bz.eq.az)) then
         print*, 'error in eval_GZ_qi... nodes a and b coincide'
         stop
      endif
c
      H1=1.D0/((bx-ax)*(by-ay)*(bz-az))
c
      IF(i.EQ.1) eval_GZ_qi=-H1*(bx-x)*(by-y)
      IF(i.EQ.2) eval_GZ_qi=-H1*(x-ax)*(by-y)
      IF(i.EQ.3) eval_GZ_qi=-H1*(x-ax)*(y-ay)
      IF(i.EQ.4) eval_GZ_qi=-H1*(bx-x)*(y-ay)
      IF(i.EQ.5) eval_GZ_qi= H1*(bx-x)*(by-y)
      IF(i.EQ.6) eval_GZ_qi= H1*(x-ax)*(by-y)
      IF(i.EQ.7) eval_GZ_qi= H1*(x-ax)*(y-ay)
      IF(i.EQ.8) eval_GZ_qi= H1*(bx-x)*(y-ay)
C
      RETURN
      END
*=======================================================



C=======================================================================C      
	SUBROUTINE test_quad (n,nx,ny,nz,nc,nmx,nmy,nmz,
     >                ax,ay,az,x,y,z,R,num,bcV,u,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,n,nc,ix(8),iy(8),iz(8),bcV(6),
     >   ixfd,iyfd,izfd,num(nc),nx,ny,nz,nmx,nmy,nmz
	DOUBLEPRECISION F(nmz,nmy,nmx),u(nmz,nmy,nmx),R,
     >           x(*),y(*),z(*),R2,ax(*),ay(*),az(*)
        doubleprecision xc,yc,zc,interpol_Q1_X,interpol_Q1_Y,
     >  interpol_Q1_Z,eval_GX_qi,eval_GY_qi,eval_GZ_qi,rho_raf
        double precision xmin,xmax,ymin,ymax,zmin,zmax,rt,rt2,w0,w1,w2,
     >  dx,dy,dz,a(5),b(5),c(5),ga(5,5,5),integ_qn
C
      R2=R**2
c
      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c

      xmin=x(ix(1)) ; xmax=x(ix(7)) ; xc=0.5d0*(xmin+xmax)
      ymin=y(iy(1)) ; ymax=y(iy(7)) ; yc=0.5d0*(ymin+ymax)
      zmin=z(iz(1)) ; zmax=z(iz(7)) ; zc=0.5d0*(zmin+zmax)
c
       call Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c
      do ks=1,8
c
       do l=1,n
        do j=1,n
          do i=1,n
            rho_raf=(a(i)-ax(1))**2+(b(j)-ay(1))**2+(c(l)-az(1))**2
!       if (rho_raf.LE.R2) then  
            ga(i,j,l)=(
     > (interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >             eval_GX_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z))     +
     > (interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >             eval_GY_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z))     +
     > (interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,a(i),b(j),c(l),x,y,z,U)*
     >  eval_GZ_qi (ks,ix,iy,iz,a(i),b(j),c(l),x,y,z)) 
     >                                                      )
!       else
!            ga(i,j,l)=0.d0
!       endif
          enddo
        enddo
       enddo
c
         ixfd=ix(ks)!-bcV(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcV(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcV(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >                F(izfd,iyfd,ixfd)+
     >                    integ_QN(n,dx,dy,dz,w0,w1,w2,ga)
c
      enddo
      enddo
c
      RETURN
      END
