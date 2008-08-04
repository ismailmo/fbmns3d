
c Author : Mourad Ismail   (ismail@ujf-grenoble.fr)

C--------------------------------------------------------------------C
C--------------------------FUNCTIONS---------------------------------C
C--------------------------------------------------------------------C
C
C======================================================================C
      function interpol_Q1 (nmx,nmy,nmz,ix,iy,iz,x,y,z,xx,yy,zz,U)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j
      integer ix(8),iy(8),iz(8)
      double precision x,y,z,xx(*),yy(*),zz(*),interpol_Q1,evalqi
      double precision U(nmz,nmy,nmx)
c
      interpol_Q1=0.D0
c
      do j=1,8
      interpol_Q1=interpol_Q1+u(iz(j),iy(j),ix(j))*
     > evalqi (j,ix,iy,iz,x,y,z,xx,yy,zz)
      enddo
c
      return
      end
C
C======================================================================C
      function interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,x,y,z,xx,yy,zz,U)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j
      integer ix(8),iy(8),iz(8)
      double precision x,y,z,xx(*),yy(*),zz(*),interpol_Q1_X,eval_GX_qi
      double precision U(nmz,nmy,nmx)
c
      interpol_Q1_X=0.D0
c
      do j=1,8
      interpol_Q1_X=interpol_Q1_X+u(iz(j),iy(j),ix(j))*
     > eval_GX_qi (j,ix,iy,iz,x,y,z,xx,yy,zz)
      enddo
c
      return
      end
C======================================================================C
      function interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,x,y,z,xx,yy,zz,U)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j
      integer ix(8),iy(8),iz(8)
      double precision x,y,z,xx(*),yy(*),zz(*),interpol_Q1_Y,eval_GY_qi
      double precision U(nmz,nmy,nmx)
c
      interpol_Q1_Y=0.D0
c
      do j=1,8
      interpol_Q1_Y=interpol_Q1_Y+u(iz(j),iy(j),ix(j))*
     > eval_GY_qi (j,ix,iy,iz,x,y,z,xx,yy,zz)
      enddo
c
      return
      end
C======================================================================C
      function interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,x,y,z,xx,yy,zz,U)
C======================================================================C
      implicit none
      integer nmx,nmy,nmz,j
      integer ix(8),iy(8),iz(8)
      double precision x,y,z,xx(*),yy(*),zz(*),interpol_Q1_Z,eval_GZ_qi
      double precision U(nmz,nmy,nmx)
c
      interpol_Q1_Z=0.D0
c
      do j=1,8
      interpol_Q1_Z=interpol_Q1_Z+u(iz(j),iy(j),ix(j))*
     > eval_GZ_qi (j,ix,iy,iz,x,y,z,xx,yy,zz)
      enddo
c
      return
      end
C========================================================================
      FUNCTION evalqi (i,ix,iy,iz,x,y,z,a,b,c)
C========================================================================
C
      IMPLICIT NONE
      INTEGER i,ix(8),iy(8),iz(8)
      DOUBLE PRECISION H1,ax,ay,az,bx,by,bz,x,y,z,a(*),b(*),c(*),evalqi
C
      evalqi=0.D0
c
      ax=a(ix(1))
      ay=b(iy(1))     
      az=c(iz(1))
      bx=a(ix(7))
      by=b(iy(7))     
      bz=c(iz(7))
c
      H1=1.D0/((bx-ax)*(by-ay)*(bz-az))
c
      IF(i.EQ.1) evalqi=H1*(bx-x)*(by-y)*(bz-z)
      IF(i.EQ.2) evalqi=H1*(x-ax)*(by-y)*(bz-z)
      IF(i.EQ.3) evalqi=H1*(x-ax)*(y-ay)*(bz-z)
      IF(i.EQ.4) evalqi=H1*(bx-x)*(y-ay)*(bz-z)
      IF(i.EQ.5) evalqi=H1*(bx-x)*(by-y)*(z-az)
      IF(i.EQ.6) evalqi=H1*(x-ax)*(by-y)*(z-az)
      IF(i.EQ.7) evalqi=H1*(x-ax)*(y-ay)*(z-az)
      IF(i.EQ.8) evalqi=H1*(bx-x)*(y-ay)*(z-az)
C
      RETURN
      END
C========================================================================
      SUBROUTINE evalGradQi (i,ix,iy,iz,x,y,z,a,b,c,GradQi)
C========================================================================
C
      IMPLICIT NONE
      INTEGER i,ix(8),iy(8),iz(8)
      DOUBLE PRECISION H1,ax,ay,az,bx,by,bz,x,y,z,a(*),b(*),c(*),
     >     GradQi(3)
C
      GradQi(1)=0.D0 ; GradQi(2)=0.D0 ; GradQi(3)=0.D0
c
      ax=a(ix(1))
      ay=b(iy(1))     
      az=c(iz(1))
      bx=a(ix(7))
      by=b(iy(7))     
      bz=c(iz(7))
c
      H1=1.D0/((bx-ax)*(by-ay)*(bz-az))
c
      IF(i.EQ.1) then
        GradQi(1)=-H1*(by-y)*(bz-z)
        GradQi(2)=-H1*(bx-x)*(bz-z)
        GradQi(3)=-H1*(bx-x)*(by-y)
        return
      endif
      IF(i.EQ.2) then 
        GradQi(1)= H1*(by-y)*(bz-z)
        GradQi(2)=-H1*(x-ax)*(bz-z)
        GradQi(3)=-H1*(x-ax)*(by-y)
        return
      endif
      IF(i.EQ.3) then 
        GradQi(1)= H1*(y-ay)*(bz-z)
        GradQi(2)= H1*(x-ax)*(bz-z)
        GradQi(3)=-H1*(x-ax)*(y-ay)
        return
      endif
      IF(i.EQ.4) then 
        GradQi(1)=-H1*(y-ay)*(bz-z)
        GradQi(2)= H1*(bx-x)*(bz-z)
        GradQi(3)=-H1*(bx-x)*(y-ay)
        return
      endif
      IF(i.EQ.5) then
        GradQi(1)=-H1*(by-y)*(z-az)
        GradQi(2)=-H1*(bx-x)*(z-az)
        GradQi(3)= H1*(bx-x)*(by-y)
        return
      endif
      IF(i.EQ.6) then 
        GradQi(1)= H1*(by-y)*(z-az)
        GradQi(2)=-H1*(x-ax)*(z-az)
        GradQi(3)= H1*(x-ax)*(by-y)
        return
      endif
      IF(i.EQ.7) then 
        GradQi(1)= H1*(y-ay)*(z-az)
        GradQi(2)= H1*(x-ax)*(z-az)
        GradQi(3)= H1*(x-ax)*(y-ay)
        return
      endif
      IF(i.EQ.8) then 
        GradQi(1)=-H1*(y-ay)*(z-az)
        GradQi(2)= H1*(bx-x)*(z-az)
        GradQi(3)= H1*(bx-x)*(y-ay)
        return
      endif
C
      RETURN
      END
c=================================================================C
C
C-----------------------------------------------------------------C
C--------------------------ROUTINES-------------------------------C
C-----------------------------------------------------------------C
C
c=================================================================C
      SUBROUTINE INTERPOL_C(np,ncs,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >             nrep,num,nums,ax,ay,az,xx,yy,zz,SS,U,uinter)
      IMPLICIT NONE
      integer nmx,nmy,nmz,nx,ny,nz
      INTEGER ncs,i,ii,j,p,np,trans1,i1,i2,ip
      INTEGER ix(8),iy(8),iz(8),nrep(*),num(*),nums(8,ncs),knrep
      DOUBLE PRECISION evalqi,x,y,z,hx,hy,hz,x1,y1,z1,x3,y3,z3,R,rho,
     > xx(*),yy(*),zz(*),SS(3,*),U(nmz,nmy,nmx),uinter(*),
     > ax(*),ay(*),az(*)
C
!Le tableau Uinter contient l'interpole de la solutions aux centres
!des elements de gamma et de gamma prime de chaque particule.
!pour une particule p, ces valeurs sont stockees de ceete maniere
!interpole sur gamma      <--- Uinter(i), i=2(p-1)ncs+1,.....,(2p-1)ncs
!interpole sur gammaPrime <--- Uinter(i), i=(2p-1)ncs+1,.....,2p*ncs
C
C!interpolation de U sur Gamma!
      DO p=1,np
         ip=2*(p-1)*ncs
         i1=ip+1 ; i2=ip+ncs
      DO i=i1,i2
       ii=i-ip
       x1=SS(1,nums(1,ii))
       y1=SS(2,nums(1,ii))
       z1=SS(3,nums(1,ii))
       x3=SS(1,nums(3,ii))
       y3=SS(2,nums(3,ii))
       z3=SS(3,nums(3,ii))
c
      x=0.5D0*(x1+x3) ; y=0.5D0*(y1+y3) ; z=0.5D0*(z1+z3)
      rho=dsqrt(x**2+y**2+z**2)
      x=ax(p)+x*R/rho ; y=ay(p)+y*R/rho ; z=az(p)+z*R/rho
c
       uinter(i)=0.
       trans1=int(ax(p)/hx)+int(ay(p)/hy)*(nx-1)+
     >        int(az(p)/hz)*(nx-1)*(ny-1)
c
	 knrep=num(nrep(i-ip)+trans1)
         call Coord (knrep,nx,ny,ix,iy,iz)
c        
!       if((xx(ix(1)).le.x.and.x.le.xx(ix(7))).and. !??????????????
!     >       (yy(iy(1)).le.y.and.y.le.yy(iy(7))).and.
!     >       (zz(iz(1)).le.z.and.z.le.zz(iz(7))) ) then
       DO j=1,8
        Uinter(i)=Uinter(i)+u(iz(j),iy(j),ix(j))*
     >                     evalqi(j,ix,iy,iz,x,y,z,xx,yy,zz)
       ENDDO
!      else !??????????????
!         print*,'Possible error in interpolate_c routine.'
!         print*,'Bad localisation in element number ', i, '?'
!         print*,'X--->',xx(ix(1)),x,xx(ix(7))
!         print*,'Y--->',yy(iy(1)),y,yy(iy(7))
!         print*,'Z--->',zz(iz(1)),z,zz(iz(7))
!         stop
!      endif !??????????????
      ENDDO
      ENDDO
C!interpolation de U sur Gamma prime!
      DO p=1,np
         ip=(2*p-1)*ncs
         i1=ip+1 ; i2=ip+ncs
      DO i=i1,i2
       Uinter(i)=0.
       ii=i-ip
c
       x1=SS(1,nums(1,ii))
       y1=SS(2,nums(1,ii))
       z1=SS(3,nums(1,ii))
       x3=SS(1,nums(3,ii))
       y3=SS(2,nums(3,ii))
       z3=SS(3,nums(3,ii))
c
      x=0.5D0*(x1+x3) ; y=0.5D0*(y1+y3) ; z=0.5D0*(z1+z3)
      rho=dsqrt(x**2+y**2+z**2)
      x=ax(p)+x*R/rho ; y=ay(p)+y*R/rho ; z=az(p)+z*R/rho
c
       trans1=int(ax(p)/hx)+int(ay(p)/hy)*(nx-1)+
     >        int(az(p)/hz)*(nx-1)*(ny-1)       
c
       knrep=num(nrep(i-ip)+trans1)
        call Coord (knrep,nx,ny,ix,iy,iz)
       DO j=1,8
        Uinter(i)=Uinter(i)+u(iz(j),iy(j),ix(j))*
     > evalqi(j,ix,iy,iz,x,y,z,xx,yy,zz)
       ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C====================================================================C
        SUBROUTINE U_Grad (nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,U,
     >                      D_U_X,D_U_Y,D_U_Z )
C====================================================================C
        IMPLICIT NONE
        INTEGER nx,ny,nz,nmx,nmy,nmz,i,j,k
        DOUBLEPRECISION U(nmz,nmy,nmz),D_U_X(nmz,nmy,nmz),
     >   D_U_Y(nmz,nmy,nmz),D_U_Z(nmz,nmy,nmz),hx,hy,hz
C
!      do k=1,nz
!         do j=1,ny
!            do i=1,nx
!             if (i.eq.1) then
!               D_U_X(k,j,i)=(U( k , j ,i+1)-U( k , j , i ))/hx
!              elseif (i.eq.nx) then
!               D_U_X(k,j,i)=(U( k , j ,i  )-U( k , j ,i-1))/hx
!               else
!              D_U_X(k,j,i)=(U( k , j ,i+1)-U( k , j ,i-1))/(2*hx)
!            endif
!             if (j.eq.1) then
!               D_U_Y(k,j,i)=(U( k , j+1 ,i)-U( k ,  j , i ))/hy
!              elseif (j.eq.ny) then
!               D_U_Y(k,j,i)=(U( k , j ,i  )-U( k , j-1 ,i ))/hy
!               else
!               D_U_Y(k,j,i)=(U( k ,j+1, i )-U( k ,j-1, i ))/(2*hy)
!             endif
!             if (k.eq.1) then
!               D_U_Z(k,j,i)=(U( k+1 , j ,i)-U( k , j , i ))/hz
!              elseif (k.eq.nz) then
!               D_U_Z(k,j,i)=(U( k , j ,i  )-U( k-1 , j ,i))/hz
!               else
!               D_U_Z(k,j,i)=(U(k+1, j , i )-U(k-1, j , i ))/(2*hz)
!             endif
!            end do
!         end do
!      end do
C
      do k=1,nz-1
         do j=2,ny-1
            do i=2,nx-1
               D_U_X(k,j,i)=(U( k , j ,i+1)-U( k , j ,i-1))/(2*hx)
               D_U_Y(k,j,i)=(U( k ,j+1, i )-U( k ,j-1, i ))/(2*hy)
               if (k.eq.1) then
               D_U_Z(k,j,i)=(U(2, j , i )-U(1, j , i ))/(hz)
               else
               D_U_Z(k,j,i)=(U(k+1, j , i )-U(k-1, j , i ))/(2*hz)
               endif
            end do
         end do
      end do

        RETURN
        END
C===================================================================C
      SUBROUTINE PROJECT_GF (nx,ny,nz,px,py,pz,nmx,nmy,nmz,numP,
     >                           x,y,z,xp,yp,zp,Ug,Uf)
C========================================================================
C
      IMPLICIT NONE
      integer nx,ny,nz,px,py,pz,nmx,nmy,nmz
      integer i1,i2,j1,j2,k1,k2,NcP
      INTEGER i,j,k,knrep,ix(8),iy(8),iz(8),numP(*)
      DOUBLE PRECISION x(*),y(*),z(*),interpol_Q1,a,b,c
      double precision xp(*),yp(*),zp(*)
      double precision Uf(nmz,nmy,nmx),Ug(nmz,nmy,nmx)
C
      NcP=(px-1)*(py-1)*(pz-1)
c
      do knrep=1,NcP
         call Coord (numP(knrep),px,py,ix,iy,iz)
         k1=2*iz(1)-1 ; k2=2*iz(1)+1
         j1=2*iy(1)-1 ; j2=2*iy(1)+1
         i1=2*ix(1)-1 ; i2=2*ix(1)+1
c
         do k=k1,k2
            c=z(k)
            do j=j1,j2
               b=y(j)
               do i=i1,i2
                  a=x(i)
          Uf(k,j,i)=interpol_Q1(nmx,nmy,nmz,ix,iy,iz,
     >           a,b,c,xp,yp,zp,Ug)
               enddo
            enddo
         enddo
       enddo
c
      RETURN
      END
C===========================================================================C
      Function Integ_G(ncs, area, DnV)
C===========================================================================C
      IMPLICIT NONE
      INTEGER ncs,k
      DOUBLE PRECISION Integ_G, area(ncs), DnV(*)
C
         Integ_G = 0.
      do k = 1 , ncs
         Integ_G = Integ_G + area(k) * DnV(k)
      enddo
C
      RETURN
      END
C===========================================================================C
      subroutine Integ_GP(ncs, area, nums, SS, R, Pinter,
     >  Integ_GPX,Integ_GPY,Integ_GPZ)
C===========================================================================C
      IMPLICIT NONE
      INTEGER ncs,k,nums(8,ncs)
      doubleprecision SS(3,*)
      doubleprecision x,y,z,rho,x1,y1,z1,x3,y3,z3,R
      DOUBLE PRECISION Integ_GPX,Integ_GPY,Integ_GPZ
      doubleprecision  area(ncs), Pinter(*)
C
         Integ_GPX = 0. ; Integ_GPY = 0. ; Integ_GPZ = 0.
      do k = 1 , ncs
       x1=SS(1,nums(1,k))
       y1=SS(2,nums(1,k))
       z1=SS(3,nums(1,k))
       x3=SS(1,nums(3,k))
       y3=SS(2,nums(3,k))
       z3=SS(3,nums(3,k))
      x=0.5D0*(x1+x3) ; y=0.5D0*(y1+y3) ; z=0.5D0*(z1+z3)
      rho=dsqrt(x**2+y**2+z**2)
      
!      x=ax(1)+x*R/rho ; y=ay(1)+y*R/rho ; z=az(1)+z*R/rho
      x=x*R/rho ; y=y*R/rho ; z=z*R/rho
c
         Integ_GPX = Integ_GPX + area(k) * Pinter(k) * x
         Integ_GPY = Integ_GPY + area(k) * Pinter(k) * y
         Integ_GPZ = Integ_GPZ + area(k) * Pinter(k) * z
c
      enddo
C
      RETURN
      END
C===========================================================================C
      SUBROUTINE INTERPOL_S(np,ncs,nss,nx,ny,nz,nmx,nmy,nmz,hx,hy,hz,R,
     >             nrep,num,nums,ax,ay,az,xx,yy,zz,SS,U,uinter)
C===========================================================================C
      IMPLICIT NONE
      integer nmx,nmy,nmz,nx,ny,nz
      INTEGER ncs,nss,i,ii,j,p,np,trans1,i1,i2,ip
      INTEGER ix(8),iy(8),iz(8),nrep(*),num(*),nums(8,ncs),knrep
      DOUBLE PRECISION evalqi,x,y,z,hx,hy,hz,R,
     > xx(*),yy(*),zz(*),SS(3,*),U(nmz,nmy,nmx),uinter(*),
     > ax(*),ay(*),az(*)
C
      DO p=1,np
         ip=2*(p-1)*nss
         i1=ip+1 ; i2=ip+nss
         DO i=i1,i2
            ii=i-ip
            x=ax(p)+SS(1,ii)
            y=ay(p)+SS(2,ii)
            z=az(p)+SS(3,ii)
c     
            uinter(i)=0.
            trans1=int(ax(p)/hx)+int(ay(p)/hy)*(nx-1)+
     >           int(az(p)/hz)*(nx-1)*(ny-1)
c     
            knrep=num(nrep(i-ip)+trans1)
            call Coord (knrep,nx,ny,ix,iy,iz)
c     
            DO j=1,8
               Uinter(i)=Uinter(i)+u(iz(j),iy(j),ix(j))*
     >              evalqi(j,ix,iy,iz,x,y,z,xx,yy,zz)
            ENDDO
         ENDDO
      ENDDO
c
      return
      end
C=======================================================================C      
	SUBROUTINE ComputeDivU (nx,ny,nz,nc,nmx,nmy,nmz,
     >                              xp,yp,zp,num,Ux,Uy,Uz,f)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,nc,ix(8),iy(8),iz(8),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz,
     >   ixf(8),iyf(8),izf(8)
        logical switch
	DOUBLEPRECISION F(nmz,nmy,nmx),Ux(nmz,nmy,nmx),Uy(nmz,nmy,nmx),
     >    Uz(nmz,nmy,nmx),xp(*),yp(*),zp(*)
        doubleprecision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z
c

      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c
      do ks=1,8
c
         ixfd=ix(ks)!-bcP(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcP(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcP(5) !decalage du a l'elimination des sommets Dirichlet 
c
         F(izfd,iyfd,ixfd)=
     >        interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Ux) +
     >        interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uy) +
     >        interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uz)
c
!eval_GX_qi (ks,ix,iy,iz,x,y,z,xx,yy,zz)
      enddo
      enddo
c
      RETURN
      END

C=======================================================================C      
	SUBROUTINE ComputeStress (nc,nu,
     >     dxUx,dyUx,dzUx,dxUy,dyUy,dzUy,dxUz,dyUz,dzUz,
     >     norx,nory,norz,SigmaX,SigmaY,SigmaZ,diver)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER k,nc
	DOUBLEPRECISION nu,dxUx(*),dyUx(*),
     >       dzUx(*),dxUy(*),dyUy(*),
     >       dzUy(*),dxUz(*),dyUz(*),
     >       dzUz(*),diver(*),
     >       SigmaX(*),SigmaY(*),SigmaZ(*),norx(*),nory(*),norz(*)
c
        do k = 1 , nc
           SigmaX(k) = nu*( 2.d0*dxUx(k)*norx(k) +
     >          (dyUx(k) + dxUy(k))*nory(k) +
     >          (dzUx(k) + dxUz(k))*norz(k) )
c
           SigmaY(k) = nu*( 2.d0*dyUy(k)*nory(k) +
     >          (dxUy(k) + dyUx(k))*norx(k) +
     >          (dzUy(k) + dyUz(k))*norz(k) )
c
           SigmaZ(k) = nu*( 2.d0*dzUz(k)*norz(k) +
     >          (dxUz(k) + dzUx(k))*norx(k) +
     >          (dyUz(k) + dzUy(k))*nory(k) )
!!!!!!!!!
         diver(k)=dxUx(k) + dyUy(k) + dzUz(k)
        enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE xyzDerivatives (nx,ny,nz,nc,nmx,nmy,nmz,
     >     xp,yp,zp,num,Ux,Uy,Uz,
     >     dxUx,dyUx,dzUx,dxUy,dyUy,dzUy,dxUz,dyUz,dzUz)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,nc,ix(8),iy(8),iz(8),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz,
     >   ixf(8),iyf(8),izf(8)
        logical switch
	DOUBLEPRECISION Ux(nmz,nmy,nmx),Uy(nmz,nmy,nmx),
     >       Uz(nmz,nmy,nmx),xp(*),yp(*),zp(*)
        doubleprecision dxUx(nmz,nmy,nmx),dyUx(nmz,nmy,nmx),
     >       dzUx(nmz,nmy,nmx)
        doubleprecision dxUy(nmz,nmy,nmx),dyUy(nmz,nmy,nmx),
     >       dzUy(nmz,nmy,nmx)
        doubleprecision dxUz(nmz,nmy,nmx),dyUz(nmz,nmy,nmx),
     >       dzUz(nmz,nmy,nmx)
        doubleprecision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z
c

      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c
      do ks=1,8
c
         ixfd=ix(ks)!-bcP(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcP(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcP(5) !decalage du a l'elimination des sommets Dirichlet 
ccc--------Ux------------
         dxUx(izfd,iyfd,ixfd)=
     >        interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Ux)
         dyUx(izfd,iyfd,ixfd)=
     >        interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Ux)
         dzUx(izfd,iyfd,ixfd)=
     >        interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Ux)
ccc--------Uy------------
         dxUy(izfd,iyfd,ixfd)=
     >        interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uy)
         dyUy(izfd,iyfd,ixfd)=
     >        interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uy)
         dzUy(izfd,iyfd,ixfd)=
     >        interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uy)
ccc--------Uz------------
         dxUz(izfd,iyfd,ixfd)=
     >        interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uz)
         dyUz(izfd,iyfd,ixfd)=
     >        interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uz)
         dzUz(izfd,iyfd,ixfd)=
     >        interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,Uz)
      enddo
      enddo
c
      RETURN
      END
C=======================================================================C      
	SUBROUTINE gradientPressure (nx,ny,nz,nc,nmx,nmy,nmz,
     >     xp,yp,zp,num,P,
     >     dxP,dyP,dzP)
C=======================================================================C	
	IMPLICIT NONE
	INTEGER i,j,k,ks,l,nc,ix(8),iy(8),iz(8),
     >   ixfd,iyfd,izfd,num(*),nx,ny,nz,nmx,nmy,nmz,
     >   ixf(8),iyf(8),izf(8)
        logical switch
	DOUBLEPRECISION xp(*),yp(*),zp(*)
        doubleprecision dxP(nmz,nmy,nmx),dyP(nmz,nmy,nmx),
     >       dzP(nmz,nmy,nmx),P(nmz,nmy,nmx)
        doubleprecision interpol_Q1_X,interpol_Q1_Y,interpol_Q1_Z
c

      DO k=1,nc
c
        call  Coord (num(k),nx,ny,ix,iy,iz)
c
      do ks=1,8
c
         ixfd=ix(ks)!-bcP(1) !decalage du a l'elimination des sommets Dirichlet
         iyfd=iy(ks)!-bcP(3) !decalage du a l'elimination des sommets Dirichlet 
         izfd=iz(ks)!-bcP(5) !decalage du a l'elimination des sommets Dirichlet 
ccc--------Ux------------
         dxP(izfd,iyfd,ixfd)=
     >        interpol_Q1_X (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,P)
         dyP(izfd,iyfd,ixfd)=
     >        interpol_Q1_Y (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,P)
         dzP(izfd,iyfd,ixfd)=
     >        interpol_Q1_Z (nmx,nmy,nmz,ix,iy,iz,
     >        xp(ixfd),yp(iyfd),zp(izfd),xp,yp,zp,P)
      enddo
      enddo
c
      RETURN
      END
