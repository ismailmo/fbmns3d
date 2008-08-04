
c Author : Mourad Ismail   (ismail@ujf-grenoble.fr)

*=======================================================================*
	SUBROUTINE MSH_OMG (nx,ny,nz,hx,hy,hz,hv,
     >                 xi,yi,zi,xl,yl,zl,x,y,z,NUM)      
*=======================================================================*
	IMPLICIT NONE
	INTEGER nx,ny,nz,i,j,k,l
	INTEGER num(*)
	DOUBLE PRECISION xl,yl,zl,xi,yi,zi
	DOUBLE PRECISION x(*),y(*),z(*), hx,hy,hz,hv
C
      x(1)=xi ;  y(1)=yi ;  z(1)=zi
c
      hx = xl/(nx-1)
      do i=1,nx-1
       x(i+1)=x(i)+hx
      end do
c
      hy = yl/(ny-1)
      do j=1,ny-1
       y(j+1)=y(j)+hy
      end do
c
      hz = zl/(nz-1)
      do j=1,nz-1
       z(j+1)=z(j)+hz
      end do
C
      hv = hx*hy*hz
      k=1
      DO l=1,nz-1
       DO j=1,ny-1
        DO i=1,nx-1
         num(k)=i+(j-1)*nx+(l-1)*nx*ny
	 k=k+1
        ENDDO
       ENDDO
      ENDDO
C
      return
      end
*=======================================================================*
	SUBROUTINE MSH_sphere (nn,nss,ncs,r,SS,NUMS)      
*=======================================================================*
	IMPLICIT NONE
	INTEGER nn,nss,ncs,i,j,k,kk,ndep,nmoins,nfin
	INTEGER nums(8,ncs)
	DOUBLE PRECISION r,pi
	DOUBLE PRECISION SS(3,2*nss)
C
	pi=dacos(-1.d0)
C-----------------------------------
	!maillage de la sphere
C-----------------------------------
C
        CALL MAILLREF (nn,nss,ncs,r,SS,nums)
        CALL ROTATION (1,pi/2,nss,nn*nn+1,nn*(2*nn-1),nn*(nn-1),SS)
        CALL RENUM    (1,nn,ncs,(nn-1)*(nn-1)+1,0,
     >                                             2*(nn-1)*(nn-1),nums)
        CALL ROTATION (1,pi,nss,nn*(2*nn-1)+1,nn*(3*nn-2),
     >                                                   2*nn*(nn-1),SS)
        CALL RENUM    (1,nn,ncs,2*(nn-1)*(nn-1)+1,0,
     >                                             3*(nn-1)*(nn-1),nums)	
        CALL ROTATION (1,3*pi/2,nss,nn*(3*nn-2)+1,4*nn*(nn-1),
     >                                                   3*nn*(nn-1),SS)
        CALL RENUM    (2,nn,ncs,3*(nn-1)*(nn-1)+1,
     >                4*(nn-1)*(nn-1)-(nn-1)+1,(nn-1)*(4*(nn-1)-1),nums)  
	DO i=1,nn-2
	 ndep  = 4*nn*(nn-1)+1+(i-1)*(nn-2)
	 nfin  = 4*nn*(nn-1)+i*(nn-2)
	 nmoins= 4*nn*(nn-1)-(nn+1)-2*(i-1)
         CALL ROTATION (2,pi/2,nss,ndep,nfin,nmoins,SS)        
     	ENDDO
	!couche1
	nums(1,4*(nn-1)*(nn-1)+1)=nn
	nums(2,4*(nn-1)*(nn-1)+1)=nn*(nn+2*(nn-1)+nn-2)	
	nums(3,4*(nn-1)*(nn-1)+1)=4*nn*(nn-1)+1
	nums(4,4*(nn-1)*(nn-1)+1)=nums(1,4*(nn-1)*(nn-1)+1)+nn
	DO k=4*(nn-1)*(nn-1)+2,4*(nn-1)*(nn-1)+nn-1
	 nums(1,k)=nums(2,k-1)
	 nums(2,k)=nums(1,k)-nn
	IF (k.EQ.4*(nn-1)*(nn-1)+nn-1) THEN
	 nums(3,k)=nums(2,k)-nn
	ELSE
	 nums(3,k)=nums(3,k-1)+1
	ENDIF
	 nums(4,k)=nums(3,k-1)
	ENDDO
	!couches internes
	kk=4*(nn-1)*(nn-1)+nn
	DO i=1,nn-3
	nums(1,kk)=nums(4,kk-(nn-1))
	nums(2,kk)=nums(3,kk-(nn-1))
	nums(3,kk)=nums(2,kk)+nn-2
	nums(4,kk)=nums(1,kk)+nn
	DO k=kk+1,kk+nn-2
	 nums(1,k)=nums(2,k-1)
	 IF (k.EQ.kk+nn-2) THEN
	 nums(2,k)=nums(3,k-(nn-1))
	 nums(3,k)=nums(2,k)-nn
	 nums(4,k)=nums(3,k-1)
	 ELSE
	 nums(2,k)=nums(1,k)+1
	 nums(3,k)=nums(2,k)+nn-2
	 nums(4,k)=nums(3,k)-1
	 ENDIF
	ENDDO	
	kk=kk+nn-1
	ENDDO
	!couche n-1
	nums(1,kk)=(nn-1)*nn
	nums(2,kk)=nums(3,kk-(nn-1))
	nums(3,kk)=nums(1,kk)+2*nn
	nums(4,kk)=nums(1,kk)+nn
	DO k=kk+1,kk+nn-2
	 nums(1,k)=nums(2,k-1)
	IF (k.EQ.kk+nn-2) THEN
	 nums(2,k)=nums(3,k-(nn-1))
	ELSE	 
	 nums(2,k)=nums(1,k)+1
	ENDIF
	 nums(3,k)=nums(3,k-1)+nn
	 nums(4,k)=nums(3,k)-nn
	ENDDO	
********************
********************
	DO i=1,nn-2
	 ndep  = 4*nn*(nn-1)+(nn-2)*(nn-2)+1+(i-1)*(nn-2)
	 nfin  = 4*nn*(nn-1)+(nn-2)*(nn-2)+i*(nn-2)
	 nmoins= (nn-2)*(nn-3)+1
        DO j=ndep,nfin
         SS(1,j) = -SS(1,j-nmoins)
         SS(2,j) = SS(2,j-nmoins)
         SS(3,j) =-SS(3,j-nmoins)
	 nmoins=nmoins+2
        ENDDO	 
     	ENDDO
	!couche1
	kk=5*(nn-1)*(nn-1)+1
	nums(1,kk)=1+nn*(4*(nn-1)-1)
	nums(2,kk)=1
	nums(3,kk)=nums(2,kk)+nn
	nums(4,kk)=4*nn*(nn-1)+(nn-2)*(nn-2)+1
	DO k=kk+1,kk+nn-2
	 nums(1,k)=nums(1,k-1)-nn
	 nums(2,k)=nums(1,k)+nn
	 nums(3,k)=nums(4,k-1)	
	IF (k.EQ.kk+nn-2) THEN
	 nums(4,k)=nums(1,k)-nn
	ELSE
	 nums(4,k)=nums(4,k-1)+1	
	ENDIF
	ENDDO
	!couches internes
	kk=kk+nn-1
	DO i=1,nn-3
	nums(1,kk)=nums(4,kk-(nn-1))
	nums(2,kk)=nums(3,kk-(nn-1))
	nums(3,kk)=nums(2,kk)+nn
	nums(4,kk)=nums(1,kk)+nn-2
	DO k=kk+1,kk+nn-2
	 IF (k.EQ.kk+nn-2) THEN
	 nums(1,k)=nums(4,k-(nn-1))
	 nums(2,k)=nums(1,k-1)
	 nums(3,k)=nums(2,k)+nn-2
	 nums(4,k)=nums(1,k)-nn
	 ELSE
 	 nums(1,k)=nums(1,k-1)+1
	 nums(2,k)=nums(1,k-1)
	 nums(3,k)=nums(4,k-1)
	 nums(4,k)=nums(3,k)+1
	 ENDIF
	ENDDO	
	kk=kk+nn-1
	ENDDO
	!!couche n-1
	nums(1,kk)=nums(4,kk-(nn-1))
	nums(2,kk)=nums(2,kk-(nn-1))+nn
	nums(3,kk)=nums(2,kk)+nn
	nums(4,kk)=nums(3,kk)+nn
	DO k=kk+1,kk+nn-2
	 nums(1,k)=nums(4,k-(nn-1))
	 nums(2,k)=nums(3,k-(nn-1))
	 nums(3,k)=nums(4,k-1)
	 nums(4,k)=nums(3,k)+nn
	ENDDO	
********************
********************	
	RETURN
	END	
C
C
c=========================================================================
        SUBROUTINE ROTATION (axe,angle,ns,ndep,nfin,nmoins,S)
c
        IMPLICIT NONE
        INTEGER ns,ndep,nfin,nmoins,axe,i
        DOUBLE PRECISION    angle,S(3,2*ns)
c
        IF (axe.EQ.1) THEN
        DO i=ndep,nfin
         S(1,i) = S(1,i-nmoins)
         S(2,i) = Dcos(angle)*S(2,i-nmoins)+Dsin(angle)*S(3,i-nmoins)
         S(3,i) =-Dsin(angle)*S(2,i-nmoins)+Dcos(angle)*S(3,i-nmoins)
        ENDDO
        ELSEIF (axe.EQ.2) THEN
        DO i=ndep,nfin
         S(1,i) = Dcos(angle)*S(1,i-nmoins)+Dsin(angle)*S(3,i-nmoins)
         S(2,i) = S(2,i-nmoins)
         S(3,i) =-Dsin(angle)*S(1,i-nmoins)+Dcos(angle)*S(3,i-nmoins)
        ENDDO
        ELSE
        PRINT*, 'ERREUR DANS LE SOUS-PROGRAMME ROTATION'
        STOP
        ENDIF
c
        RETURN
        END
c=========================================================================
        SUBROUTINE MAILLREF (nn,ns,nc,r,S,num)
        IMPLICIT NONE
C
        INTEGER nn,ns,nc,i,j,ii,k,num(8,nc)
        DOUBLE PRECISION    r,div1,div2,rac,S(3,2*ns),h
C
	h=((2.D0/Dsqrt(3.D0))*r)/(nn-1)
C MAILLAGE de la face (z=1) du cube
        ii=1
          DO j=1,nn
           DO i=1,nn
	    S(1,ii)=-(1.D0/Dsqrt(3.D0))*r+(i-1)*h
	    S(2,ii)=-(1.D0/Dsqrt(3.D0))*r+(j-1)*h
	    S(3,ii)=(1.D0/sqrt(3.D0))*r
            ii=ii+1
           ENDDO
          ENDDO
        k=1
         DO j=1,nn-1
          DO i=1,nn-1
            num(1,k)=i+(j-1)*nn
            num(2,k)=num(1,k)+1
            num(3,k)=num(2,k)+nn
            num(4,k)=num(3,k)-1
            k=k+1
          ENDDO
         ENDDO
C PROJECTION sur L'1/6 de la sphere
        DO i=1,nn*nn
         div1=S(1,i)/S(3,i)
         div2=S(2,i)/S(3,i)
         rac=DSQRT(1.D0+div1**2+div2**2)
         S(1,i)=div1*r/rac
         S(2,i)=div2*r/rac
         S(3,i)=r/rac
        ENDDO
C
        RETURN
        END
c=========================================================================
	SUBROUTINE RENUM (nrecol,nn,nc,ndcouch1,ndcouch2,nfrest,num)
        IMPLICIT NONE
C
        INTEGER nn,nc,nrecol,ndcouch1,ndcouch2,nfrest,k,num(8,nc)
C
	 DO k=ndcouch1,ndcouch1+nn-2
          num(1,k)=num(4,k-(nn-1))
          num(2,k)=num(1,k)+1
          num(3,k)=num(3,k-((nn-1)*(nn-1)))+(nn*(nn-1))
          num(4,k)=num(3,k)-1
         ENDDO
	 DO k=ndcouch1+nn-1,nfrest
            num(1,k)=num(4,k-(nn-1))
            num(2,k)=num(3,k-(nn-1))
            num(3,k)=num(3,k-((nn-1)*(nn-1)))+(nn*(nn-1))
            num(4,k)=num(3,k)-1
         ENDDO
! nrecol=2=========>2 recollements	 
	 IF (nrecol.EQ.2) THEN
          DO k=ndcouch2,ndcouch2+nn-2
           num(1,k)=num(4,k-(nn-1))
           num(2,k)=num(3,k-(nn-1))
           num(3,k)=num(2,k-(nn-1)*(4*(nn-1)-1))
           num(4,k)=num(1,k-(nn-1)*(4*(nn-1)-1))
          ENDDO   	 
	 ELSEIF (nrecol.EQ.1) THEN
	 GOTO 403
	 ENDIF
C
403     RETURN
        END

C======================================================================
      SUBROUTINE MSH_omega (nsse,ncs,epsilon,r,nums,SS)
C
      IMPLICIT NONE
      INTEGER i,j,i1,i2,ncs,nsse,nums(8,ncs)
      DOUBLE PRECISION    r,epsilon,tt,SS(3,2*nsse)
C
C----------------sommets de omega
      i1=nsse+1 ; i2=2*nsse ; tt=1.D0+(epsilon/r)
       DO i=i1,i2  
        DO j=1,3
          SS(j,i)=tt*SS(j,i-nsse)
        ENDDO
       ENDDO
C-----------------cubes de omega
      i1=1 ; i2=ncs
      DO i=i1,i2
          DO j=5,8
           nums(j,i)=nums(j-4,i)+nsse
          ENDDO
      ENDDO
C
      RETURN
      END
C======================================================================
      SUBROUTINE Surface_Gamma (ncs,nsse,nums,SS,AireS)
C
      IMPLICIT NONE
      INTEGER i,ncs,nsse,nums(8,ncs)
      DOUBLE PRECISION xc,yc,zc,xcn,ycn,zcn,SS(3,2*nsse),AireS(*)
C
! calcul des aires des elements de gamma en utilisant le produit vectoriel
! on suppose que le quadrangle est plan :-( !!!!!
      DO i=1,ncs
         xc = (SS(1,nums(1,i))+SS(1,nums(3,i)))*0.5D0
         yc = (SS(2,nums(1,i))+SS(2,nums(3,i)))*0.5D0
         zc = (SS(3,nums(1,i))+SS(3,nums(3,i)))*0.5D0
C--------------
         xcn =(SS(2,nums(2,i))+SS(2,nums(1,i)))*
     >        (SS(3,nums(4,i))+SS(3,nums(1,i)))-
     >        (SS(2,nums(4,i))+SS(2,nums(1,i)))*
     >        (SS(3,nums(2,i))+SS(3,nums(1,i)))
         ycn= (SS(1,nums(4,i))+SS(1,nums(1,i)))*
     >        (SS(3,nums(2,i))+SS(3,nums(1,i)))-
     >        (SS(1,nums(2,i))+SS(1,nums(1,i)))*
     >        (SS(3,nums(4,i))+SS(3,nums(1,i)))
         zcn= (SS(1,nums(2,i))+SS(1,nums(1,i)))*
     >        (SS(2,nums(4,i))+SS(2,nums(1,i)))-
     >        (SS(1,nums(4,i))+SS(1,nums(1,i)))*
     >        (SS(2,nums(2,i))+SS(2,nums(1,i)))
C--------------
      AireS(i)=Dabs(((xc*xcn)+(yc*ycn)+(zc*zcn))/
     >                          Dsqrt(xc**2+yc**2+zc**2))
      ENDDO
C
      RETURN
      END
C======================================================================
      SUBROUTINE Locate_G_Gp(nx,ny,nz,hx,hy,hz,xxi,yyi,zzi,ncs,nss,
     >                    nums,R,SS,nrep)
C======================================================================
	IMPLICIT NONE
C
	INTEGER i,j,k,l,ncs,nss,nx,ny,nz,nrep(*),nums(8,ncs)
	DOUBLE PRECISION hx,hy,hz,xxi,yyi,zzi,SS(3,2*nss)
	DOUBLE PRECISION x1,y1,z1,x3,y3,z3,xc,yc,zc,R,rho
C
      DO l=1,2*ncs
      nrep(l)=0
      ENDDO
CC!Reperage des points de quadrature de Gamma dans le maillage globlal!
      DO l=1,ncs
       x1=SS(1,nums(1,l))
       y1=SS(2,nums(1,l))
       z1=SS(3,nums(1,l))
       x3=SS(1,nums(3,l))
       y3=SS(2,nums(3,l))
       z3=SS(3,nums(3,l))
c
      xc=0.5D0*(x1+x3) ; yc=0.5D0*(y1+y3) ; zc=0.5D0*(z1+z3)
      rho=dsqrt(xc**2+yc**2+zc**2)
      xc=xc*R/rho ; yc=yc*R/rho ; zc=zc*R/rho
c	    
	    i=int(dabs(xc-xxi)/hx)+1
	    j=int(dabs(yc-yyi)/hy)+1
	    k=int(dabs(zc-zzi)/hz)+1
	    nrep(l)=(k-1)*(nx-1)*(ny-1)+(j-1)*(nx-1)+i
      ENDDO
CC!Reperage des points de quadrature de GammaPrime dans le maillage globlal!
      DO l=ncs+1,2*ncs
       x1=SS(1,nums(1,l-ncs)+nss)
       y1=SS(2,nums(1,l-ncs)+nss)
       z1=SS(3,nums(1,l-ncs)+nss)
       x3=SS(1,nums(3,l-ncs)+nss)
       y3=SS(2,nums(3,l-ncs)+nss)
       z3=SS(3,nums(3,l-ncs)+nss)
c
      xc=0.5D0*(x1+x3) ; yc=0.5D0*(y1+y3) ; zc=0.5D0*(z1+z3)
      rho=dsqrt(xc**2+yc**2+zc**2)
      xc=xc*R/rho ; yc=yc*R/rho ; zc=zc*R/rho
c
         i=int(dabs(xc-xxi)/hx)+1
         j=int(dabs(yc-yyi)/hy)+1
         k=int(dabs(zc-zzi)/hz)+1
c       
	    nrep(l)=(k-1)*(nx-1)*(ny-1)+(j-1)*(nx-1)+i
      ENDDO

      RETURN
      END
C======================================================================
      SUBROUTINE LocateSG(nx,ny,nz,hx,hy,hz,xxi,yyi,zzi,nss,
     >                    SS,nreps)
C======================================================================
	IMPLICIT NONE
C
	INTEGER i,j,k,l,nss,nx,ny,nz,nreps(*)
	DOUBLE PRECISION hx,hy,hz,xxi,yyi,zzi,SS(3,2*nss)
	DOUBLE PRECISION x,y,z
C
      DO l=1,nss
      nreps(l)=0
      ENDDO
CC!Reperage des sommets de Gamma dans le maillage globlal!
      DO l=1,nss
       x=SS(1,l)
       y=SS(2,l)
       z=SS(3,l)
c	    
	    i=int(dabs(x-xxi)/hx)+1
	    j=int(dabs(y-yyi)/hy)+1
	    k=int(dabs(z-zzi)/hz)+1
	    nreps(l)=(k-1)*(nx-1)*(ny-1)+(j-1)*(nx-1)+i
      ENDDO
c
	return
	end
C======================================================================
      SUBROUTINE ptsForceLocate(nx,ny,nz,hx,hy,hz,xxi,yyi,zzi,ncs,nss,
     >                    nums,R,SS,ptsForce)
C======================================================================
      IMPLICIT NONE
C
      INTEGER i,j,k,l,ncs,nss,nx,ny,nz,ptsForce(*),nums(8,ncs)
      DOUBLE PRECISION hx,hy,hz,xxi,yyi,zzi,SS(3,2*nss)
      DOUBLE PRECISION x1,y1,z1,x3,y3,z3,xc,yc,zc,R,rho,tt
c
	tt=1.5d0
c
      DO l=1,ncs
      ptsForce(l)=0
      ENDDO
CC!Reperage des points de quadrature de Gamma (gonflée) dans le maillage globlal!
      DO l=1,ncs
       x1=SS(1,nums(1,l))
       y1=SS(2,nums(1,l))
       z1=SS(3,nums(1,l))
       x3=SS(1,nums(3,l))
       y3=SS(2,nums(3,l))
       z3=SS(3,nums(3,l))
c
      xc=0.5D0*(x1+x3) ; yc=0.5D0*(y1+y3) ; zc=0.5D0*(z1+z3)
      rho=dsqrt(xc**2+yc**2+zc**2)
      xc=xc*R*tt/rho ; yc=yc*R*tt/rho ; zc=zc*R*tt/rho
c
         i=int(dabs(xc-xxi)/hx)+1
         j=int(dabs(yc-yyi)/hy)+1
         k=int(dabs(zc-zzi)/hz)+1
         ptsForce(l)=(k-1)*(nx-1)*(ny-1)+(j-1)*(nx-1)+i
      ENDDO
	return
	end
C======================================================================C
C======================================================================
      SUBROUTINE Cubes_Num_B (nmx,nmy,nmz,nc,np,num,xi,yi,zi,hx,hy,hz,R,
     >                        ax,ay,az,ncb,numb)
C
      IMPLICIT NONE
      INTEGER nmx,nmy,nmz,nc,nxx,nyy,nzz,i2,j2,k2,
     >  i,j,k,l,p,np,nsA
      INTEGER num(*),numb(nc,np),ncb(*)
      DOUBLE PRECISION xi,yi,zi,hx,hy,hz,R,ax(*),ay(*),az(*)
C
      DO p=1,np
c      !coin inferieur A1 du + petit  cube contenant la boule d'indice p
        i=int(dabs(ax(p)-xi-R)/hx)+1
        j=int(dabs(ay(p)-yi-R)/hy)+1
        k=int(dabs(az(p)-zi-R)/hz)+1
c      !son numero dans la numerotation globale
        nsA=((k-1)*nmx*nmy)+((j-1)*nmx)+i
c      !point A2 definissant le coin superieur (coin superieur=le sommet de type 7 associe a A2)
        i2=int(dabs(ax(p)-xi+R)/hx)+1
        j2=int(dabs(ay(p)-yi+R)/hy)+1
        k2=int(dabs(az(p)-zi+R)/hz)+1
c      !nbre de pts du cube dans les 3 directions
	nxx=i2+1-i+1 ; nyy=j2+1-j+1 ; nzz=k2+1-k+1
c      !Renumerotation des elements du maillage global formant ce cube
        k=1
        DO l=1,nzz-1
         DO j=1,nyy-1
          DO i=1,nxx-1
           numb(k,p)=nsA+(i-1)+((j-1)*nmx)+((l-1)*nmx*nmy)
	   k=k+1
          ENDDO
         ENDDO
        ENDDO
	!nbre d'elements formant le + petit cube contenant la boule p
        ncb(p)=k-1
      ENDDO
C
      RETURN
      END
C======================================================================
      subroutine area_spherical_quads (ncs,nums,SS,R,area)
C======================================================================
      implicit none
      integer i,k,ncs,nums(8,ncs)
      double precision SS(3,*),area(*),pi,r,r2,angle,prodsca
      double precision o1(3),o2(3),o3(3),o4(3),vect1(3),vect2(3)
c
      pi=dacos(-1.d0) ; r2=r**2
c
      do k=1,ncs
	angle=0.d0
c 
        do i=1,3
         O1(i)=SS(i,nums(1,k))
         O2(i)=SS(i,nums(2,k))
         O3(i)=SS(i,nums(3,k))
         O4(i)=SS(i,nums(4,k))
        enddo
c computing the first spherical angle of the quad
      call n_vect3(o1,o2,vect1)
      call n_vect3(o2,o3,vect2)
      angle=dacos(prodsca(3,vect1,vect2))
c computing the second spherical angle and add it to the first one
      call n_vect3(o2,o3,vect1)
      call n_vect3(o3,o4,vect2)
      angle=angle+dacos(prodsca(3,vect1,vect2))
c computind the third  spherical angle and add it to previous
      call n_vect3(o3,o4,vect1)
      call n_vect3(o4,o1,vect2)
      angle=angle+dacos(prodsca(3,vect1,vect2))
c computing the forth spherical angle and make the sum of the four ones
      call n_vect3(o4,o1,vect1)
      call n_vect3(o1,o2,vect2)
      angle=angle+dacos(prodsca(3,vect1,vect2))
c THE AREA  
      area(k)=R2*dabs(angle-(2.d0*pi))
      enddo
c
      return
      end
C======================================================================
      subroutine normalVectors (ncs,np,nums,SS,R,ax,ay,az,
     >    norx,nory,norz)
C======================================================================
      implicit none
      integer i,k,p,np,ncs,nums(8,ncs)
      double precision SS(3,*),r,xc,yc,zc,rho,nory(*)
      double precision ax(*),ay(*),az(*),norx(*),norz(*),xden
c
	open (40, file="normal.dat")
	DO p=1,np
	   do k=1,ncs
c 
	      xc = (ax(p)+SS(1,nums(1,k))+ax(p)+SS(1,nums(3,k)))*0.5d0
	      yc = (ay(p)+SS(2,nums(1,k))+ay(p)+SS(2,nums(3,k)))*0.5d0
	      zc = (az(p)+SS(3,nums(1,k))+az(p)+SS(3,nums(3,k)))*0.5d0
	      rho=dsqrt(xc**2+yc**2+zc**2)
c       
	      xc = xc * R/rho
	      yc = yc * R/rho
	      zc = zc * R/rho
c
	    !  rho=dsqrt(xc**2+yc**2+zc**2)
c           !
	    !  norx(k) = xc/rho
	    !  nory(k) = yc/rho
	    !  norz(k) = zc/rho
c
	      xden = r**2 - xc**2 - yc**2
	    print*, (2*(r**2)-3*(xc**2)-3*(yc**2) +
     >              (xc**2+yc**2)*(-r**2 - 2*(xc**2) - yc**2)/(r**2) -
     >              2*(xc**2)*(yc**2)/((r**2)*dsqrt(xden))
     >)/xden
	      write(40,*) xc,yc,zc,norx(k),nory(k),norz(k)
	   enddo
	enddo
	close(40)
      return
      end
