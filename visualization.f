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

C=======================================================================C
      subroutine Visu_mtv (nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     > xi,yi,zi,xl,yl,zl,f,U)
C=======================================================================C
C
      implicit none
      integer nmx,nmy,nmz,i,j,k,ii,jj,kk,nx,ny,nz,i_temps
      double precision xi,yi,zi,xl,yl,zl,dt
      double precision f(nmz,nmy,nmx),U(nmz,nmy,nmx)
C
      ii=(nx/2)+1 ; jj=(ny/2)+1 ; kk=(nz/2)+1
C
      print*, 'Writing in file computed_2D_Y.mtv ...'
      OPEN (UNIT=1,FILE='computed_2D_Y.mtv',STATUS='UNKNOWN')
      WRITE(1,*)'$DATA = CONTOUR'
      WRITE(1,*)'%meshplot = off'
      WRITE(1,*)'%contstyle=2'
      WRITE(1,*)'%xmin =', xi,' xmax=', xi+xl
      WRITE(1,*)'%ymin =', zi,' ymax=', zi+zl
      WRITE(1,*)'%nx =  ',nx
      WRITE(1,*)'%ny =  ',nz
      WRITE(1,*)'%nsteps =',50

      do k=1,nz
         do j=ny/2,ny/2
            do i=1,nx
                  write(1,*) f(k,j,i)
            end do
         end do
      end do
      print*, 'done'
      print*, ' '
!      CLOSE(1)
C
      print*, 'Writing in file computed_2D_X.mtv ...'
      OPEN (UNIT=31,FILE='computed_2D_X.mtv',STATUS='UNKNOWN')
      WRITE(31,*)'$DATA = CONTOUR'
      WRITE(31,*)'%meshplot = off'
      WRITE(31,*)'%contstyle=2'
      WRITE(31,*)'%xmin =', yi,' xmax=', yi+yl
      WRITE(31,*)'%ymin =', zi,' ymax=', zi+zl
      WRITE(31,*)'%nx =  ',ny
      WRITE(31,*)'%ny =  ',nz
      WRITE(31,*)'%nsteps =',50

      do k=1,nz
         do j=1,ny
            do i=nx/2,nx/2
                  write(31,*) f(k,j,i)
            end do
         end do
      end do
      print*, 'done'
      print*, ' '

C
      print*, 'Writing in file computed_2D_Z.mtv ...'
      OPEN (UNIT=32,FILE='computed_2D_Z.mtv',STATUS='UNKNOWN')
      WRITE(32,*)'$DATA = CONTOUR'
      WRITE(32,*)'%meshplot = off'
      WRITE(32,*)'%contstyle=2'
      WRITE(32,*)'%xmin =', xi,' xmax=', xi+xl
      WRITE(32,*)'%ymin =', yi,' ymax=', yi+yl
      WRITE(32,*)'%nx =  ',nx
      WRITE(32,*)'%ny =  ',ny
      WRITE(32,*)'%nsteps =',50

      do j=1,ny
         do i=1,nx
            write(32,*) f(nz,j,i)
         end do
      end do
      
      print*, 'done'
      print*, ' '


C
!      print*, 'Writing in file computed_3D.mtv ...'
!      OPEN (UNIT=3,FILE='computed_3D.mtv',STATUS='UNKNOWN')
!        write(3,*)'$DATA = GRID4D'
!        write(3,*)'%toplabel="t=',i_temps*dt,'"'
!        write(3,*)'%axisscale=f'
!        write(3,*)'%meshplot=off'
!	write(3,*)'% vxmin=',xi,' ',' vxmax=',xi+xl
!	write(3,*)'% vymin=',yi,' ',' vymax=',yi+yl
!	write(3,*)'% vzmin=',zi,' ',' vzmax=',zi+zl
!        write(3,*)'% xmin =',xi,' ',' xmax=',xi+xl,' ','nx=',nx
!        write(3,*)'% zmin =',zi,' ',' zmax=',zi+zl,' ','nz=',nz
!        write(3,*)'% ymin =',yi,',',
!     > ' ymax=',  (2*yi+yl)/2  ,' ','ny=',(ny/2)
!        WRITE(3,*)'%nsteps =  ',30
!        do k=1,nz
!         do j=1,(ny/2)
!            do i=1,nx
!                  write(3,*) f(k,j,i)
!            end do
!         end do
!      end do
!      print*, 'done'
!      print*, ' '
!      close(3)
c
!      print*, 'Writing in file error.mtv ...'
!      OPEN (UNIT=2,FILE='error.mtv',STATUS='UNKNOWN')
!      WRITE(2,*)'$DATA = CONTOUR'
!      WRITE(2,*)'%meshplot = on'
!      WRITE(2,*)'%contstyle=2'
!      WRITE(2,*)'%xmin =', xi,' xmax=', xi+xl
!      WRITE(2,*)'%ymin =', yi,' ymax=', yi+yl
!      WRITE(2,*)'%nx =  ',nx
!      WRITE(2,*)'%ny =  ',ny
!      do k=kk,kk
!         do j=1,ny
!            do i=1,nx
!               WRITE(2,*) dabs(f(k,j,i)-u(k,j,i))
!            end do
!         end do
!      end do
!      print*, 'done'
!      print*, ' '
!      CLOSE(2)
C
      return
      end
C=======================================================================C
      subroutine view_sphere (np,ncs,SS,ax,ay,az,nums)
C=======================================================================C
C
      implicit none
      integer np,ncs,p,i,ii,nums(8,ncs)
      double precision SS(3,*),ax(*),ay(*),az(*)
C
       print*, 'Writing in file sphere.mtv ...'
        OPEN(UNIT=9,FILE='sphere.mtv',STATUS='unknown')
        write(9,*) '$DATA = CURVE3D'
        write(9,*) '%dfilltype=1'
        write(9,*) '%pointid=False'
        write(9,*) '%hiddenline=False'
        write(9,*) '%fillcolor=1'
C
      DO p=1,np
	DO i=1,ncs
        WRITE(9,*) ax(p)+SS(1,nums(1,i)), ay(p)+SS(2,nums(1,i)),
     >              az(p)+SS(3,nums(1,i))
        WRITE(9,*) ax(p)+SS(1,nums(2,i)), ay(p)+SS(2,nums(2,i)),
     >              az(p)+SS(3,nums(2,i))
        WRITE(9,*) ax(p)+SS(1,nums(3,i)), ay(p)+SS(2,nums(3,i)),
     >              az(p)+SS(3,nums(3,i))
        WRITE(9,*) ax(p)+SS(1,nums(4,i)), ay(p)+SS(2,nums(4,i)),
     >              az(p)+SS(3,nums(4,i))
        WRITE(9,*)
        ENDDO
      ENDDO
      print*, 'done'
      print*, ' '
        CLOSE(9)
C
      return
      end
C=======================================================================C
      subroutine sphere4D (nss,ncs,SS,ax,ay,az,nums,P)
C=======================================================================C
      implicit none
      integer ncs,nss,i,nums(8,ncs),ii
      double precision SS(3,*),ax(*),ay(*),az(*),P(*)
c
      open(9,FILE="sphere4D.dx", status='unknown')
c
        write(9,*)'object 1 class array type float rank 1' 
        write(9,*)'shape 3 items ',nss,' data follows'

	DO i=1,nss
           WRITE(9,*) ax(1)+SS(1,i), ay(1)+SS(2,i),
     >          az(1)+SS(3,i)
	ENDDO

        write(9,*)'object 2 class array type int rank 1'
	write(9,*)'shape 4 items ',ncs,' data follows'
c

        do i=1,ncs
           WRITE(9,*) nums(1,i)-1,nums(2,i)-1,
     >                nums(4,i)-1,nums(3,i)-1
	enddo
c
        write(9,*)'attribute ','"','element type','"',' string '
     > ,' "quads"'
	write(9,*)'attribute ','"','ref','"',' string ','"','positions',
     >                                                               '"'
	write(9,*)'object 3 class array type float rank 0  items '
	write(9,*)nss,' data follows '
c
	DO i=1,nss
	   WRITE(9,*) P(i)
	ENDDO
c
	write(9,*)'attribute ','"','dep','"',' string ','"','positions',
     >                                                               '"'
        write(9,*)'object ','"','irregular positions', 
     >                       ' irregular connections','"',' class field'
        write(9,*)'component ','"','positions','"',' value 1'
        write(9,*)'component ','"','connections','"',' value 2'
	write(9,*)'component ','"','data','"',' value 3'
	write(9,*)'end'
        print*,'done'
        print*, ' '
c
        close(9)
c
        return
        end

C=======================================================================C
      subroutine sphere4DCenter (nss,ncs,SS,ax,ay,az,nums,P)
C=======================================================================C
      implicit none
      integer ncs,nss,i,nums(8,ncs),ii
      double precision SS(3,*),ax(*),ay(*),az(*),P(*)
c
      open(9,FILE="sphere4DCenter.dx", status='unknown')
c
        write(9,*)'object 1 class array type float rank 1' 
        write(9,*)'shape 3 items ',nss,' data follows'

	DO i=1,nss
           WRITE(9,*) ax(1)+SS(1,i), ay(1)+SS(2,i),
     >          az(1)+SS(3,i)
	ENDDO

        write(9,*)'object 2 class array type int rank 1'
	write(9,*)'shape 4 items ',ncs,' data follows'
c

        do i=1,ncs
           WRITE(9,*) nums(1,i)-1,nums(2,i)-1,
     >                nums(4,i)-1,nums(3,i)-1
	enddo
c
        write(9,*)'attribute ','"','element type','"',' string '
     > ,' "quads"'
	write(9,*)'attribute ','"','ref','"',' string ','"','positions',
     >                                                               '"'
	write(9,*)'object 3 class array type float rank 0  items '
	write(9,*)ncs,' data follows '
c
	DO i=1,ncs
	   WRITE(9,*) P(i)
	ENDDO
c
	write(9,*)'attribute ','"','dep','"',' string ','"','connections'
     >                                                              ,'"'
        write(9,*)'object ','"','irregular positions', 
     >                       ' irregular connections','"',' class field'
        write(9,*)'component ','"','positions','"',' value 1'
        write(9,*)'component ','"','connections','"',' value 2'
	write(9,*)'component ','"','data','"',' value 3'
	write(9,*)'end'
        print*,'done'
        print*, ' '
c
        close(9)
c
        return
        end
C=======================================================================C
      subroutine View_gradient (nx,ny,nz,nmx,nmy,nmz,i_temps,dt,
     > xi,yi,zi,xl,yl,zl,x,y,z,v_X,V_Y,V_Z)
C=======================================================================C
C
      implicit none
      integer nmx,nmy,nmz,i,j,k,ii,jj,kk,nx,ny,nz,i_temps
      double precision xi,yi,zi,xl,yl,zl,dt,x(*),y(*),z(*)
      doubleprecision V_X(nmz,nmy,nmx),V_Y(nmz,nmy,nmx),V_Z(nmz,nmy,nmx)
C
      ii=(nx/2)+1 ; jj=(ny/2)+1 ; kk=(nz/2)+1
C
      print*, 'Writing in file gradient.mtv ...'
      OPEN (UNIT=8,FILE='gradient.mtv',STATUS='UNKNOWN')
C 
      WRITE(8,*)'$DATA = VECTOR'
      write(8,*)'%axisscale=f'
      WRITE(8,*)'%xmin =', xi,' xmax=', xi+xl
      WRITE(8,*)'%ymin =', yi,' ymax=', yi+yl
      WRITE(8,*)'%zmin =', zi,' zmax=', zi+zl
      WRITE(8,*)'%vscale = 0.1'
C
      do k=1,nz!,2
         do j=jj,jj
            do i=1,nx!,2
       write(8,*) x(i),y(j),z(k),V_X(k,j,i),V_Y(k,j,i),V_Z(k,j,i)
            end do
         end do
      end do
      print*, 'done'
      print*, ' '
!      CLOSE(8)
c
      return
      end
C=======================================================================C
      subroutine Visu_DX (l,nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,vfn,
     > x,y,z,uo)
C=======================================================================C
      implicit none
      integer nmx,nmy,nmz,nx,ny,nz,ns,nc,i_temps,i,j,k,l
      integer num(*),numloc(8)
      double precision x(*),y(*),z(*)
      double precision Uo(nmz,nmy,nmx)
      character*4 vfn*11
C
      if     (i_temps.lt.10)then
        write(vfn(11:11),'(I1)')i_temps
      else if(i_temps.lt.100) then
        write(vfn(10:11),'(I2)')i_temps
      else if(i_temps.lt.1000) then
        write(vfn( 9:11),'(I3)')i_temps
      else
        print * , 'i_temps must be in ]0 , 1000[ '
        return
        stop
      endif
C
      print*, 'Writing in file ',vfn, ' ...'
      open(4,FILE=vfn)
      if (l.eq.0) then
        write(4,*)'object 1 class array type float rank 1' 
        write(4,*)'shape 3 items ',ns,' data follows'
	DO k=1,nz
         do j=1,ny
            do i=1,nx
	   WRITE(4,*)x(i),y(j),z(k)
            enddo
         enddo
	ENDDO
        write(4,*)'object 2 class array type int rank 1'
	write(4,*)'shape 8 items ',nc,' data follows'
        DO i=1,nc
	   numloc(1)=num(i)
	   numloc(2)=num(i)+1
           numloc(4)=num(i)+nx+1
	   numloc(3)=num(i)+nx
	   numloc(5)=num(i)+(nx*ny)
	   numloc(6)=num(i)+(nx*ny)+1
           numloc(8)=num(i)+(nx*ny)+1+nx
	   numloc(7)=num(i)+(nx*ny)+nx
        WRITE(4,*) numloc(1)-1,numloc(2)-1,numloc(3)-1,numloc(4)-1,
     >              numloc(5)-1,numloc(6)-1,numloc(7)-1,numloc(8)-1
	enddo
        write(4,*)'attribute ','"','element type','"',' string '
     > ,' "cubes"'
	write(4,*)'attribute ','"','ref','"',' string ','"','positions',
     >                                                               '"'
	write(4,*)'object 3 class array type float rank 0  items '
	write(4,*)ns,' data follows '
	DO k=1,nz
         do j=1,ny
            do i=1,nx
	   WRITE(4,*)uo(k,j,i)!???? indices a l'envers !! voir ligne 226
            enddo
         enddo
	ENDDO
	write(4,*)'attribute ','"','dep','"',' string ','"','positions',
     >                                                               '"'
        write(4,*)'object ','"','irregular positions', 
     >                       ' irregular connections','"',' class field'
        write(4,*)'component ','"','positions','"',' value 1'
        write(4,*)'component ','"','connections','"',' value 2'
	write(4,*)'component ','"','data','"',' value 3'
	write(4,*)'end'
        print*,'done'
        print*, ' '
        else
	DO i=1,nx
         do j=1,ny
            do k=1,nz
	   WRITE(4,*)uo(k,j,i)
            enddo
         enddo
	ENDDO
        endif
        close(4)
C
        return
        end
C=======================================================================C
      subroutine Vof (nx,ny,nz,nmx,nmy,nmz,ns,nc,num,i_temps,
     >vfn,cfn,x,y,z,u1,u3)
C=======================================================================C
      implicit none
      integer nmx,nmy,nmz,nx,ny,nz,ns,nc,i_temps,i,j,k,l,nshow
      integer num(*),numloc(8)
      double precision x(*),y(*),z(*)
      double precision U1(nmz,nmy,nmx),U3(nmz,nmy,nmx)
      character*4 vfn*11,cfn*11
C
      nshow=ny/2
      l=nx*nz
c
      if     (i_temps.lt.10)then
        write(vfn(11:11),'(I1)')i_temps
	   if(i_temps.eq.1)write(cfn(11:11),'(I1)')i_temps
      else if(i_temps.lt.100) then
        write(vfn(10:11),'(I2)')i_temps
	   if(i_temps.eq.1)write(cfn(10:11),'(I2)')i_temps
      else if(i_temps.lt.1000) then
        write(vfn( 9:11),'(I3)')i_temps
	   if(i_temps.eq.1)write(cfn( 9:11),'(I3)')i_temps
      else
        print * , 'i_temps must be in ]0 , 1000[ '
        return
        stop
      endif
C
      print*, 'Writing in files ',vfn,cfn, ' ...'
      open(41,FILE=vfn)
	   if(i_temps.eq.1)open(42,FILE=cfn)
       WRITE(41,*) l
       WRITE(41,*) 0.5
       	   if(i_temps.eq.1)WRITE(42,*) l
      l=1
	DO i=1,nx
         do j=nshow,nshow
            do k=1,nz
	   WRITE(41,1234)l,u1(k,j,i),u3(k,j,i)
	   if(i_temps.eq.1) WRITE(42,1234)l,x(i),z(k)
           l=l+1
            enddo
         enddo
	ENDDO
 1234   FORMAT(I8,D15.7,D15.7)
        close(41)
        close(42)
C
        return
        end


