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

c Authors: Tuomo Rossi   (tro@math.jyu.fi),
c          Jari Toivanen (tene@math.jyu.fi)

C--------------------------------------------------------------------C
C--------------------------ROUTINES----------------------------------C
C--------------------------------------------------------------------C
c************************************************************************
c
c Subroutine: dcq2d
c
c Purpose:
c   A fast direct method for solving the block tridiagonal
c   linear system
c 
c   Au = f,
c
c   where the matrix A is separable, that is, the system
c   can be expressed in the form:
c
c   (A1(x)M2 + M1(x)A2 + ch*M1(x)M2)u = f.
c
c   The notation "(x)" denotes the tensor product of the matrices,
c   that is, A(x)B = {a_ij*B}.
c   A1 and A2 are symmetric tridiagonal matrices of dimension
c   n1 and n2, respectively. M1 and M2 are tridiagonal matrices of the
c   same dimension. Restriction: the matrix M1 must be positive
c   definite and the matrix A must be nonsingular.
c
c   The above system can be written in a block form
c
c   C_i u_i-1 + D_i u_i + C_i+1 u_i+1 = f_i,                      (1)
c
c   where u_i and f_i are the i:th blocks of length n2 of the vectors
c   u and f, respectively. Here C_i = a1(i)*M2, i=2,...,n1, and
c   D_i = (b1(i) + ch*d1(i))*M2 + d1(i)*A2, i=1,...,n1.
c
c
c Version: 1.0
c
c Date: 29 Jan 1997
c
c Parameters:
c
c   Input:
c          n1     - The dimension of the matrices A1 and M1
c          n2     - The dimension of the matrices A2 and M2
c          ldf    - The leading dimension of the two-dimensional 
c                   array f; ldf should be at least n2
c          a1     - The codiagonal of the matrix A1;
c                   the first components is in position a1(2)
c          b1     - The diagonal of the matrix A1
c          c1     - The codiagonal of the matrix M1;
c                   the first components is in position c1(2)
c          d1     - The diagonal of the matrix M1
c          a2     - The codiagonal of the matrix A2;
c                   the first components is in position a2(2)
c          b2     - The diagonal of the matrix A2
c          c2     - The codiagonal of the matrix M2;
c                   the first components is in position c2(2)
c          d2     - The diagonal of the matrix M2
c          ch     - The coefficient in the equation
c          ldw    - The length of double precision workspace;
c                   the minimum value is 6*nl*n1 + max(9*n1, 10*n2),
c                   where nl = 1 + max(int(log4(n1)), 0)
c          liw    - The length of integer workspace;
c                   the minimum value is
c                   (4**nl - 1)/3 + 2*nl + 5*n1 + 4,
c                   where nl is the same as previously
c          init   - Flag which indicates whether the purpose of
c                   the call is too initialize the data structures
c                   in the workspace (init = .true.) or
c                   to solve the problem (init = .false.);
c                   first the subroutine should be initialized
c                   with init = .true.;
c                   after that the subroutine can be called
c                   several times with init = .false.
c                   provided that the values of n1, n2, ldf,
c                   a1, b1, d1, a2, b2, d2, dw and iw are unchanged
c
c   Input/output:
c          f      - On entry f contains the right hand side vector;
c                   on successful exit f contains the solution;
c                   The components are stored according to the
c                   block representation (1), that is, the first
c                   n2 components of f contain the block f_1 and so on.
c                   The solution is returned in the same order.
c
c   Output:
c          ierr   - Error flag indicating failure.
c                   Possible return values:
c                   0  no error
c                   1  n1 < 1 or n2 < 1 or ldf < n2
c                   2  double precision workspace too short;
c                      the required amount of workspace can
c                      be found in iw(1) if liw > 1
c                   3  integer workspace too short;
c                      the required amount of workspace can
c                      be found in iw(2) if liw > 1
c                   4  failure in the LAPACK subroutine dstebz
c                      while solving eigenvalues;
c                      possibly some of the arrays a1, b1 or d1
c                      is incorrect
c                   5  failure in the LAPACK subroutine dstein
c                      while solving eigenvectors;
c                      possibly some of the arrays a1, b1 or d1
c                      is incorrect
c
c   Workspace:
c          dw     - Double precision workspace, length at least ldw
c          iw     - Integer workspace, length at least liw
c
c
c Subroutines called:
c   dpbstf, dsbgst, dstebz and dstein from the LAPACK library.
c   daxpy, dcopy, dnrm2 and dscal from the BLAS1 library. These subroutines
c   can be obtained from
c   the NETLIB archive.
c
c
c Language: FORTRAN
c
c Portability: FORTRAN-77 with do-enddo extension
c
c
c Algorithm:
c   Divide & conquer algorithm for linear systems with
c   separable block tridiagonal matrices.
c
c Complexity estimate: about 44*n1*n2*log4(n1+1) - 41*n1*n2 flops.
c
c
c References:
c
c   T. Rossi, J. Toivanen:
c   A parallel fast direct solver for block tridiagonal systems with
c   separable matrices of arbitrary dimension,
c   Report 21/96, Laboratory of Scientific Computing,
c   University of Jyvaskyla, 1996.
c
c   T. Rossi, J. Toivanen:
c   DC2D and DC3D Version 1.0 User Manual,
c   Reports on Applied Mathematics and Computing, No. 1,
c   University of Jyvaskyla, 1997.
c
c
c Authors: Tuomo Rossi   (tro@math.jyu.fi),
c          Jari Toivanen (tene@math.jyu.fi)
c
c Address: University of Jyvaskyla
c          Department of Mathematics
c          Laboratory of Scientific Computing
c          P.O. Box 35
c          FIN-40351 Jyvaskyla
c          Finland
c
c Copyright: Tuomo Rossi and Jari Toivanen, 1997
c
c************************************************************************
C========================================================================C
      subroutine dcq2d(n1,n2,f,ldf,a1,b1,c1,d1,a2,b2,c2,d2,ch,
     &                dw,ldw,iw,liw,init,ierr)
C========================================================================C
      integer n1, n2, ldf, ldw, liw, iw(liw), ierr
      double precision f(ldf*n1), a1(n1), b1(n1), c1(n1), d1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2), ch, dw(ldw)
      logical init
c
      integer ieig, iwev, iv1, iv3, ig, ir, ix, itri
      integer ip4, iiwev, isplit
      integer iep, ilb, ilen, iloc, isp, iub, j
      integer k, level, ll, m, nl, ns
      double precision c
c
      ierr = 0
c
      if (n1.lt.1.or.n2.lt.1.or.ldf.lt.n2) then
         ierr = 1
         return
      end if
      if (liw.lt.2) then
         ierr = 3
         return
      end if
c
      nl = 1 + max(int(log(dble(n1))/log(4.d0)),0)
c
c Pointers to the real work space
c
      ieig  = 1
      iwev  = ieig + 6*nl*n1
      iv1   = ieig + 6*nl*n1
      iv3   = iv1  + n2
      ig    = iv3  + n2
      ir    = ig   + 3*n2
      ix    = ir   + 3*n2
      itri  = ix   + n2
      iw(1) = max(iwev+9*n1,itri+n2) - 1
c
      if (iw(1).gt.ldw) then
         ierr = 2
         return
      end if
c
c Pointers to the integer work space
c
      ip4    = 3
      iiwev  = ip4    + nl + 1
      isplit = iiwev  + 5*n1
      iw(2)  = isplit + nl + (4**nl - 1)/3 - 1
c
      if (iw(2).gt.liw) then
         ierr = 3
         return
      end if
c
      if (init) then
         iw(ip4) = 1
         do k=1,nl
            iw(ip4+k) = 4*iw(ip4+k-1)
         end do
c
c Make the division into strips
c
         call inispl(n1,iw(isplit),nl,iw(ip4))
c
c Compute the eigenvalues and eigenvectors for the partial solution problems
c
         do level=1,nl
            m = iw(ip4+level-1)
            do iloc=1,m
               call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
               ilen = iub - ilb - 1
               if (ilen.ge.1) then
                  isp = isplit + level + (iw(ip4+level) - 1)/3
     &                + 4*(iloc - 1)
                  if (level.eq.nl) isp = isplit
                  iep = ieig + 6*(ilb + n1*(level - 1))
                  call eigval(ilen,a1(ilb+1),b1(ilb+1),c1(ilb+1),
     &                        d1(ilb+1),dw(iep),dw(iwev),iw(iiwev),ierr)
                  if (ierr.ne.0) return
                  call eigvec(ilen,a1(ilb+1),b1(ilb+1),c1(ilb+1),
     &                        d1(ilb+1),dw(iep),iw(isp),dw(iwev),
     &                        iw(iiwev),ierr)
                  if (ierr.ne.0) return
               end if
            end do
         end do
c
         return
      end if
c
c First recursion, bottom level
c
      level = nl
      m = iw(ip4+level-1)
c
      do iloc=1,m
c
c Find the bounds for the strip
c
         call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
         ilen = iub - ilb - 1
c
         if (ilen.gt.0) then
            iep = ieig + 6*(ilb + n1*(level - 1))
            ll = ilb*ldf + 1
c
            if (ilen.eq.1) then
c
c Problem with one grid column
c
               c = dw(iep+1)**2
               do k=0,n2-1
                  dw(ix+k) = c*f(ll+k)
               end do
               c = dw(iep) + ch
               call soltri(n2,a2,b2,c,c2,d2,dw(ix),dw(itri))
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(ix),dw(ix),.false.)
            else if (ilen.eq.2) then
c
c Problem with two grid columns
c
               call soldbl(n2,dw(iep),a2,b2,c2,d2,ch,
     &                     f(ll),ldf,dw(ir),n2,dw(itri))
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(ir),dw(ir+n2),.false.)
            else
c
c Problem with three grid columns
c
               call soltrb(n2,dw(iep),a2,b2,c2,d2,ch,
     &                     f(ll),ldf,dw(ir),n2,dw(itri),.false.)
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(ir),dw(ir+2*n2),.false.)
            end if
c
            if (ilb.ne.0) then
               ll = (ilb - 1)*ldf + 1
               call daxpy(n2,1.d0,dw(iv1),1,f(ll),1)
            end if
            if (iub.ne.n1+1) then
               ll = (iub - 1)*ldf + 1
               call daxpy(n2,1.d0,dw(iv3),1,f(ll),1)
            end if
         end if
      end do
c
c First recursion, levels through bottom - 1 to top + 1
c
      do level=nl-1,2,-1
         m = iw(ip4+level-1)
c
         do iloc=1,m
c
c Find the bounds for the strip
c
            call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
            ilen = iub - ilb - 1
            ns = min(ilen/4,3)
            if (ilen.gt.3.and.ns.gt.0) then
c
c Problem with 'ns' grid columns
c
               isp = isplit + level + (iw(ip4+level) - 1)/3
     &             + 4*(iloc - 1) + 1
               do k=0,ns-1
                  ll = (iw(isp+k) - 1)*ldf + 1
                  call dcopy(n2,f(ll),1,dw(ig+k*n2),1)
               end do
               do k=iv1,iv1+2*n2-1
                  dw(k) = 0.d0
               end do
c
c Go through eigenvalues one by one
c
               do j=ilb+1,iub-1
                  iep = ieig + 6*(j - 1 + n1*(level - 1))
c
                  call ftrans(ns,n2,dw(ix),dw(ig),dw(iep+2))
c
                  c = dw(iep) + ch
                  call soltri(n2,a2,b2,c,c2,d2,dw(ix),dw(itri))
c
                  if (ilb.ne.0)
     &               call daxpy(n2,dw(iep+1),dw(ix),1,dw(iv1),1)
                  if (iub.ne.n1+1)
     &               call daxpy(n2,dw(iep+5),dw(ix),1,dw(iv3),1)
               end do
c
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(iv1),dw(iv3),.false.)
c
               if (ilb.ne.0) then
                  ll = (ilb - 1)*ldf + 1
                  call daxpy(n2,1.d0,dw(iv1),1,f(ll),1)
               end if
               if (iub.ne.n1+1) then
                  ll = (iub - 1)*ldf + 1
                  call daxpy(n2,1.d0,dw(iv3),1,f(ll),1)
               end if
            end if
         end do
      end do
c
c Second recursion, levels through top to bottom - 1
c      
      do level=1,nl-1
         m = iw(ip4+level-1)
         do iloc=1,m
c
c Find the bounds for the strip
c
            call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
            ilen = iub - ilb - 1
            ns = min(ilen/4,3)
            if (ilen.gt.3.and.ns.gt.0) then
c
c Problem with 'ns' grid columns
c
               isp = isplit + level + (iw(ip4+level) - 1)/3
     &             + 4*(iloc - 1) + 1
               do k=0,ns-1
                  ll = (iw(isp+k) - 1)*ldf + 1
                  call dcopy(n2,f(ll),1,dw(ig+k*n2),1)
               end do
               if (ilb.ne.0) then
                  ll = (ilb - 1)*ldf + 1
                  call dcopy(n2,f(ll),1,dw(iv1),1)
               end if
               if (iub.ne.n1+1) then
                  ll = (iub - 1)*ldf + 1
                  call dcopy(n2,f(ll),1,dw(iv3),1)
               end if
c
c Set the nonhomogenous boundary conditions
c
               call upforc(n1,n2,ilb,iub,dw(iv1),dw(iv3),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(iv1),dw(iv3),.false.)
c
               do k=ir,ir+ns*n2-1
                  dw(k) = 0.d0
               end do
c
c Go through eigenvalues one by one
c
               do j=ilb+1,iub-1
                  iep = ieig + 6*(j - 1 + n1*(level - 1))
c
                  call ftrans(ns,n2,dw(ix),dw(ig),dw(iep+2))
c
                  if (ilb.ne.0)
     &               call daxpy(n2,dw(iep+1),dw(iv1),1,dw(ix),1)
                  if (iub.ne.n1+1)
     &               call daxpy(n2,dw(iep+5),dw(iv3),1,dw(ix),1)
c
                  c = dw(iep) + ch
                  call soltri(n2,a2,b2,c,c2,d2,dw(ix),dw(itri))
c     
                  do k=0,ns-1
                     call daxpy(n2,dw(iep+k+2),dw(ix),1,dw(ir+k*n2),1)
                  end do
               end do
c
c Update the solution
c
               do k=0,ns-1
                  ll = (iw(isp+k) - 1)*ldf + 1
                  call dcopy(n2,dw(ir+k*n2),1,f(ll),1)
               end do
            end if
         end do
      end do
c
c Second recursion, bottom level
c
      level = nl
      m = iw(ip4+level-1)
c
      do iloc=1,m
c
c Find the bounds for the strip
c
         call getbnd(level,iloc,ilb,iub,iw(isplit),iw(ip4))
         ilen = iub - ilb - 1
c
         if (ilen.gt.0) then
            if (ilb.ne.0) then
               ll = (ilb - 1)*ldf + 1
               call dcopy(n2,f(ll),1,dw(iv1),1)
            end if
            if (iub.ne.n1+1) then
               ll = (iub - 1)*ldf + 1
               call dcopy(n2,f(ll),1,dw(iv3),1)
            end if
c
            iep = ieig + 6*(ilb + n1*(level - 1))
            ll = ilb*ldf + 1
c
            if (ilen.eq.1) then
c
c Problem with one grid column
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(iv1),dw(iv3),.true.)
               c = dw(iep+1)**2
               call dscal(n2,c,f(ll),1)
               c = dw(iep) + ch
               call soltri(n2,a2,b2,c,c2,d2,f(ll),dw(itri))
            else if (ilen.eq.2) then
c
c Problem with two grid columns
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll+ldf),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(iv1),dw(iv3),.true.)
               call soldbl(n2,dw(iep),a2,b2,c2,d2,ch,
     &                     f(ll),ldf,f(ll),ldf,dw(itri))
            else
c
c Problem with three grid columns
c
               call upforc(n1,n2,ilb,iub,f(ll),f(ll+2*ldf),a1,c1,
     &                     a2,b2,c2,d2,ch,dw(iv1),dw(iv3),.true.)
               call soltrb(n2,dw(iep),a2,b2,c2,d2,ch,
     &                     f(ll),ldf,f(ll),ldf,dw(itri),.true.)
            end if
         end if
      end do
c
      return
      end
c
c************************************************************************
c
c Initialization of the data structure containing
c the division into strips
c
C========================================================================C
      subroutine inispl(n,split,level,p4)
C========================================================================C
      integer n, split(*), level, p4(*)
c
      integer i, ipp, icp, iend, ilen, id, im, j, k
c
      split(1) = 0
      split(2) = n + 1
      ipp   = 1
      icp   = 3
      level = 1
c
 100  split(icp) = split(ipp)
      icp = icp + 1
      do i=1,p4(level)
         ipp  = ipp + 1
         iend = split(ipp)
         ilen = iend - split(ipp-1)
         k = min((ilen - 1)/4 + 1, 4)
         id = ilen/k
         im = mod(ilen,k)
         do j=1,4
            k = split(icp-1) + id + min(im,1)
            split(icp) = min(k,iend)
            im = max(im-1,0)
            icp  = icp + 1
         end do
      end do
      ipp   = ipp + 1
      level = level + 1
      if (split(ipp+1)-split(ipp).gt.4) goto 100
c
      return
      end
c
c************************************************************************
c
c Find the bounds for a given strip from the data structure
c
C========================================================================C
      subroutine getbnd(level,loc,ilb,iub,split,p4)
C========================================================================C
      integer level, loc, ilb, iub, split(*), p4(*)
c
      integer i
c
      i = level + (p4(level) - 1)/3 + loc
      ilb = split(i-1)
      iub = split(i)
c
      return
      end
c
c************************************************************************
c
c Compute the eigenvalues for the generalized eigensystem of length n
c
C========================================================================C
      subroutine eigval(n,a,b,c,d,eigen,dw,iw,ierr)
C========================================================================C
      integer n, iw(5*n), ierr
      double precision a(n), b(n), c(n), d(n), eigen(6,n), dw(9*n)
c
      integer i, id, ic, ie, iu, iiw1, iiw2, iiw3, m
      double precision cc, ix
c
c Pointers to the workspace
c
      id = 1
      ic = id + 2*n
      ie = ic + 2*n
      iu = ie + n
c
      iiw1 = 1
      iiw2 = iiw1 + n
      iiw3 = iiw2 + n
c
c     store the matrix M into a working array in uplo='U' form, 
c     and compute the split Cholesky factorization of M.
c
      dw(id) = 0.d0
      call dcopy(n,d,1,dw(id+1),2)
      call dcopy(n-1,c(2),1,dw(id+2),2)
c
      call dpbstf('U',n,1,dw(id),2,ierr)
      if (ierr.ne.0) then
         ierr = 4
         return
      end if
c
c     store the matrix A into a working array in uplo='U' form, 
c     and form matrix C.
c
      dw(ic) = 0.d0
      call dcopy(n,b,1,dw(ic+1),2)
      call dcopy(n-1,a(2),1,dw(ic+2),2)
c
      call dsbgst('N','U',n,1,1,dw(ic),2,dw(id),2,ix,1,
     &            dw(iu),ierr)
      if (ierr.ne.0) then
         ierr = 4
         return
      end if
c
c  Copy the transformed matrix
c
      call dcopy(n,dw(ic+1),2,dw(id),1)
      call dcopy(n-1,dw(ic+2),2,dw(id+n),1)
c
c  Compute the eigenvalues
c
      call dstebz('A','E',n,cc,cc,i,i,0.d0,dw(id),dw(id+n),m,i,
     &            dw(ie),iw(iiw1),iw(iiw2),dw(iu),iw(iiw3),ierr)
c
      if (ierr.ne.0) then
         ierr = 4
         return
      end if
c
      do i=1,n
         eigen(1,i) = dw(ie+i-1)
      end do
c
      return
      end
C
c************************************************************************
c
c Compute the required components of eigenvectors
c for the generalized eigensystem of length n
c
C========================================================================C
      subroutine eigvec(n,a,b,c,d,eigen,isp,dw,iw,ierr)
C========================================================================C
      integer n, isp(*), iw(3*n), ierr
      double precision a(n), b(n), c(n), d(n), eigen(6,*), dw(9*n)
c
      double precision s, mnorma
      integer i, j, k, ipos, iz, idw
c
c Pointers to the workspace
c
      iz  = 1
      idw = iz + n
c
      do j=1,n
         call ginvit(n,1,eigen(1,j),dw(iz),b,a,d,c,dw(idw))
c
c Normalize the eigenvector
c
         s = 1.d0/mnorma(n,d,c,dw(iz),dw(idw))
         do i=1,n
            dw(iz+i-1) = dw(iz+i-1)*s
         end do
c
c Copy the required components
c
         if (n.le.3) then
            do k=1,n
               eigen(k+1,j) = dw(iz+k-1)
            end do
         else
            eigen(2,j) = dw(iz)
            eigen(6,j) = dw(iz+n-1)
            do k=1,3
               ipos = isp(k+1) - isp(1)
               if (ipos.lt.n) then
                  eigen(k+2,j) = dw(iz+ipos-1)
               else
                  eigen(k+2,j) = 0.d0
               end if
            end do
         end if
      end do
c
      return
      end
c************************************************************************
c
c Solve a tridiagonal linear system
c
C========================================================================C
      subroutine soltri(n,a,b,s,c,d,x,w)
C========================================================================C
      integer n
      double precision a(n), b(n), s, c(n), d(n), x(n), w(n)
c
      integer i
      double precision dd, ai, an, xp, wp
c
      dd = 1.d0/(b(1) + s*d(1))
      xp = x(1)*dd
      x(1) = xp
      if (n.eq.1) return
      an = a(2) + s*c(2)
      wp = -an*dd
      w(1) = wp
      do i=2,n-1
         ai = an
         dd = -1.d0/(b(i) + s*d(i) + ai*wp)
         xp = (ai*xp - x(i))*dd
         x(i) = xp
         an = a(i+1) + s*c(i+1)
         wp = an*dd
         w(i) = wp
      end do
      dd = -1.d0/(b(n) + s*d(n) + an*wp)
      xp = (an*xp - x(n))*dd
      x(n) = xp
      do i=n-1,1,-1
         xp = x(i) + w(i)*xp
         x(i) = xp
      end do
c
      return
      end
c
c************************************************************************
c
c Update a force vector
c
C========================================================================C
      subroutine upforc(n1,n2,ilb,iub,fl,fu,a1,c1,a2,b2,c2,d2,ch,
     &                  vl,vu,add)
C========================================================================C
      integer n1, n2, ilb, iub
      double precision fl(n2), fu(n2), a1(n1), c1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision ch, vl(n2), vu(n2)
      logical add
c
      integer k
      double precision cc1, cc2, vl1, vl2, vu1, vu2
c
      if (ilb.ne.0) then
         cc1 = -a1(ilb+1) - ch*c1(ilb+1)
         cc2 = -c1(ilb+1)
         if (add) then
            vl1 = vl(1)
            vl2 = vl(2)
            fl(1) = fl(1) + (cc1*d2(1) + cc2*b2(1))*vl1 +
     &                      (cc1*c2(2) + cc2*a2(2))*vl2
            do k=2,n2-1
               fl(k) = fl(k) + (cc1*c2(k)   + cc2*a2(k))*vl1 +
     &                         (cc1*d2(k)   + cc2*b2(k))*vl2 +
     &                         (cc1*c2(k+1) + cc2*a2(k+1))*vl(k+1)
               vl1 = vl2
               vl2 = vl(k+1)
            end do
            fl(n2) = fl(n2) + (cc1*c2(n2) + cc2*a2(n2))*vl1 +
     &                        (cc1*d2(n2) + cc2*b2(n2))*vl2
         else
            vl1 = vl(1)
            vl2 = vl(2)
            fl(1) = (cc1*d2(1) + cc2*b2(1))*vl1 +
     &              (cc1*c2(2) + cc2*a2(2))*vl2
            do k=2,n2-1
               fl(k) = (cc1*c2(k)   + cc2*a2(k))*vl1 +
     &                 (cc1*d2(k)   + cc2*b2(k))*vl2 +
     &                 (cc1*c2(k+1) + cc2*a2(k+1))*vl(k+1)
               vl1 = vl2
               vl2 = vl(k+1)
            end do
            fl(n2) = (cc1*c2(n2) + cc2*a2(n2))*vl1 +
     &               (cc1*d2(n2) + cc2*b2(n2))*vl2
         end if
      end if
c
      if (iub.ne.n1+1) then
         cc1 = -a1(iub) - ch*c1(iub)
         cc2 = -c1(iub)
         if (add) then
            vu1 = vu(1)
            vu2 = vu(2)
            fu(1) = fu(1) + (cc1*d2(1) + cc2*b2(1))*vu1 +
     &                      (cc1*c2(2) + cc2*a2(2))*vu2
            do k=2,n2-1
               fu(k) = fu(k) + (cc1*c2(k)   + cc2*a2(k))*vu1 +
     &                         (cc1*d2(k)   + cc2*b2(k))*vu2 +
     &                         (cc1*c2(k+1) + cc2*a2(k+1))*vu(k+1)
               vu1 = vu2
               vu2 = vu(k+1)
            end do
            fu(n2) = fu(n2) + (cc1*c2(n2) + cc2*a2(n2))*vu1 +
     &                        (cc1*d2(n2) + cc2*b2(n2))*vu2
         else
            vu1 = vu(1)
            vu2 = vu(2)
            fu(1) = (cc1*d2(1) + cc2*b2(1))*vu1 +
     &              (cc1*c2(2) + cc2*a2(2))*vu2
            do k=2,n2-1
               fu(k) = (cc1*c2(k)   + cc2*a2(k))*vu1 +
     &                 (cc1*d2(k)   + cc2*b2(k))*vu2 +
     &                 (cc1*c2(k+1) + cc2*a2(k+1))*vu(k+1)
               vu1 = vu2
               vu2 = vu(k+1)
            end do
            fu(n2) = (cc1*c2(n2) + cc2*a2(n2))*vu1 +
     &               (cc1*d2(n2) + cc2*b2(n2))*vu2
         end if
      end if
c
      return
      end
c************************************************************************
c
c Solve a coupled problem with 3 columns and n rows
c using separation technique
c
C========================================================================C
      subroutine soltrb(n,eigen,a,b,c,d,ch,f,ldf,u,ldu,w,midcol)
C========================================================================C
      integer n, ldf, ldu
      double precision eigen(6,3), a(n), b(n), c(n), d(n)
      double precision ch, f(ldf,3), u(ldu,3), w(n)
      logical midcol
c
      integer i
      double precision ev11, ev12, ev13, ev21, ev22, ev23
      double precision ev31, ev32, ev33, u1, u2, u3, cc
c
      ev11 = eigen(2,1)
      ev12 = eigen(3,1)
      ev13 = eigen(4,1)
      ev21 = eigen(2,2)
      ev22 = eigen(3,2)
      ev23 = eigen(4,2)
      ev31 = eigen(2,3)
      ev32 = eigen(3,3)
      ev33 = eigen(4,3)
c
c First Fourier transform
c
      do i=1,n
         u1     = f(i,1)
         u2     = f(i,2)
         u3     = f(i,3)
         u(i,1) = ev11*u1 + ev12*u2 + ev13*u3
         u(i,2) = ev21*u1 + ev22*u2 + ev23*u3
         u(i,3) = ev31*u1 + ev32*u2 + ev33*u3
      end do
c
c Solve the tridiagonal systems
c
      do i=1,3
         cc = eigen(1,i) + ch
         call soltri(n,a,b,cc,c,d,u(1,i),w)
      end do
c
c Second Fourier transform
c
      if (midcol) then
         do i=1,n
            u1     = u(i,1)
            u2     = u(i,2)
            u3     = u(i,3)
            u(i,1) = ev11*u1 + ev21*u2 + ev31*u3
            u(i,2) = ev12*u1 + ev22*u2 + ev32*u3
            u(i,3) = ev13*u1 + ev23*u2 + ev33*u3
         end do
      else
         do i=1,n
            u1     = u(i,1)
            u2     = u(i,2)
            u3     = u(i,3)
            u(i,1) = ev11*u1 + ev21*u2 + ev31*u3
            u(i,3) = ev13*u1 + ev23*u2 + ev33*u3
         end do
      end if
c
      return
      end
c************************************************************************
c
c Do Fourier transform for 1, 2 or 3 columns
c
C========================================================================C
      subroutine ftrans(ns,n,x,g,evec)
C========================================================================C
      integer ns, n
      double precision x(n), g(n,ns), evec(ns)
c
      integer k
      double precision e1, e2, e3
c
      if (ns.eq.3) then
         e1 = evec(1)
         e2 = evec(2)
         e3 = evec(3)
         do k=1,n
            x(k) = e1*g(k,1) + e2*g(k,2) + e3*g(k,3)
         end do
      else if (ns.eq.2) then
         e1 = evec(1)
         e2 = evec(2)
         do k=1,n
            x(k) = e1*g(k,1) + e2*g(k,2)
         end do
      else
         e1 = evec(1)
         do k=1,n
            x(k) = e1*g(k,1)
         end do
      end if
c
      return
      end
c************************************************************************
c
c Solve a coupled problem with 2 columns and n rows
c using separation technique
c
C========================================================================C
      subroutine soldbl(n,eigen,a,b,c,d,ch,f,ldf,u,ldu,w)
C========================================================================C
      integer n, ldf, ldu
      double precision eigen(6,2), a(n), b(n), c(n), d(n)
      double precision ch, f(ldf,2), u(ldu,2), w(n)
c
      integer i
      double precision ev11, ev12, ev21, ev22, u1, u2, cc
c
      ev11 = eigen(2,1)
      ev12 = eigen(3,1)
      ev21 = eigen(2,2)
      ev22 = eigen(3,2)
c
c First Fourier transform
c
      do i=1,n
         u1     = f(i,1)
         u2     = f(i,2)
         u(i,1) = ev11*u1 + ev12*u2
         u(i,2) = ev21*u1 + ev22*u2
      end do
c
c Solve the tridiagonal systems
c
      do i=1,2
         cc = eigen(1,i) + ch
         call soltri(n,a,b,cc,c,d,u(1,i),w)
      end do
c
c Second Fourier transform
c
      do i=1,n
         u1     = u(i,1)
         u2     = u(i,2)
         u(i,1) = ev11*u1 + ev21*u2
         u(i,2) = ev12*u1 + ev22*u2
      end do
c
      return
      end
C
C========================================================================C
      subroutine ginvit(n,m,e,v,d,b,dd,bb,w)
C========================================================================C
      integer n, m
      double precision e(m), v(n,m), d(n), b(n), dd(n), bb(n), w(4*n)
c
      double precision c
      double precision s, t, epsmat, dasum
      integer i, j, it
c
      t = sqrt(dble(n))*epsmat(n,d,b)
c
c Inverse iteration
c
      do j=1,m
         call formlu(n,d,b,dd,bb,e(j),w)
         c = t
         do i=1,n
            v(i,j) = c
         end do
         call backlu(n,v(1,j),w)
         it = 0
c
 100     it = it + 1
         call forwlu(n,v(1,j),b,bb,e(j),w)
         call backlu(n,v(1,j),w)
c
         s = dasum(n,v(1,j),1)
         if (it.le.5.and.s.lt.1.d0) then
            c = t/s
            call dscal(n,c,v(1,j),1)
            goto 100
         end if
c
         c = 1.d0/s
         call dscal(n,c,v(1,j),1)
      end do
c
      return
      end
C
C========================================================================C
      subroutine formlu(n,d,b,dd,bb,e,w)
C========================================================================C
      integer n
      double precision d(n), b(n), dd(n), bb(n), e, w(4,n)
c
      double precision u, v, xu
      double precision epsmat
      integer i
c
      u = d(1) - e*dd(1)
      if (n.gt.1) v = b(2) - e*bb(2)
c
      do i=2,n
         if (abs(b(i)-e*bb(i)).lt.abs(u)) then
            xu = (b(i)-e*bb(i))/u
            w(1,i) = xu
            w(2,i-1) = u
            w(3,i-1) = v
            w(4,i-1) = 0.d0
            u = d(i) - e*dd(i) - xu*v
            if (i.ne.n) v = b(i+1) - e*bb(i+1)
         else
            xu = u/(b(i) - e*bb(i))
            w(1,i) = xu
            w(2,i-1) = b(i) - e*bb(i)
            w(3,i-1) = d(i) - e*dd(i)
            if (i.ne.n) then
               w(4,i-1) = b(i+1) - e*bb(i+1)
            else
               w(4,i-1) = 0.d0
            end if
            u = v - xu*w(3,i-1)
            v = -xu*w(4,i-1)
         end if
      end do
c
      if (u.eq.0.d0) then
         w(2,n) = epsmat(n,d,b)
      else
         w(2,n) = u
      end if
c
      return
      end
c
c
c
C========================================================================C
      subroutine forwlu(n,x,b,bb,e,w)
C========================================================================C
      integer n
      double precision x(n), b(n), bb(n), e, w(4,n)
c
      double precision u
      integer i
c
      do i=2,n
         u = x(i)
         if (w(2,i-1).eq.b(i)-e*bb(i)) then
            u = x(i-1)
            x(i-1) = x(i)
         end if
         x(i) = u - w(1,i)*x(i-1)
      end do
c
      return
      end
c
c
c
C========================================================================C
      subroutine backlu(n,x,w)
C========================================================================C
      integer n
      double precision x(n), w(4,n)
c
      integer i
c
      x(n) = x(n)/w(2,n)
      if (n.gt.1) x(n-1) = (x(n-1) - w(3,n-1)*x(n))/w(2,n-1)
c
      do i=n-2,1,-1
         x(i) = (x(i) - w(3,i)*x(i+1) - w(4,i)*x(i+2))/w(2,i)
      end do
c
      return
      end
C
C--------------------------------------------------------------------C
C--------------------------FUNCTIONS---------------------------------C
C--------------------------------------------------------------------C
C========================================================================C
      double precision function mnorma(n,d,c,v,dw)
C========================================================================C
      integer n
      double precision d(n), c(n), v(n), dw(n)
c
      integer i
      double precision ddot
c
      dw(1) = d(1)*v(1) + c(2)*v(2)
      do i=2,n-1
         dw(i) = d(i)*v(i) + c(i)*v(i-1) + c(i+1)*v(i+1)
      end do
      dw(n) = d(n)*v(n) + c(n)*v(n-1)
c
      mnorma = sqrt(ddot(n,v,1,dw,1))
c
      return
      end
c
c
C========================================================================C
      double precision function epsmat(n,d,b)
C========================================================================C
      integer n
      double precision d(n), b(n)
c
      double precision a, c, norm
      integer i
c
      a = 4.d0/3.d0
      c = a - 1.d0
      a = c + c + c
      epsmat = abs(a - 1.d0)
      if (epsmat.eq.0.d0) stop 'function epsmat failed'
c
      norm = abs(d(1))
      do i=2,n
         norm = max(abs(d(i)) + abs(b(i)), norm)
      end do
c
      epsmat = epsmat*norm
c
      return
      end
c

