c************************************************************************
c
c Subroutine: dcq3d
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
c   (A1(x)M2(x)M3 + M1(x)A2(x)M3 + M1(x)M2(x)A3 + ch*M1(x)M2(x)M3)u = f.
c
c   The notation "(x)" denotes the tensor product of the matrices,
c   that is, A(x)B = {a_ij*B}.
c   A1, A2 and A3 are symmetric tridiagonal matrices of dimension
c   n1, n2 and n3, respectively. M1, M2 and M3 are tridiagonal matrices of the
c   same dimension. Restriction: the matrices M1 and M2 must be positive
c   definite and the matrix A must be nonsingular.
c
c   The above system can be written in a block form
c
c   C_i u_i-1 + D_i u_i + C_i+1 u_i+1 = f_i,                      (1)
c
c   where u_i and f_i are the i:th blocks of length n2*n3 of the vectors
c   u and f, respectively. Here C_i = (a1(i) + ch*c(i))*M2(x)M3, i=2,...,n1,
c   and D_i = (b1(i) + ch*d1(i))*M2(x)M3 + d1(i)*(A2(x)M3+M2(x)A3), i=1,...,n1.
c
c Version: 0.2
c
c Date: 29 Jan 1997
c
c Parameters:
c
c   Input:
c          n1     - The dimension of the matrices A1 and M1
c          n2     - The dimension of the matrices A2 and M2
c          n3     - The dimension of the matrices A3 and M3
c          ldf3   - The first leading dimension of the three-dimensional 
c                   array f; ldf3 should be at least n3. The right-hand
c                   side vector f is assumed to be stored into an
c                   array f(ldf3,ldf2,*) in the order 3-2-1.
c          ldf2   - The second leading dimension of the three-dimensional 
c                   array f; ldf2 should be at least n2
c          a1     - The codiagonal of the matrix A1;
c                   the first component is in position a1(2)
c          b1     - The diagonal of the matrix A1
c          c1     - The codiagonal of the matrix M1;
c                   the first component is in position c1(2)
c          d1     - The diagonal of the matrix M1
c          a2     - The codiagonal of the matrix A2;
c                   the first component is in position a2(2)
c          b2     - The diagonal of the matrix A2
c          c2     - The codiagonal of the matrix M2;
c                   the first component is in position c2(2)
c          d2     - The diagonal of the matrix M2
c          a3     - The codiagonal of the matrix A3;
c                   the first component is in position a3(2)
c          b3     - The diagonal of the matrix A3
c          c3     - The codiagonal of the matrix M3;
c                   the first component is in position c3(2)
c          d3     - The diagonal of the matrix M3
c          ch     - The coefficient in the equation
c          ldw    - The length of double precision workspace;
c                   the minimum value is
c                   6*(n1*nl1 + n2*nl2) + max(9*n1,7*n2*n3) +
c                   max(9*n2,10*n3),
c                   where
c                   nl1 = 1 + max(int(log(dble(n1))/log(4.d0)),0) and
c                   nl2 = 1 + max(int(log(dble(n2))/log(4.d0)),0)
c          liw    - The length of integer workspace;
c                   the minimum value is
c                   (4**nl1 + 4**nl2 - 2)/3 + 2*(nl1 + nl2) + 5*n1 + 5*n2 + 7,
c                   where nl1 and nl2 are the same as previously
c          init   - Flag which indicates whether the purpose of
c                   the call is to initialize the data structures
c                   in the workspace (init = .true.) or
c                   to solve the problem (init = .false.);
c                   first the subroutine should be initialized
c                   with init = .true.;
c                   after that the subroutine can be called
c                   several times with init = .false.
c                   provided that the values of n1, n2, n3, ldf2, ldf3,
c                   a1, b1, d1, a2, b2, d2, a3, b3, d3, dw and iw
c                   are unchanged
c
c   Input/output:
c          f      - On entry f contains the right hand side vector;
c                   on successful exit f contains the solution;
c                   The components are stored according to the
c                   block representation (1), that is, the first
c                   n2*ldf3 components of f contain the block f_1 and so on.
c                   The solution is returned in the same order.
c
c   Output:
c          ierr   - Error flag indicating failure.
c                   Possible return values:
c                   0   no error
c                   1   n1 < 1 or n2 < 1 or n3 < 1 or ldf2 < n2
c                       or ldf3 < n3
c                   2   double precision workspace too short;
c                       the required amount of workspace can
c                       be found in iw(1) if liw >= 2
c                   3   integer workspace too short;
c                       the required amount of workspace can
c                       be found in iw(2) if liw >= 2
c                   4   failure in the LAPACK subroutine dstebz
c                       while solving eigenvalues;
c                       possibly some of the arrays a1, b1 or d1
c                       is incorrect
c                   5   failure in the LAPACK subroutine dstein
c                       while solving eigenvectors;
c                       possibly some of the arrays a1, b1 or d1
c                       is incorrect
c                   6-8 internal error; this should never happen
c                   9   failure in the LAPACK subroutine dstebz
c                       while solving eigenvalues of two-dimensional
c                       problems; possibly some of the arrays a2, b2 or d2
c                       is incorrect
c                   10  failure in the LAPACK subroutine dstein
c                       while solving eigenvectors of two-dimensional
c                       problems; possibly some of the arrays a2, b2 or d2
c                       is incorrect 
c   Workspace:
c          dw     - Double precision workspace, length at least ldw
c          iw     - Integer workspace, length at least liw
c
c
c Subroutines called:
c   inispl, getbnd, eigval, eigvec and dcq2d from dcq2d.f
c   daxpy and dscal from the BLAS1 library. This library
c   can be obtained from the NETLIB archive.
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
c Complexity estimate: about 
c
c           104*n1*n2*n3*log4 ((n1 + 1)/2)*log4 ((n2 + 1)/2)
c
c flops.
c
c
c References:
c   A.A. Abakumov, A.Yu. Yeremin, Yu.A. Kuznetsov:
c   Efficient fast direct method of solving Poisson's equation on
c   a parallelepiped and its implementation in array processor.
c   Sov. J. Numer. Anal. Math. Modelling 3(1988), 1-20.
c
c   T. Rossi:
c   Fictitious domain methods with separable preconditioners.
c   Ph.D. thesis, Report 69, University of Jyvaskyla,
c   Department of Mathematics, Jyvaskyla, 1995.
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
c************************************************************************
      subroutine dcq3d(n1,n2,n3,f,ldf2,ldf3,a1,b1,c1,d1,a2,b2,c2,d2,
     &                a3,b3,c3,d3,ch,dw,ldw,iw,liw,init,ierr)
      integer n1, n2, n3, ldf2, ldf3, ldw, liw, iw(liw), ierr
      double precision f(ldf2*ldf3*n1), a1(n1), b1(n1), c1(n1), d1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision ch, dw(ldw)
      logical init
c
      integer i, idw2d, ieig, iiwev, iiw2d, ilb, ilen, iloc, ip4
      integer ir, irk, isp, isplit, iub, iv1, iv3, iv5, iwev, ix, j
      integer iep, k, ldw2d, level, liw2d, ll, m, nl, nl2, ns
      double precision c
c
      ierr = 0
c
      if (n1.lt.1.or.n2.lt.1.or.n3.lt.1.or.ldf2.lt.n2.or.
     &    ldf3.lt.n3) then
         ierr = 1
         return
      end if
      if (liw.lt.2) then
         ierr = 3
         return
      end if
c
      nl  = 1 + max(int(log(dble(n1))/log(4.d0)),0)
      nl2 = 1 + max(int(log(dble(n2))/log(4.d0)),0)
c
c Pointers to the real work space
c
      ldw2d = 6*n2*nl2 + max(9*n2, 10*n3)
c
      ieig  = 1
      iwev  = ieig  + 6*n1*nl
      iv1   = ieig  + 6*n1*nl
      iv3   = iv1   + n2*n3
      iv5   = iv3   + n2*n3
      ir    = iv5   + n2*n3
      ix    = ir    + 3*n2*n3
      idw2d = ix    + n2*n3
      iw(1) = idw2d + ldw2d - 1
c
      if (iw(1).gt.ldw) then
         ierr = 2
         return
      end if
c
c Pointers to the integer work space
c
      liw2d  = (4**nl2 - 1)/3 + 2*nl2 + 5*n2 + 4
      ip4    = 3
      iiwev  = ip4    + nl + 1
      isplit = iiwev  + 5*n1
      iiw2d  = isplit + nl + (4**nl - 1)/3
      iw(2)  = iiw2d  + liw2d - 1

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
     &                  + 4*(iloc - 1)
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
         call dcq2d(n2,n3,f,ldf3,a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              dw(idw2d),ldw2d,iw(iiw2d),liw2d,.true.,ierr)
         if (ierr.ne.0) then
            ierr = ierr + 5
         end if
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
            ll = ilb*ldf2*ldf3 + 1
c
            if (ilen.eq.1) then
c
c Problem with one grid plane
c
               c  = dw(iep+1)**2
               do j=0,n2-1
                  do k=0,n3-1
                     dw(ix+j*n3+k) = c*f(ll+j*ldf3+k)
                  end do
               end do
c
               c = dw(iep) + ch
               call dcq2d(n2,n3,dw(ix),n3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &              dw(idw2d),ldw2d,iw(iiw2d),liw2d,.false.,ierr)
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,a1,c1,
     &              a2,b2,c2,d2,a3,b3,c3,d3,ch,dw,n3,ix-1,ix-1,.true.)
c
            else if (ilen.eq.2) then
c
c Problem with two grid planes
c
               call slcp3d(n2,n3,dw(iep),a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f(ilb*ldf2*ldf3+1),ldf2,ldf3,dw(iv1),n2,n3,
     &              dw(idw2d),ldw2d,iw(iiw2d),liw2d)
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,a1,c1,a2,b2,
     &              c2,d2,a3,b3,c3,d3,ch,dw,n3,iv1-1,iv3-1,.true.)
c
c Problem with three grid planes
c
            else
               call sltr3d(n2,n3,dw(iep),a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f(ilb*ldf2*ldf3+1),ldf2,ldf3,dw(iv1),n2,n3,
     &              dw(idw2d),ldw2d,iw(iiw2d),liw2d,.false.)
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,a1,c1,a2,b2,
     &              c2,d2,a3,b3,c3,d3,ch,dw,n3,iv1-1,iv5-1,.true.)
            end if
         end if
      end do
c
c First recursion, levels through bottom - 1 to top + 1
c
      do level=nl-1,2,-1
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
     &               + 4*(iloc - 1)
c
               call dscal(n2*n3,0.d0,dw(iv1),1)
               call dscal(n2*n3,0.d0,dw(iv3),1)
c
c Go through interior columns one by one
c
               do j=ilb+1,iub-1
                  iep = ieig + 6*(j - 1 + n1*(level - 1))
c
                  call ftrn3d(ns,n1,n2,n3,dw(ix),f,ldf2,ldf3,iw(isp+1),
     &                        dw(iep+2))
c     
                  c = dw(iep) + ch
                  call dcq2d(n2,n3,dw(ix),n3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &                      dw(idw2d),ldw2d,iw(iiw2d),liw2d,.false.,
     &                      ierr)
c
                  if (ilb.ne.0) 
     &                 call daxpy(n2*n3,dw(iep+1),dw(ix),1,dw(iv1),1)
                  if (iub.ne.n1+1)
     &                 call daxpy(n2*n3,dw(iep+5),dw(ix),1,dw(iv3),1)
               end do
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,a1,c1,a2,b2,
     &              c2,d2,a3,b3,c3,d3,ch,dw,n3,iv1-1,iv3-1,.true.)
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
c Problem with 'ns' grid planes
c
               isp = isplit + level + (iw(ip4+level) - 1)/3
     &             + 4*(iloc - 1)
c
c Set the nonhomogenous boundary conditions
c
               call upfr3d(n1,n2,n3,ilb,iub,dw,n2,n3,iv1-1,iv3-1,
     &              a1,c1,a2,b2,c2,d2,a3,b3,c3,d3,ch,f,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,.false.)
c
               call dscal(ns*n2*n3,0.d0,dw(ir),1)
c
c Go through the interior grid planes one by one
c
               do j=ilb+1,iub-1
                  iep = ieig + 6*(j - 1 + n1*(level - 1))
c
                  call ftrn3d(ns,n1,n2,n3,dw(ix),f,ldf2,ldf3,iw(isp+1),
     &                        dw(iep+2))
c
                  if (ilb.ne.0)
     &                 call daxpy(n2*n3,dw(iep+1),dw(iv1),1,dw(ix),1)
                  if (iub.ne.n1+1)
     &                 call daxpy(n2*n3,dw(iep+5),dw(iv3),1,dw(ix),1)
c
                  c = dw(iep) + ch
                  call dcq2d(n2,n3,dw(ix),n3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &                      dw(idw2d),ldw2d,iw(iiw2d),liw2d,.false.,
     &                      ierr)
c
                  do k=1,ns
                     irk = (k - 1)*n2*n3
                     call daxpy(n2*n3,dw(iep+k+1),dw(ix),1,dw(ir+irk),1)
                  end do
               end do
c
c Update the solution
c
               do k=1,ns
                  irk = (k - 1)*n2*n3
                  ll  = (iw(isp+k) - 1)*ldf2*ldf3 + 1
                  do i=1,n2
                     do j=0,n3-1
                        f(ll+(i-1)*ldf3+j) = dw(ir+irk+(i-1)*n3+j)
                     end do
                  end do
               end do
            end if
         end do
      end do
c
c Second recursion, bottom level
c
      level = nl
      m     = iw(ip4+level-1)
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
            ll  = ilb*ldf2*ldf3 + 1
c            
            if (ilen.eq.1) then
c
c Problem with one grid plane
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,ll-1,ll-1,
     &              a1,c1,a2,b2,c2,d2,a3,b3,c3,d3,ch,f,ldf3,
     &              (ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,.true.)
c
               c = dw(iep+1)**2
               do j=1,n2
                  call dscal(n3,c,f(ll+(j-1)*ldf3),1)
               end do
c
               c = dw(iep) + ch
               call dcq2d(n2,n3,f(ll),ldf3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &                    dw(idw2d),ldw2d,iw(iiw2d),liw2d,.false.,ierr)
c
            else if (ilen.eq.2) then
c
c Problem with two grid planes
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,ll-1,
     &              (ilb+1)*ldf2*ldf3,a1,c1,a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f,ldf3,(ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,.true.)
c
               call slcp3d(n2,n3,dw(iep),a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f(ll),ldf2,ldf3,f(ll),ldf2,ldf3,dw(idw2d),ldw2d,
     &              iw(iiw2d),liw2d)
            else if (ilen.eq.3) then
c
c Problem with three grid planes
c
               call upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,ll-1,
     &              (ilb+2)*ldf2*ldf3,a1,c1,a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f,ldf3,(ilb-1)*ldf2*ldf3,(iub-1)*ldf2*ldf3,.true.)
c
               call sltr3d(n2,n3,dw(iep),a2,b2,c2,d2,a3,b3,c3,d3,ch,
     &              f(ll),ldf2,ldf3,f(ll),ldf2,ldf3,dw(idw2d),
     &              ldw2d,iw(iiw2d),liw2d,.true.)
            end if
         end if
      end do
c
      return
      end
c
c************************************************************************
c
c Solve a coupled problem with 2 grid planes
c using the separation technique
c
      subroutine slcp3d(n2,n3,eigen,a2,b2,c2,d2,a3,b3,c3,d3,ch,f,
     &                  ldf2,ldf3,u,ldu2,ldu3,dw,ldw,iw,liw)
      integer n2, n3, ldf2, ldf3, ldu2, ldu3, ldw, liw, iw(liw)
      double precision eigen(6,2), a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision ch, f(ldf3,ldf2,2), u(ldu3,ldu2,2), dw(ldw)
c
      integer j, k, ierr
      double precision ev11, ev12, ev21, ev22, u1, u2, c
c
      ev11 = eigen(2,1)
      ev12 = eigen(3,1)
      ev21 = eigen(2,2)
      ev22 = eigen(3,2)
c
c First Fourier transform
c
      do j=1,n2
         do k=1,n3
            u1 = f(k,j,1)
            u2 = f(k,j,2)
            u(k,j,1) = ev11*u1 + ev12*u2
            u(k,j,2) = ev21*u1 + ev22*u2
         end do
      end do
c
c Solve the block tridiagonal systems
c
      do j=1,2
         c = eigen(1,j) + ch
         call dcq2d(n2,n3,u(1,1,j),ldu3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &             dw,ldw,iw,liw,.false.,ierr)
      end do
c
c Second Fourier transform
c
      do j=1,n2
         do k=1,n3
            u1 = u(k,j,1)
            u2 = u(k,j,2)
            u(k,j,1) = ev11*u1 + ev21*u2
            u(k,j,2) = ev12*u1 + ev22*u2
         end do
      end do
c
      return
      end
c
c************************************************************************
c
c Solve a coupled problem with 3 columns and n rows
c using separation technique
c
      subroutine sltr3d(n2,n3,eigen,a2,b2,c2,d2,a3,b3,c3,d3,ch,f,
     &                  ldf2,ldf3,u,ldu2,ldu3,dw,ldw,iw,liw,midcol)
      integer n2, n3, ldf2, ldf3, ldu2, ldu3, ldw, liw, iw(liw)
      double precision eigen(6,3), a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3)
      double precision ch, f(ldf3,ldf2,3), u(ldu3,ldu2,3), dw(ldw)
      logical midcol
c
      integer j, k, ierr
      double precision ev11, ev12, ev13, ev21, ev22, ev23
      double precision ev31, ev32, ev33, u1, u2, u3, c
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
      do j=1,n2
         do k=1,n3
            u1 = f(k,j,1)
            u2 = f(k,j,2)
            u3 = f(k,j,3)
            u(k,j,1) = ev11*u1 + ev12*u2 + ev13*u3
            u(k,j,2) = ev21*u1 + ev22*u2 + ev23*u3
            u(k,j,3) = ev31*u1 + ev32*u2 + ev33*u3
         end do
      end do
c
c Solve the block tridiagonal systems
c
      do j=1,3
         c = eigen(1,j) + ch
         call dcq2d(n2,n3,u(1,1,j),ldu3,a2,b2,c2,d2,a3,b3,c3,d3,c,
     &             dw,ldw,iw,liw,.false.,ierr)
      end do
c
c Second Fourier transform
c
      if (midcol) then
         do j=1,n2
            do k=1,n3
               u1 = u(k,j,1)
               u2 = u(k,j,2)
               u3 = u(k,j,3)
               u(k,j,1) = ev11*u1 + ev21*u2 + ev31*u3
               u(k,j,2) = ev12*u1 + ev22*u2 + ev32*u3
               u(k,j,3) = ev13*u1 + ev23*u2 + ev33*u3
            end do
         end do
      else
         do j=1,n2
            do k=1,n3
               u1 = u(k,j,1)
               u2 = u(k,j,2)
               u3 = u(k,j,3)
               u(k,j,1) = ev11*u1 + ev21*u2 + ev31*u3
               u(k,j,3) = ev13*u1 + ev23*u2 + ev33*u3
            end do
         end do
      end if
c
      return
      end
c
c************************************************************************
c
c Do Fourier transform for 1, 2 or 3 columns
c
      subroutine ftrn3d(ns,n1,n2,n3,x,f,ldf2,ldf3,isplit,evec)
      integer ns, n1, n2, n3, ldf2, ldf3, isplit(3)
      double precision x(n2*n3), f(ldf2*ldf3*n1), evec(ns)
c
      integer j,k, l1, l2, l3
      double precision e1, e2, e3
c
      if (ns.eq.3) then
         l1 = (isplit(1) - 1)*ldf2*ldf3
         l2 = (isplit(2) - 1)*ldf2*ldf3
         l3 = (isplit(3) - 1)*ldf2*ldf3
         e1 = evec(1)
         e2 = evec(2)
         e3 = evec(3)
         do j=1,n2
            do k=1,n3
               x((j-1)*n3+k) = e1*f(l1+(j-1)*ldf3+k) +
     &              e2*f(l2+(j-1)*ldf3+k) + e3*f(l3+(j-1)*ldf3+k)
            end do
         end do
      else if (ns.eq.2) then
         l1 = (isplit(1) - 1)*ldf2*ldf3
         l2 = (isplit(2) - 1)*ldf2*ldf3
         e1 = evec(1)
         e2 = evec(2)
         do j=1,n2
            do k=1,n3
               x((j-1)*n3+k) = e1*f(l1+(j-1)*ldf3+k) +
     &              e2*f(l2+(j-1)*ldf3+k)
            end do
         end do
      else
         l1 = (isplit(1) - 1)*ldf2*ldf3
         e1 = evec(1)
         do j=1,n2
            do k=1,n3
               x((j-1)*n3+k) = e1*f(l1+(j-1)*ldf3+k)
            end do
         end do
      end if
c
      return
      end
c
c************************************************************************
c
c Update a force vector
c
      subroutine upfr3d(n1,n2,n3,ilb,iub,f,ldf2,ldf3,ifl,ifu,a1,c1,
     &                  a2,b2,c2,d2,a3,b3,c3,d3,ch,v,ldv,ivl,ivu,add)
      integer n1, n2, n3, ldf2, ldf3, ilb, iub
      integer ifl, ifu, ldv, ivl, ivu
      double precision f(ldf2*ldf3*n1), a1(n1), c1(n1)
      double precision a2(n2), b2(n2), c2(n2), d2(n2)
      double precision a3(n3), b3(n3), c3(n3), d3(n3), ch, v(*)
      logical add
c
      integer j, k
      double precision c
c
      if (ilb.ne.0) then
         c = -a1(ilb+1)
         if (add) then
            call tmatxv(n2,n3,c2,d2,c3,d3,-a1(ilb+1)-ch*c1(ilb+1),
     &           f(ifl+1),ldf3,v(ivl+1),ldv)
            call tmatxv(n2,n3,a2,b2,c3,d3,-c1(ilb+1),f(ifl+1),ldf3,
     &           v(ivl+1),ldv)
            call tmatxv(n2,n3,c2,d2,a3,b3,-c1(ilb+1),f(ifl+1),ldf3,
     &           v(ivl+1),ldv)
         else
            do j=1,n2
               do k=1,n3
                  f(ifl+(j-1)*ldf3+k) = 0.d0
               end do
            end do
            call tmatxv(n2,n3,c2,d2,c3,d3,-a1(ilb+1)-ch*c1(ilb+1),
     &           f(ifl+1),ldf3,v(ivl+1),ldv)
            call tmatxv(n2,n3,a2,b2,c3,d3,-c1(ilb+1),f(ifl+1),ldf3,
     &           v(ivl+1),ldv)
            call tmatxv(n2,n3,c2,d2,a3,b3,-c1(ilb+1),f(ifl+1),ldf3,
     &           v(ivl+1),ldv)
         end if
      end if
c
      if (iub.ne.n1+1) then
         c = -a1(iub)
         if (add) then
            call tmatxv(n2,n3,c2,d2,c3,d3,-a1(iub)-ch*c1(iub),
     &           f(ifu+1),ldf3,v(ivu+1),ldv)
            call tmatxv(n2,n3,a2,b2,c3,d3,-c1(iub),f(ifu+1),ldf3,
     &           v(ivu+1),ldv)
            call tmatxv(n2,n3,c2,d2,a3,b3,-c1(iub),f(ifu+1),ldf3,
     &           v(ivu+1),ldv)
         else
            do j=1,n2
               do k=1,n3
                  f(ifu+(j-1)*ldf3+k) = 0.d0
               end do
            end do
            call tmatxv(n2,n3,c2,d2,c3,d3,-a1(iub)-ch*c1(iub),
     &           f(ifu+1),ldf3,v(ivu+1),ldv)
            call tmatxv(n2,n3,a2,b2,c3,d3,-c1(iub),f(ifu+1),ldf3,
     &           v(ivu+1),ldv)
            call tmatxv(n2,n3,c2,d2,a3,b3,-c1(iub),f(ifu+1),ldf3,
     &           v(ivu+1),ldv)
         end if
      end if
c
      return
      end
c
c
c
      subroutine tmatxv(n2,n3,c2,d2,c3,d3,ch,u,ldu,f,ldf)
      integer n2, n3, ldu, ldf
      double precision c2(n2), d2(n2), c3(n3), d3(n3)
      double precision u(ldu,n2), f(ldf,n2), ch
c
      integer i, j
      u(1,1) = u(1,1) + ch*(d2(1)*(d3(1)*f(1,1)+c3(2)*f(2,1)) +
     &         c2(2)*(d3(1)*f(1,2)+c3(2)*f(2,2)))
      do j=2,n3-1
         u(j,1) = u(j,1) + ch*(d2(1)*(c3(j)*f(j-1,1)+d3(j)*f(j,1) +
     &        c3(j+1)*f(j+1,1)) + c2(2)*(c3(j)*f(j-1,2) +
     &        d3(j)*f(j,2)+c3(j+1)*f(j+1,2)))
      end do
      u(n3,1) = u(n3,1) + ch*(d2(1)*(c3(n3)*f(n3-1,1)+d3(n3)*f(n3,1)) +
     &     c2(2)*(c3(n3)*f(n3-1,2)+d3(n3)*f(n3,2)))
      do i=2,n2-1
         u(1,i) = u(1,i) + ch*(c2(i)*(d3(1)*f(1,i-1)+c3(2)*f(2,i-1)) +
     &        d2(i)*(d3(1)*f(1,i)+c3(2)*f(2,i)) +
     &        c2(i+1)*(d3(1)*f(1,i+1)+c3(2)*f(2,i+1)))
         do j=2,n3-1
            u(j,i) = u(j,i) + ch*(c2(i)*(c3(j)*f(j-1,i-1) +
     &           d3(j)*f(j,i-1) + c3(j+1)*f(j+1,i-1)) + 
     &           d2(i)*(c3(j)*f(j-1,i)+d3(j)*f(j,i) +
     &           c3(j+1)*f(j+1,i)) + c2(i+1)*(c3(j)*f(j-1,i+1) +
     &           d3(j)*f(j,i+1) + c3(j+1)*f(j+1,i+1)))
         end do
         u(n3,i) = u(n3,i) + ch*(c2(i)*(c3(n3)*f(n3-1,i-1) +
     &        d3(n3)*f(n3,i-1)) + d2(i)*(c3(n3)*f(n3-1,i) +
     &        d3(n3)*f(n3,i)) + c2(i+1)*(c3(n3)*f(n3-1,i+1) +
     &        d3(n3)*f(n3,i+1)))
      end do
      u(1,n2) = u(1,n2) + ch*(c2(n2)*(d3(1)*f(1,n2-1) +
     &     c3(2)*f(2,n2-1)) + d2(n2)*(d3(1)*f(1,n2)+c3(2)*f(2,n2)))
      do j=2,n3-1
         u(j,n2) = u(j,n2) + ch*(c2(n2)*(c3(j)*f(j-1,n2-1) +
     &        d3(j)*f(j,n2-1) + c3(j+1)*f(j+1,n2-1)) + 
     &        d2(n2)*(c3(j)*f(j-1,n2)+d3(j)*f(j,n2) +
     &        c3(j+1)*f(j+1,n2)))
      end do
      u(n3,n2) = u(n3,n2) + ch*(c2(n2)*(c3(n3)*f(n3-1,n2-1) +
     &     d3(n3)*f(n3,n2-1)) + d2(n2)*(c3(n3)*f(n3-1,n2) +
     &     d3(n3)*f(n3,n2)))
c
      return
      end
