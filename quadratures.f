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

c========================================================================c
      subroutine Quad_coord(n,xmin,xmax,ymin,ymax,zmin,zmax,
     >                      rt,rt2,w0,w1,w2,dx,dy,dz,a,b,c)
c========================================================================c
c
      implicit none
      integer n
      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision rt,rt2,w0,w1,w2,dx,dy,dz
      double precision a(5),b(5),c(5)
c
      if (n.eq.2) then
c
c------------------------------------------------------------------------------c
c       Coordonnes des 8 points de quadrature
c------------------------------------------------------------------------------c
c
      rt = dsqrt(3.d0)/3.d0
c
      dx = (xmax-xmin)/2
      a(3) = (xmax+xmin)/2
      a(1) =a(3)-dx*rt
      a(2) =a(3)+dx*rt      
C
      dy   = (ymax-ymin)/2
      b(3) = (ymax+ymin)/2
      b(1) =b(3)-dy*rt
      b(2) =b(3)+dy*rt      
c
      dz   = (zmax-zmin)/2
      c(3) = (zmax+zmin)/2
      c(1) =c(3)-dz*rt
      c(2) =c(3)+dz*rt    
c
      elseif (n.eq.3) then
c
c------------------------------------------------------------------------------c
c       Coordonnes des 27 points de quadrature
c------------------------------------------------------------------------------c
c
      rt = dsqrt(15.d0)/5.d0
c
      w0  = 8.d0/9.d0
      w1  = 5.d0/9.d0
c
      dx = (xmax-xmin)/2
      a(2) = (xmax+xmin)/2
      a(1) = a(2)-dx*rt
      a(3) = a(2)+dx*rt      
C
      dy   = (ymax-ymin)/2
      b(2) = (ymax+ymin)/2
      b(1) =  b(2)-dy*rt
      b(3) =  b(2)+dy*rt      
c
      dz   = (zmax-zmin)/2
      c(2) = (zmax+zmin)/2
      c(1) =  c(2)-dz*rt
      c(3) =  c(2)+dz*rt    
c
      elseif (n.eq.5) then
c
c------------------------------------------------------------------------------c
c       Coordonnees des 125 points de quadrature
c------------------------------------------------------------------------------c
c
      rt = dsqrt(245.d0-14.d0*dsqrt(70.d0))/21.d0
      rt2 = dsqrt(245.d0+14.d0*dsqrt(70.d0))/21.d0
c
      w0  = 128.d0/225.d0
      w1  = (322.d0+13.d0*dsqrt(70.d0))/900.d0
      w2  = (322.d0-13.d0*dsqrt(70.d0))/900.d0
c
      dx = (xmax-xmin)/2
      a(3) = (xmax+xmin)/2
      a(4) = a(3)-dx*rt2
      a(1) = a(3)-dx*rt
      a(2) = a(3)+dx*rt      
      a(5) = a(3)+dx*rt2      
C
      dy = (ymax-ymin)/2
      b(3) = (ymax+ymin)/2
      b(4) = b(3)-dy*rt2
      b(1) = b(3)-dy*rt
      b(2) = b(3)+dy*rt      
      b(5) = b(3)+dy*rt2      
C
      dz = (zmax-zmin)/2
      c(3) = (zmax+zmin)/2
      c(4) = c(3)-dz*rt2
      c(1) = c(3)-dz*rt
      c(2) = c(3)+dz*rt      
      c(5) = c(3)+dz*rt2      
c
      else
       print*,' '
       print*,'ABORTING.......'
       print*,'Error:',n,'= quadature points number not implemented'
       print*,'n must be 2 3 or 5'
       print*,' '
       stop
c
      endif
c
      return
      end
c========================================================================c
      Function Integ_QN(n,dx,dy,dz,w0,w1,w2,ga)
c========================================================================c
c
      implicit none
      integer n,i,j,k,l
      double precision w0,w1,w2,dx,dy,dz,integ_QN
      double precision ga(5,5,5)
c
      if (n.eq.2) then
c
c------------------------------------------------------------------------------c
c       Formule de quadrature de Gauss-Legendre a 8 points: 2 points dans 
c        dans chaque direction. Exacte pour les polynomes de degres 3 par
c        rapport a chacune des variables
c        par exple l'integrale par rapport a x es donnee par:
c   integ=dx*(ga(x1,y,z)+ga(x2,y,z))
c------------------------------------------------------------------------------c
c
      integ_QN=
     >dx*(
     >     (dy*(
     >           (dz*(ga(1,1,1)+ga(1,1,2)))
     >                                      +
     >           (dz*(ga(1,2,1)+ga(1,2,2)))
     >            ))
     >               +
     >     (dy*(
     >           (dz*(ga(2,1,1)+ga(2,1,2)))
     >                                      +
     >           (dz*(ga(2,2,1)+ga(2,2,2)))
     >            ))
     >      )
c
      elseif (n.eq.3) then
c
c------------------------------------------------------------------------------c
c       Formule de quadrature de Gauss-Legendre a 27 points: 3 points dans 
c        dans chaque direction. Exacte pour les polynomes de degres 5 par
c        rapport a chacune des variables
c        par exple l'integrale par rapport a x es donnee par:
c   integ=dx*(w0*ga(x0,y,z)+w1*(ga(x1,y,z)+ga(x2,y,z)))
c------------------------------------------------------------------------------c
c
      integ_QN=
     >dx*(w0*(dy*(w0* 
     >              (dz*(w0*(ga(2,2,2))+w1*((ga(2,2,1))+(ga(2,2,3)))))+
     >              w1*(
     >              (dz*(w0*(ga(2,3,2))+w1*((ga(2,3,1))+(ga(2,3,3)))))+
     >              (dz*(w0*(ga(2,1,2))+w1*((ga(2,1,1))+(ga(2,1,3))))) 
     >                   ) )) +
     >    w1*((dy*(w0* 
     >              (dz*(w0*(ga(3,2,2))+w1*((ga(3,2,1))+(ga(3,2,3)))))+
     >              w1*(
     >              (dz*(w0*(ga(3,3,2))+w1*((ga(3,3,1))+(ga(3,3,3)))))+
     >              (dz*(w0*(ga(3,1,2))+w1*((ga(3,1,1))+(ga(3,1,3))))) 
     >                  ) )) +
     >        (dy*(w0* 
     >              (dz*(w0*(ga(1,2,2))+w1*((ga(1,2,1))+(ga(1,2,3)))))+
     >              w1*(
     >              (dz*(w0*(ga(1,3,2))+w1*((ga(1,3,1))+(ga(1,3,3)))))+
     >              (dz*(w0*(ga(1,1,2))+w1*((ga(1,1,1))+(ga(1,1,3))))) 
     >                  ) ))
     >         ))
c
      elseif (n.eq.5) then
c
c------------------------------------------------------------------------------c
c       Formule de quadrature de Gauss-Legendre a 125 points: 5 points dans 
c        dans chaque direction. Exacte pour les polynomes de degres 9 par
c        rapport a chacune des variables
c        par exple l'integrale par rapport a x es donnee par:
c   integ=dx*(w0*ga(x0,y,z)+w1*(ga(x1,y,z)+ga(x2,y,z))+w2*(ga(x4,y,z)+ga(x5,y,z)))
c------------------------------------------------------------------------------c
c
      integ_QN=dx*(w0*
     >              (dy*(w0*
     >(dz*(w0*ga(3,3,3)+w1*(ga(3,3,1)+ga(3,3,2))+
     >                  w2*(ga(3,3,4)+ga(3,3,5))))
     >                      +w1*(
     >(dz*(w0*ga(3,1,3)+w1*(ga(3,1,1)+ga(3,1,2))+
     >                  w2*(ga(3,1,4)+ga(3,1,5))))
     >                           +
     >(dz*(w0*ga(3,2,3)+w1*(ga(3,2,1)+ga(3,2,2))+
     >                  w2*(ga(3,2,4)+ga(3,2,5))))
     >                              )+w2*(
     >(dz*(w0*ga(3,4,3)+w1*(ga(3,4,1)+ga(3,4,2))+
     >                  w2*(ga(3,4,4)+ga(3,4,5))))
     >                                    +
     >(dz*(w0*ga(3,5,3)+w1*(ga(3,5,1)+ga(3,5,2))+
     >                  w2*(ga(3,5,4)+ga(3,5,5))))
     >                                      )))
     >         +
     >          w1*(
     >                (dy*(w0*
     >(dz*(w0*ga(1,3,3)+w1*(ga(1,3,1)+ga(1,3,2))+
     >                  w2*(ga(1,3,4)+ga(1,3,5))))
     >                         +w1*(
     >(dz*(w0*ga(1,1,3)+w1*(ga(1,1,1)+ga(1,1,2))+
     >                  w2*(ga(1,1,4)+ga(1,1,5))))
     >                                +
     >(dz*(w0*ga(1,2,3)+w1*(ga(1,2,1)+ga(1,2,2))+
     >                  w2*(ga(1,2,4)+ga(1,2,5))))
     >                                )+w2*(
     >(dz*(w0*ga(1,4,3)+w1*(ga(1,4,1)+ga(1,4,2))+
     >                  w2*(ga(1,4,4)+ga(1,4,5))))
     >                                        +
     >(dz*(w0*ga(1,5,3)+w1*(ga(1,5,1)+ga(1,5,2))+
     >                  w2*(ga(1,5,4)+ga(1,5,5))))  )))
     >                 +
     >           (dy*(w0*
     >(dz*(w0*ga(2,3,3)+w1*(ga(2,3,1)+ga(2,3,2))+
     >                  w2*(ga(2,3,4)+ga(2,3,5))))
     >                   +w1*(
     >(dz*(w0*ga(2,1,3)+w1*(ga(2,1,1)+ga(2,1,2))+
     >                  w2*(ga(2,1,4)+ga(2,1,5))))
     >                           +
     >(dz*(w0*ga(2,2,3)+w1*(ga(2,2,1)+ga(2,2,2))+
     >                  w2*(ga(2,2,4)+ga(2,2,5))))
     >                          )+w2*(
     >(dz*(w0*ga(2,4,3)+w1*(ga(2,4,1)+ga(2,4,2))+
     >                  w2*(ga(2,4,4)+ga(2,4,5))))
     >                                 +
     >(dz*(w0*ga(2,5,3)+w1*(ga(2,5,1)+ga(2,5,2))+
     >                  w2*(ga(2,5,4)+ga(2,5,5))))
     >                                  )))
     >                 )+w2*(
     >           (dy*(w0*
     >(dz*(w0*ga(4,3,3)+w1*(ga(4,3,1)+ga(4,3,2))+
     >                  w2*(ga(4,3,4)+ga(4,3,5))))
     >                   +w1*(
     >(dz*(w0*ga(4,1,3)+w1*(ga(4,1,1)+ga(4,1,2))+
     >                  w2*(ga(4,1,4)+ga(4,1,5))))
     >                         +
     >(dz*(w0*ga(4,2,3)+w1*(ga(4,2,1)+ga(4,2,2))+
     >                  w2*(ga(4,2,4)+ga(4,2,5))))
     >                         )+w2*(
     >(dz*(w0*ga(4,4,3)+w1*(ga(4,4,1)+ga(4,4,2))+
     >                  w2*(ga(4,4,4)+ga(4,4,5))))
     >                               +
     >(dz*(w0*ga(4,5,3)+w1*(ga(4,5,1)+ga(4,5,2))+
     >                  w2*(ga(4,5,4)+ga(4,5,5))))
     >                                 )))
     >                         +
     >            (dy*(w0*
     >(dz*(w0*ga(5,3,3)+w1*(ga(5,3,1)+ga(5,3,2))+
     >                  w2*(ga(5,3,4)+ga(5,3,5))))
     >                   +w1*(
     >(dz*(w0*ga(5,1,3)+w1*(ga(5,1,1)+ga(5,1,2))+
     >                  w2*(ga(5,1,4)+ga(5,1,5))))
     >                         +
     >(dz*(w0*ga(5,2,3)+w1*(ga(5,2,1)+ga(5,2,2))+
     >                  w2*(ga(5,2,4)+ga(5,2,5))))
     >                           )+w2*(
     >(dz*(w0*ga(5,4,3)+w1*(ga(5,4,1)+ga(5,4,2))+
     >                  w2*(ga(5,4,4)+ga(5,4,5))))
     >                                  +
     >(dz*(w0*ga(5,5,3)+w1*(ga(5,5,1)+ga(5,5,2))+
     >                  w2*(ga(5,5,4)+ga(5,5,5))))
     >                                    )))
     >                                       ))  
      else
       print*,' '
       print*,'ABORTING.......'
       print*,'Error:',n,'= quadature points number not implemented'
       print*,'n must be 2 3 or 5'
       print*,' '
       stop
cc
      endif
c
      return
      end

c========================================================================c
c     integration sur les quads de z=zf pour le clacul du debit
      Function Integ_QNZF(n,dx,dy,dz,w0,w1,w2,ga)
c========================================================================c
c
      implicit none
      integer n,i,j,k,l
      double precision w0,w1,w2,dx,dy,dz,integ_QNZF
      double precision ga(5,5,5)
c     
      if (n.eq.2) then
c     
c------------------------------------------------------------------------------c
c     Formule de quadrature de Gauss-Legendre a 8 points: 2 points dans 
c     dans chaque direction. Exacte pour les polynomes de degres 3 par
c     rapport a chacune des variables
c     par exple l'integrale par rapport a x es donnee par:
c     integ=dx*(ga(x1,y,z)+ga(x2,y,z))
c------------------------------------------------------------------------------c
c
         integ_QNZF=
     >   dx*(
     >        (dy*( ga(1,1,2)+ga(1,2,2) ))
     >                                    +
     >        (dy*( ga(2,1,2)+ga(2,2,2) ))
     >       )         
c     
        else
          print*,' '
          print*,'ABORTING.......'
          print*,'Error:',n,'= quadature points number not implemented'
          print*,'n must be 2'
          print*,' '
          stop
cc
       endif
c
       return
       end
