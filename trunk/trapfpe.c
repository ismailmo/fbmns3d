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
c Author : Pascal Hav�   (have@ann.jussieu.fr)

/* ce petit bout code li� � une appli quelconque d�clanchera un SIGFPE pour toutes op�rations invalides 
 * suivant les crit�res d�finis ci-dessous. (Surtout les NaN !!!)
 *
 * Mes flags en glibc > 2.0, sinon modif plus bas dans __setfpucw */

//#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INEXACT)
//#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)//
#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO)

/* D�finition de la version de la glic */
#include <features.h>

/* fichier contenant les constantes */
#if __GLIBC_MINOR__ > 0
#include <fenv.h>
#else
#include <fpu_control.h>
#endif

static void __attribute__ ((constructor)) trapfpe ()
{
#if __GLIBC_MINOR__ > 0
  /* Version glibc > 2.0 */
  /* Sur intel correspond � 
   *  unsigned short cw;
   *  asm volatile ("fnstcw %0":"=m" (cw)); // Get flags
   *  cw &= ~traps;                         // Enable new flags
   *  asm volatile ("fldcw %0"::"m" (cw));  // Set flags
   */
  feenableexcept(MY_FLAGS);
#else
  /* Version glibc 2.0 */
  __setfpucw (_FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM));
#endif
}
