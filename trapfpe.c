/* ce petit bout code lié à une appli quelconque déclanchera un SIGFPE pour toutes opérations invalides 
 * suivant les critéres définis ci-dessous. (Surtout les NaN !!!)
 *
 * Mes flags en glibc > 2.0, sinon modif plus bas dans __setfpucw */

//#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INEXACT)
//#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)//
#define MY_FLAGS (FE_INVALID | FE_DIVBYZERO)

/* Définition de la version de la glic */
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
  /* Sur intel correspond à 
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
