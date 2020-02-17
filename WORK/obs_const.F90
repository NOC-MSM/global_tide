MODULE obs_const
   !!=====================================================================
   !!                       ***  MODULE obs_const  ***
   !! Observation diagnostics: Constants used by many modules
   !!===================================================================== 
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: obs_const.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

   !! * Modules used
   USE par_kind, ONLY : & ! Precision variables
      & sp         
   IMPLICIT NONE

   !! * Routine/type accessibility
   PUBLIC

   REAL(kind=sp), PARAMETER :: obfillflt=99999.

END MODULE obs_const

