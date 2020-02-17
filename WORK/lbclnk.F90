MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! NEMO        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)  Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment  
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco)  optimization of BDY comm. via lbc_bdy_lnk and lbc_obc_lnk
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie, G. Reffray)  add a C1D case  
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi  
   !!            4.0  ! 2017-03  (G. Madec) automatique allocation of array size (use with any 3rd dim size)
   !!             -   ! 2017-04  (G. Madec) remove duplicated routines (lbc_lnk_2d_9, lbc_lnk_2d_multiple, lbc_lnk_3d_gather)
   !!             -   ! 2017-05  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------
#if defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!           define the generic interfaces of lib_mpp routines
   !!----------------------------------------------------------------------
   !!   lbc_lnk       : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_bdy_lnk   : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean dynamics and tracers   
   USE lib_mpp        ! distributed memory computing library
   USE lbcnfd         ! north fold

   INTERFACE lbc_lnk
      MODULE PROCEDURE   mpp_lnk_2d      , mpp_lnk_3d      , mpp_lnk_4d
   END INTERFACE
   INTERFACE lbc_lnk_ptr
      MODULE PROCEDURE   mpp_lnk_2d_ptr  , mpp_lnk_3d_ptr  , mpp_lnk_4d_ptr
   END INTERFACE
   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE   lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
   END INTERFACE
   !
   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
   END INTERFACE
   !
   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb
   END INTERFACE

   PUBLIC   lbc_lnk       ! ocean/ice lateral boundary conditions
   PUBLIC   lbc_lnk_multi ! modified ocean/ice lateral boundary conditions
   PUBLIC   lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC   lbc_lnk_icb   ! iceberg lateral boundary conditions

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lbclnk.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#else
   !!----------------------------------------------------------------------
   !!   Default option                              shared memory computing
   !!----------------------------------------------------------------------
   !!                routines setting the appropriate values
   !!         on first and last row and column of the global domain
   !!----------------------------------------------------------------------
   !!   lbc_lnk_sum_3d: compute sum over the halos on a 3D variable on ocean mesh
   !!   lbc_lnk_sum_3d: compute sum over the halos on a 2D variable on ocean mesh
   !!   lbc_lnk       : generic interface for lbc_lnk_3d and lbc_lnk_2d
   !!   lbc_lnk_3d    : set the lateral boundary condition on a 3D variable on ocean mesh
   !!   lbc_lnk_2d    : set the lateral boundary condition on a 2D variable on ocean mesh
   !!   lbc_bdy_lnk   : set the lateral BDY boundary condition
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers   
   USE dom_oce        ! ocean space and time domain 
   USE in_out_manager ! I/O manager
   USE lbcnfd         ! north fold

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk
      MODULE PROCEDURE   lbc_lnk_2d      , lbc_lnk_3d      , lbc_lnk_4d
   END INTERFACE
   INTERFACE lbc_lnk_ptr
      MODULE PROCEDURE   lbc_lnk_2d_ptr  , lbc_lnk_3d_ptr  , lbc_lnk_4d_ptr
   END INTERFACE
   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE   lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
   END INTERFACE
   !
   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE lbc_bdy_lnk_2d, lbc_bdy_lnk_3d
   END INTERFACE
   !
   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE lbc_lnk_2d_icb
   END INTERFACE
   
   PUBLIC   lbc_lnk       ! ocean/ice  lateral boundary conditions
   PUBLIC   lbc_lnk_multi ! modified ocean/ice lateral boundary conditions
   PUBLIC   lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC   lbc_lnk_icb   ! iceberg lateral boundary conditions
   
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lbclnk.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!======================================================================
   !!   Default option                           3D shared memory computing
   !!======================================================================
   !!          routines setting land point, or east-west cyclic,
   !!             or north-south cyclic, or north fold values
   !!         on first and last row and column of the global domain
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!                   ***  routine lbc_lnk_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in lbc_lnk_... routines
   !!                ptab   :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat :   nature of array grid-points
   !!                psgn   :   sign used across the north fold boundary
   !!                kfld   :   optional, number of pt3d arrays
   !!                cd_mpp :   optional, fill the overlap area only
   !!                pval   :   optional, background value (used at closed boundaries)
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_LNK           lbc_lnk_2d
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           lbc_lnk_2d_ptr
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_LNK           lbc_lnk_3d
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           lbc_lnk_3d_ptr
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_LNK           lbc_lnk_4d
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           lbc_lnk_4d_ptr
#     include "lbc_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_4d
   
   !!======================================================================
   !!   identical routines in both C1D and shared memory computing
   !!======================================================================

   !!----------------------------------------------------------------------
   !!                   ***  routine lbc_bdy_lnk_(2,3)d  ***
   !!
   !!   wrapper rountine to 'lbc_lnk_3d'. This wrapper is used
   !!   to maintain the same interface with regards to the mpp case
   !!----------------------------------------------------------------------
   
   SUBROUTINE lbc_bdy_lnk_3d( pt3d, cd_type, psgn, ib_bdy )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pt3d      ! 3D array on which the lbc is applied
      CHARACTER(len=1)          , INTENT(in   ) ::   cd_type   ! nature of pt3d grid-points
      REAL(wp)                  , INTENT(in   ) ::   psgn      ! sign used across north fold 
      INTEGER                   , INTENT(in   ) ::   ib_bdy    ! BDY boundary set
      !!----------------------------------------------------------------------
      CALL lbc_lnk_3d( pt3d, cd_type, psgn)
   END SUBROUTINE lbc_bdy_lnk_3d


   SUBROUTINE lbc_bdy_lnk_2d( pt2d, cd_type, psgn, ib_bdy )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d      ! 3D array on which the lbc is applied
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! nature of pt3d grid-points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! sign used across north fold 
      INTEGER                 , INTENT(in   ) ::   ib_bdy    ! BDY boundary set
      !!----------------------------------------------------------------------
      CALL lbc_lnk_2d( pt2d, cd_type, psgn)
   END SUBROUTINE lbc_bdy_lnk_2d


!!gm  This routine should be removed with an optional halos size added in argument of generic routines

   SUBROUTINE lbc_lnk_2d_icb( pt2d, cd_type, psgn, ki, kj )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d      ! 2D array on which the lbc is applied
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! nature of pt3d grid-points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! sign used across north fold 
      INTEGER                 , INTENT(in   ) ::   ki, kj    ! sizes of extra halo (not needed in non-mpp)
      !!----------------------------------------------------------------------
      CALL lbc_lnk_2d( pt2d, cd_type, psgn )
   END SUBROUTINE lbc_lnk_2d_icb
!!gm end

#endif

   !!======================================================================
   !!   identical routines in both distributed and shared memory computing
   !!======================================================================

   !!----------------------------------------------------------------------
   !!                   ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Dummy Argument :
   !!       in    ==>   ptab       ! array to be loaded (2D, 3D or 4D)
   !!                   cd_nat     ! nature of pt2d array grid-points
   !!                   psgn       ! sign used across the north fold boundary
   !!       inout <=>   ptab_ptr   ! array of 2D, 3D or 4D pointers
   !!                   cdna_ptr   ! nature of ptab array grid-points
   !!                   psgn_ptr   ! sign used across the north fold boundary
   !!                   kfld       ! number of elements that has been attributed
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!                  ***   lbc_lnk_(2,3,4)d_multi   ***
   !!                     ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Argument : dummy argument use in lbc_lnk_multi_... routines
   !!
   !!----------------------------------------------------------------------

#  define DIM_2d
#     define ROUTINE_MULTI          lbc_lnk_2d_multi
#     define ROUTINE_LOAD           load_ptr_2d
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_2d


#  define DIM_3d
#     define ROUTINE_MULTI          lbc_lnk_3d_multi
#     define ROUTINE_LOAD           load_ptr_3d
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_3d


#  define DIM_4d
#     define ROUTINE_MULTI          lbc_lnk_4d_multi
#     define ROUTINE_LOAD           load_ptr_4d
#     include "lbc_lnk_multi_generic.h90"
#     undef ROUTINE_MULTI
#     undef ROUTINE_LOAD
#  undef DIM_4d

   !!======================================================================
END MODULE lbclnk

