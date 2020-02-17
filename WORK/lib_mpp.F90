MODULE lib_mpp
   !!======================================================================
   !!                       ***  MODULE  lib_mpp  ***
   !! Ocean numerics:  massively parallel processing library
   !!=====================================================================
   !! History :  OPA  !  1994  (M. Guyon, J. Escobar, M. Imbard)  Original code
   !!            7.0  !  1997  (A.M. Treguier)  SHMEM additions
   !!            8.0  !  1998  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
   !!                 !  1998  (J.M. Molines) Open boundary conditions
   !!   NEMO     1.0  !  2003  (J.M. Molines, G. Madec)  F90, free form
   !!                 !  2003  (J.M. Molines) add mpp_ini_north(_3d,_2d)
   !!             -   !  2004  (R. Bourdalle Badie)  isend option in mpi
   !!                 !  2004  (J.M. Molines) minloc, maxloc
   !!             -   !  2005  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
   !!             -   !  2005  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
   !!             -   !  2005  (R. Benshila, G. Madec)  add extra halo case
   !!             -   !  2008  (R. Benshila) add mpp_ini_ice
   !!            3.2  !  2009  (R. Benshila) SHMEM suppression, north fold in lbc_nfd
   !!            3.2  !  2009  (O. Marti)    add mpp_ini_znl
   !!            4.0  !  2011  (G. Madec)  move ctl_ routines from in_out_manager
   !!            3.5  !  2012  (S.Mocavero, I. Epicoco) Add mpp_lnk_bdy_3d/2d routines to optimize the BDY comm.
   !!            3.5  !  2013  (C. Ethe, G. Madec)  message passing arrays as local variables 
   !!            3.5  !  2013  (S.Mocavero, I.Epicoco - CMCC) north fold optimizations
   !!            3.6  !  2015  (O. Tintó and M. Castrillo - BSC) Added '_multiple' case for 2D lbc and max
   !!            4.0  !  2017  (G. Madec) automatique allocation of array argument (use any 3rd dimension)
   !!             -   !  2017  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ctl_stop      : update momentum and tracer Kz from a tke scheme
   !!   ctl_warn      : initialization, namelist read, and parameters control
   !!   ctl_opn       : Open file and check if required file is available.
   !!   ctl_nam       : Prints informations when an error occurs while reading a namelist
   !!   get_unit      : give the index of an unused logical unit
   !!----------------------------------------------------------------------
#if   defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lib_mpp_alloc : allocate mpp arrays
   !!   mynode        : indentify the processor unit
   !!   mpp_lnk       : interface (defined in lbclnk) for message passing of 2d or 3d arrays (mpp_lnk_2d, mpp_lnk_3d)
   !!   mpp_lnk_icb   : interface for message passing of 2d arrays with extra halo for icebergs (mpp_lnk_2d_icb)
   !!   mpprecv       :
   !!   mppsend       :
   !!   mppscatter    :
   !!   mppgather     :
   !!   mpp_min       : generic interface for mppmin_int , mppmin_a_int , mppmin_real, mppmin_a_real
   !!   mpp_max       : generic interface for mppmax_int , mppmax_a_int , mppmax_real, mppmax_a_real
   !!   mpp_sum       : generic interface for mppsum_int , mppsum_a_int , mppsum_real, mppsum_a_real
   !!   mpp_minloc    :
   !!   mpp_maxloc    :
   !!   mppsync       :
   !!   mppstop       :
   !!   mpp_ini_north : initialisation of north fold
   !!   mpp_lbc_north_icb : alternative to mpp_nfd for extra outer halo with icebergs
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE lbcnfd         ! north fold treatment
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   INTERFACE mpp_nfd
      MODULE PROCEDURE   mpp_nfd_2d      , mpp_nfd_3d      , mpp_nfd_4d
      MODULE PROCEDURE   mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
   END INTERFACE

   ! Interface associated to the mpp_lnk_... routines is defined in lbclnk
   PUBLIC   mpp_lnk_2d      , mpp_lnk_3d      , mpp_lnk_4d
   PUBLIC   mpp_lnk_2d_ptr, mpp_lnk_3d_ptr, mpp_lnk_4d_ptr
   !
!!gm  this should be useless
   PUBLIC   mpp_nfd_2d    , mpp_nfd_3d    , mpp_nfd_4d
   PUBLIC   mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
!!gm end
   !
   PUBLIC   ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam
   PUBLIC   mynode, mppstop, mppsync, mpp_comm_free
   PUBLIC   mpp_ini_north
   PUBLIC   mpp_lnk_2d_icb
   PUBLIC   mpp_lbc_north_icb
   PUBLIC   mpp_min, mpp_max, mpp_sum, mpp_minloc, mpp_maxloc
   PUBLIC   mpp_max_multiple
   PUBLIC   mppscatter, mppgather
   PUBLIC   mpp_ini_ice, mpp_ini_znl
   PUBLIC   mppsize
   PUBLIC   mppsend, mpprecv                          ! needed by TAM and ICB routines
   PUBLIC   mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
   PUBLIC   mpprank
   
   !! * Interfaces
   !! define generic interface for these routine as they are called sometimes
   !! with scalar arguments instead of array arguments, which causes problems
   !! for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_sum
      MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real,   &
         &             mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_minloc
      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
   END INTERFACE
   INTERFACE mpp_maxloc
      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
   END INTERFACE
   INTERFACE mpp_max_multiple
      MODULE PROCEDURE mppmax_real_multiple
   END INTERFACE

   !! ========================= !!
   !!  MPI  variable definition !!
   !! ========================= !!
!$AGRIF_DO_NOT_TREAT
   INCLUDE 'mpif.h'
!$AGRIF_END_DO_NOT_TREAT

   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .TRUE.    !: mpp flag

   INTEGER, PARAMETER         ::   nprocmax = 2**10   ! maximun dimension (required to be a power of 2)

   INTEGER ::   mppsize        ! number of process
   INTEGER ::   mpprank        ! process number  [ 0 - size-1 ]
!$AGRIF_DO_NOT_TREAT
   INTEGER, PUBLIC ::   mpi_comm_oce   ! opa local communicator
!$AGRIF_END_DO_NOT_TREAT

   INTEGER :: MPI_SUMDD

   ! variables used in case of sea-ice
   INTEGER, PUBLIC ::   ncomm_ice       !: communicator made by the processors with sea-ice (public so that it can be freed in icethd)
   INTEGER         ::   ngrp_iworld     !  group ID for the world processors (for rheology)
   INTEGER         ::   ngrp_ice        !  group ID for the ice processors (for rheology)
   INTEGER         ::   ndim_rank_ice   !  number of 'ice' processors
   INTEGER         ::   n_ice_root      !  number (in the comm_ice) of proc 0 in the ice comm
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   nrank_ice     ! dimension ndim_rank_ice

   ! variables used for zonal integration
   INTEGER, PUBLIC ::   ncomm_znl       !: communicator made by the processors on the same zonal average
   LOGICAL, PUBLIC ::   l_znl_root      !: True on the 'left'most processor on the same row
   INTEGER         ::   ngrp_znl        !  group ID for the znl processors
   INTEGER         ::   ndim_rank_znl   !  number of processors on the same zonal average
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   nrank_znl  ! dimension ndim_rank_znl, number of the procs into the same znl domain

   ! North fold condition in mpp_mpi with jpni > 1 (PUBLIC for TAM)
   INTEGER, PUBLIC ::   ngrp_world        !: group ID for the world processors
   INTEGER, PUBLIC ::   ngrp_opa          !: group ID for the opa processors
   INTEGER, PUBLIC ::   ngrp_north        !: group ID for the northern processors (to be fold)
   INTEGER, PUBLIC ::   ncomm_north       !: communicator made by the processors belonging to ngrp_north
   INTEGER, PUBLIC ::   ndim_rank_north   !: number of 'sea' processor in the northern line (can be /= jpni !)
   INTEGER, PUBLIC ::   njmppmax          !: value of njmpp for the processors of the northern line
   INTEGER, PUBLIC ::   north_root        !: number (in the comm_opa) of proc 0 in the northern comm
   INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE ::   nrank_north   !: dimension ndim_rank_north

   ! Type of send : standard, buffered, immediate
   CHARACTER(len=1), PUBLIC ::   cn_mpi_send        !: type od mpi send/recieve (S=standard, B=bsend, I=isend)
   LOGICAL         , PUBLIC ::   l_isend = .FALSE.  !: isend use indicator (T if cn_mpi_send='I')
   INTEGER         , PUBLIC ::   nn_buffer          !: size of the buffer in case of mpi_bsend

   REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE ::   tampon   ! buffer in case of bsend

   LOGICAL, PUBLIC ::   ln_nnogather                !: namelist control of northfold comms
   LOGICAL, PUBLIC ::   l_north_nogather = .FALSE.  !: internal control of northfold comms

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lib_mpp.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION mynode( ldtxt, ldname, kumnam_ref, kumnam_cfg, kumond, kstop, localComm )
      !!----------------------------------------------------------------------
      !!                  ***  routine mynode  ***
      !!
      !! ** Purpose :   Find processor unit
      !!----------------------------------------------------------------------
      CHARACTER(len=*),DIMENSION(:), INTENT(  out) ::   ldtxt        !
      CHARACTER(len=*)             , INTENT(in   ) ::   ldname       !
      INTEGER                      , INTENT(in   ) ::   kumnam_ref   ! logical unit for reference namelist
      INTEGER                      , INTENT(in   ) ::   kumnam_cfg   ! logical unit for configuration namelist
      INTEGER                      , INTENT(inout) ::   kumond       ! logical unit for namelist output
      INTEGER                      , INTENT(inout) ::   kstop        ! stop indicator
      INTEGER         , OPTIONAL   , INTENT(in   ) ::   localComm    !
      !
      INTEGER ::   mynode, ierr, code, ji, ii, ios
      LOGICAL ::   mpi_was_called
      !
      NAMELIST/nammpp/ cn_mpi_send, nn_buffer, jpni, jpnj, jpnij, ln_nnogather
      !!----------------------------------------------------------------------
      !
      ii = 1
      WRITE(ldtxt(ii),*)                                                                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'mynode : mpi initialisation'                                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~ '                                                        ;   ii = ii + 1
      !
      REWIND( kumnam_ref )              ! Namelist nammpp in reference namelist: mpi variables
      READ  ( kumnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nammpp in reference namelist', lwp )
      !
      REWIND( kumnam_cfg )              ! Namelist nammpp in configuration namelist: mpi variables
      READ  ( kumnam_cfg, nammpp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nammpp in configuration namelist', lwp )
      !
      !                              ! control print
      WRITE(ldtxt(ii),*) '   Namelist nammpp'                                             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      mpi send type          cn_mpi_send = ', cn_mpi_send       ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      size exported buffer   nn_buffer   = ', nn_buffer,' bytes';   ii = ii + 1
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN
         jpni  = Agrif_Parent(jpni )
         jpnj  = Agrif_Parent(jpnj )
         jpnij = Agrif_Parent(jpnij)
      ENDIF
#endif
      !
      IF( jpnij < 1 ) THEN         ! If jpnij is not specified in namelist then we calculate it
         jpnij = jpni * jpnj       ! this means there will be no land cutting out.
      ENDIF

      IF( jpni < 1 .OR. jpnj < 1  ) THEN
         WRITE(ldtxt(ii),*) '      jpni, jpnj and jpnij will be calculated automatically' ;   ii = ii + 1
      ELSE
         WRITE(ldtxt(ii),*) '      processor grid extent in i         jpni = ',jpni       ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '      processor grid extent in j         jpnj = ',jpnj       ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '      number of local domains           jpnij = ',jpnij      ;   ii = ii + 1
      ENDIF

      WRITE(ldtxt(ii),*) '      avoid use of mpi_allgather at the north fold  ln_nnogather = ', ln_nnogather  ; ii = ii + 1

      CALL mpi_initialized ( mpi_was_called, code )
      IF( code /= MPI_SUCCESS ) THEN
         DO ji = 1, SIZE(ldtxt)
            IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
         END DO
         WRITE(*, cform_err)
         WRITE(*, *) 'lib_mpp: Error in routine mpi_initialized'
         CALL mpi_abort( mpi_comm_world, code, ierr )
      ENDIF

      IF( mpi_was_called ) THEN
         !
         SELECT CASE ( cn_mpi_send )
         CASE ( 'S' )                ! Standard mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Standard blocking mpi send (send)'             ;   ii = ii + 1
         CASE ( 'B' )                ! Buffer mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Buffer blocking mpi send (bsend)'              ;   ii = ii + 1
            IF( Agrif_Root() )   CALL mpi_init_oce( ldtxt, ii, ierr )
         CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
            WRITE(ldtxt(ii),*) '           Immediate non-blocking send (isend)'           ;   ii = ii + 1
            l_isend = .TRUE.
         CASE DEFAULT
            WRITE(ldtxt(ii),cform_err)                                                    ;   ii = ii + 1
            WRITE(ldtxt(ii),*) '           bad value for cn_mpi_send = ', cn_mpi_send     ;   ii = ii + 1
            kstop = kstop + 1
         END SELECT
         !
      ELSEIF ( PRESENT(localComm) .AND. .NOT. mpi_was_called ) THEN
         WRITE(ldtxt(ii),*) ' lib_mpp: You cannot provide a local communicator '          ;   ii = ii + 1
         WRITE(ldtxt(ii),*) '          without calling MPI_Init before ! '                ;   ii = ii + 1
         kstop = kstop + 1
      ELSE
         SELECT CASE ( cn_mpi_send )
         CASE ( 'S' )                ! Standard mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Standard blocking mpi send (send)'             ;   ii = ii + 1
            CALL mpi_init( ierr )
         CASE ( 'B' )                ! Buffer mpi send (blocking)
            WRITE(ldtxt(ii),*) '           Buffer blocking mpi send (bsend)'              ;   ii = ii + 1
            IF( Agrif_Root() )   CALL mpi_init_oce( ldtxt, ii, ierr )
         CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
            WRITE(ldtxt(ii),*) '           Immediate non-blocking send (isend)'           ;   ii = ii + 1
            l_isend = .TRUE.
            CALL mpi_init( ierr )
         CASE DEFAULT
            WRITE(ldtxt(ii),cform_err)                                                    ;   ii = ii + 1
            WRITE(ldtxt(ii),*) '           bad value for cn_mpi_send = ', cn_mpi_send     ;   ii = ii + 1
            kstop = kstop + 1
         END SELECT
         !
      ENDIF

      IF( PRESENT(localComm) ) THEN
         IF( Agrif_Root() ) THEN
            mpi_comm_oce = localComm
         ENDIF
      ELSE
         CALL mpi_comm_dup( mpi_comm_world, mpi_comm_oce, code)
         IF( code /= MPI_SUCCESS ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in routine mpi_comm_dup'
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF
      ENDIF

#if defined key_agrif
      IF( Agrif_Root() ) THEN
         CALL Agrif_MPI_Init(mpi_comm_oce)
      ELSE
         CALL Agrif_MPI_set_grid_comm(mpi_comm_oce)
      ENDIF
#endif

      CALL mpi_comm_rank( mpi_comm_oce, mpprank, ierr )
      CALL mpi_comm_size( mpi_comm_oce, mppsize, ierr )
      mynode = mpprank

      IF( mynode == 0 ) THEN
         CALL ctl_opn( kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. , 1 )
         WRITE(kumond, nammpp)      
      ENDIF
      !
      CALL MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
      !
   END FUNCTION mynode

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_lnk_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_lnk_... routines
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
#     define ROUTINE_LNK           mpp_lnk_2d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_2d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_LNK           mpp_lnk_3d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_3d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_LNK           mpp_lnk_4d
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     define MULTI
#     define ROUTINE_LNK           mpp_lnk_4d_ptr
#     include "mpp_lnk_generic.h90"
#     undef ROUTINE_LNK
#     undef MULTI
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_nfd_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_nfd_... routines
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
#     define ROUTINE_NFD           mpp_nfd_2d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_2d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_NFD           mpp_nfd_3d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_3d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_NFD           mpp_nfd_4d
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     define MULTI
#     define ROUTINE_NFD           mpp_nfd_4d_ptr
#     include "mpp_nfd_generic.h90"
#     undef ROUTINE_NFD
#     undef MULTI
#  undef DIM_4d


   !!----------------------------------------------------------------------
   !!                   ***  routine mpp_lnk_bdy_(2,3,4)d  ***
   !!
   !!   * Argument : dummy argument use in mpp_lnk_... routines
   !!                ptab   :   array or pointer of arrays on which the boundary condition is applied
   !!                cd_nat :   nature of array grid-points
   !!                psgn   :   sign used across the north fold boundary
   !!                kb_bdy :   BDY boundary set
   !!                kfld   :   optional, number of pt3d arrays
   !!----------------------------------------------------------------------
   !
   !                       !==  2D array and array of 2D pointer  ==!
   !
#  define DIM_2d
#     define ROUTINE_BDY           mpp_lnk_bdy_2d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_2d
   !
   !                       !==  3D array and array of 3D pointer  ==!
   !
#  define DIM_3d
#     define ROUTINE_BDY           mpp_lnk_bdy_3d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_3d
   !
   !                       !==  4D array and array of 4D pointer  ==!
   !
#  define DIM_4d
#     define ROUTINE_BDY           mpp_lnk_bdy_4d
#     include "mpp_bdy_generic.h90"
#     undef ROUTINE_BDY
#  undef DIM_4d

   !!----------------------------------------------------------------------
   !!
   !!   load_array  &   mpp_lnk_2d_9    à generaliser a 3D et 4D
   
   
   !!    mpp_lnk_sum_2d et 3D   ====>>>>>>   à virer du code !!!!
   
   
   !!----------------------------------------------------------------------



   SUBROUTINE mppsend( ktyp, pmess, kbytes, kdest, md_req )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsend  ***
      !!
      !! ** Purpose :   Send messag passing array
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   pmess(*)   ! array of real
      INTEGER , INTENT(in   ) ::   kbytes     ! size of the array pmess
      INTEGER , INTENT(in   ) ::   kdest      ! receive process number
      INTEGER , INTENT(in   ) ::   ktyp       ! tag of the message
      INTEGER , INTENT(in   ) ::   md_req     ! argument for isend
      !!
      INTEGER ::   iflag
      !!----------------------------------------------------------------------
      !
      SELECT CASE ( cn_mpi_send )
      CASE ( 'S' )                ! Standard mpi send (blocking)
         CALL mpi_send ( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce        , iflag )
      CASE ( 'B' )                ! Buffer mpi send (blocking)
         CALL mpi_bsend( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce        , iflag )
      CASE ( 'I' )                ! Immediate mpi send (non-blocking send)
         ! be carefull, one more argument here : the mpi request identifier..
         CALL mpi_isend( pmess, kbytes, mpi_double_precision, kdest , ktyp, mpi_comm_oce, md_req, iflag )
      END SELECT
      !
   END SUBROUTINE mppsend


   SUBROUTINE mpprecv( ktyp, pmess, kbytes, ksource )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpprecv  ***
      !!
      !! ** Purpose :   Receive messag passing array
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout) ::   pmess(*)   ! array of real
      INTEGER , INTENT(in   ) ::   kbytes     ! suze of the array pmess
      INTEGER , INTENT(in   ) ::   ktyp       ! Tag of the recevied message
      INTEGER, OPTIONAL, INTENT(in) :: ksource    ! source process number
      !!
      INTEGER :: istatus(mpi_status_size)
      INTEGER :: iflag
      INTEGER :: use_source
      !!----------------------------------------------------------------------
      !
      ! If a specific process number has been passed to the receive call,
      ! use that one. Default is to use mpi_any_source
      use_source = mpi_any_source
      IF( PRESENT(ksource) )   use_source = ksource
      !
      CALL mpi_recv( pmess, kbytes, mpi_double_precision, use_source, ktyp, mpi_comm_oce, istatus, iflag )
      !
   END SUBROUTINE mpprecv


   SUBROUTINE mppgather( ptab, kp, pio )
      !!----------------------------------------------------------------------
      !!                   ***  routine mppgather  ***
      !!
      !! ** Purpose :   Transfert between a local subdomain array and a work
      !!     array which is distributed following the vertical level.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)      , INTENT(in   ) ::   ptab   ! subdomain input array
      INTEGER                           , INTENT(in   ) ::   kp     ! record length
      REAL(wp), DIMENSION(jpi,jpj,jpnij), INTENT(  out) ::   pio    ! subdomain input array
      !!
      INTEGER :: itaille, ierror   ! temporary integer
      !!---------------------------------------------------------------------
      !
      itaille = jpi * jpj
      CALL mpi_gather( ptab, itaille, mpi_double_precision, pio, itaille     ,   &
         &                            mpi_double_precision, kp , mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppgather


   SUBROUTINE mppscatter( pio, kp, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppscatter  ***
      !!
      !! ** Purpose :   Transfert between awork array which is distributed
      !!      following the vertical level and the local subdomain array.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpnij)  ::   pio    ! output array
      INTEGER                             ::   kp     ! Tag (not used with MPI
      REAL(wp), DIMENSION(jpi,jpj)        ::   ptab   ! subdomain array input
      !!
      INTEGER :: itaille, ierror   ! temporary integer
      !!---------------------------------------------------------------------
      !
      itaille = jpi * jpj
      !
      CALL mpi_scatter( pio, itaille, mpi_double_precision, ptab, itaille     ,   &
         &                            mpi_double_precision, kp  , mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppscatter

   !!----------------------------------------------------------------------
   !!    ***  mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
   SUBROUTINE mppmax_a_int( ktab, kdim, kcom )
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                  ::   kdim   ! size of array
      INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab   ! input array
      INTEGER , INTENT(in   ), OPTIONAL        ::   kcom   !
      INTEGER :: ierror, ilocalcomm   ! temporary integer
      INTEGER, DIMENSION(kdim) ::   iwork
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer, mpi_max, ilocalcomm, ierror )
      ktab(:) = iwork(:)
   END SUBROUTINE mppmax_a_int
   !!
   SUBROUTINE mppmax_int( ktab, kcom )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(inout)           ::   ktab   ! ???
      INTEGER, INTENT(in   ), OPTIONAL ::   kcom   ! ???
      INTEGER ::   ierror, iwork, ilocalcomm   ! temporary integer
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_max, ilocalcomm, ierror )
      ktab = iwork
   END SUBROUTINE mppmax_int
   !!
   SUBROUTINE mppmax_a_real( ptab, kdim, kcom )
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(kdim), INTENT(inout) ::   ptab
      INTEGER                  , INTENT(in   ) ::   kdim
      INTEGER , OPTIONAL       , INTENT(in   ) ::   kcom
      INTEGER :: ierror, ilocalcomm
      REAL(wp), DIMENSION(kdim) ::  zwork
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ptab, zwork, kdim, mpi_double_precision, mpi_max, ilocalcomm, ierror )
      ptab(:) = zwork(:)
   END SUBROUTINE mppmax_a_real
   !!
   SUBROUTINE mppmax_real( ptab, kcom )
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout)           ::   ptab   ! ???
      INTEGER , INTENT(in   ), OPTIONAL ::   kcom   ! ???
      INTEGER  ::   ierror, ilocalcomm
      REAL(wp) ::   zwork
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom!
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_max, ilocalcomm, ierror )
      ptab = zwork
   END SUBROUTINE mppmax_real


   !!----------------------------------------------------------------------
   !!    ***  mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real  ***
   !!   
   !!----------------------------------------------------------------------
   !!
   SUBROUTINE mppmin_a_int( ktab, kdim, kcom )
      !!----------------------------------------------------------------------
      INTEGER , INTENT( in  )                  ::   kdim   ! size of array
      INTEGER , INTENT(inout), DIMENSION(kdim) ::   ktab   ! input array
      INTEGER , INTENT( in  ), OPTIONAL        ::   kcom   ! input array
      !!
      INTEGER ::   ierror, ilocalcomm   ! temporary integer
      INTEGER, DIMENSION(kdim) ::   iwork
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer, mpi_min, ilocalcomm, ierror )
      ktab(:) = iwork(:)
   END SUBROUTINE mppmin_a_int
   !!
   SUBROUTINE mppmin_int( ktab, kcom )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(inout) ::   ktab      ! ???
      INTEGER , INTENT( in  ), OPTIONAL        ::   kcom        ! input array
      !!
      INTEGER ::  ierror, iwork, ilocalcomm
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_min, ilocalcomm, ierror )
      ktab = iwork
   END SUBROUTINE mppmin_int
   !!
   SUBROUTINE mppmin_a_real( ptab, kdim, kcom )
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                  ::   kdim
      REAL(wp), INTENT(inout), DIMENSION(kdim) ::   ptab
      INTEGER , INTENT(in   ), OPTIONAL        ::   kcom
      INTEGER :: ierror, ilocalcomm
      REAL(wp), DIMENSION(kdim) ::   zwork
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ptab, zwork, kdim, mpi_double_precision, mpi_min, ilocalcomm, ierror )
      ptab(:) = zwork(:)
   END SUBROUTINE mppmin_a_real
   !!
   SUBROUTINE mppmin_real( ptab, kcom )
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(inout)           ::   ptab        !
      INTEGER , INTENT(in   ), OPTIONAL :: kcom
      INTEGER  ::   ierror, ilocalcomm
      REAL(wp) ::   zwork
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_min, ilocalcomm, ierror )
      ptab = zwork
   END SUBROUTINE mppmin_real


   !!----------------------------------------------------------------------
   !!    ***  mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real  ***
   !!   
   !!   Global sum of 1D array or a variable (integer, real or complex)
   !!----------------------------------------------------------------------
   !!
   SUBROUTINE mppsum_a_int( ktab, kdim )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   )                   ::   kdim   ! ???
      INTEGER, INTENT(inout), DIMENSION (kdim) ::   ktab   ! ???
      INTEGER :: ierror
      INTEGER, DIMENSION (kdim) ::  iwork
      !!----------------------------------------------------------------------
      CALL mpi_allreduce( ktab, iwork, kdim, mpi_integer, mpi_sum, mpi_comm_oce, ierror )
      ktab(:) = iwork(:)
   END SUBROUTINE mppsum_a_int
   !!
   SUBROUTINE mppsum_int( ktab )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(inout) ::   ktab
      INTEGER :: ierror, iwork
      !!----------------------------------------------------------------------
      CALL mpi_allreduce( ktab, iwork, 1, mpi_integer, mpi_sum, mpi_comm_oce, ierror )
      ktab = iwork
   END SUBROUTINE mppsum_int
   !!
   SUBROUTINE mppsum_a_real( ptab, kdim, kcom )
      !!-----------------------------------------------------------------------
      INTEGER                  , INTENT(in   ) ::   kdim   ! size of ptab
      REAL(wp), DIMENSION(kdim), INTENT(inout) ::   ptab   ! input array
      INTEGER , OPTIONAL       , INTENT(in   ) ::   kcom   ! specific communicator
      INTEGER  ::   ierror, ilocalcomm    ! local integer
      REAL(wp) ::   zwork(kdim)           ! local workspace
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ptab, zwork, kdim, mpi_double_precision, mpi_sum, ilocalcomm, ierror )
      ptab(:) = zwork(:)
   END SUBROUTINE mppsum_a_real
   !!
   SUBROUTINE mppsum_real( ptab, kcom )
      !!-----------------------------------------------------------------------
      REAL(wp)          , INTENT(inout)           ::   ptab   ! input scalar
      INTEGER , OPTIONAL, INTENT(in   ) ::   kcom
      INTEGER  ::   ierror, ilocalcomm
      REAL(wp) ::   zwork
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL mpi_allreduce( ptab, zwork, 1, mpi_double_precision, mpi_sum, ilocalcomm, ierror )
      ptab = zwork
   END SUBROUTINE mppsum_real
   !!
   SUBROUTINE mppsum_realdd( ytab, kcom )
      !!-----------------------------------------------------------------------
      COMPLEX(wp)          , INTENT(inout) ::   ytab    ! input scalar
      INTEGER    , OPTIONAL, INTENT(in   ) ::   kcom
      INTEGER     ::   ierror, ilocalcomm
      COMPLEX(wp) ::   zwork
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL MPI_ALLREDUCE( ytab, zwork, 1, MPI_DOUBLE_COMPLEX, MPI_SUMDD, ilocalcomm, ierror )
      ytab = zwork
   END SUBROUTINE mppsum_realdd
   !!
   SUBROUTINE mppsum_a_realdd( ytab, kdim, kcom )
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   kdim   ! size of ytab
      COMPLEX(wp), DIMENSION(kdim), INTENT(inout) ::   ytab   ! input array
      INTEGER    , OPTIONAL       , INTENT(in   ) ::   kcom
      INTEGER:: ierror, ilocalcomm    ! local integer
      COMPLEX(wp), DIMENSION(kdim) :: zwork     ! temporary workspace
      !!-----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      CALL MPI_ALLREDUCE( ytab, zwork, kdim, MPI_DOUBLE_COMPLEX, MPI_SUMDD, ilocalcomm, ierror )
      ytab(:) = zwork(:)
   END SUBROUTINE mppsum_a_realdd
   

   SUBROUTINE mppmax_real_multiple( pt1d, kdim, kcom  )
      !!----------------------------------------------------------------------
      !!                  ***  routine mppmax_real  ***
      !!
      !! ** Purpose :   Maximum across processor of each element of a 1D arrays
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(kdim), INTENT(inout) ::   pt1d   ! 1D arrays
      INTEGER                  , INTENT(in   ) ::   kdim
      INTEGER , OPTIONAL       , INTENT(in   ) ::   kcom   ! local communicator
      !!
      INTEGER  ::   ierror, ilocalcomm
      REAL(wp), DIMENSION(kdim) ::  zwork
      !!----------------------------------------------------------------------
      ilocalcomm = mpi_comm_oce
      IF( PRESENT(kcom) )   ilocalcomm = kcom
      !
      CALL mpi_allreduce( pt1d, zwork, kdim, mpi_double_precision, mpi_max, ilocalcomm, ierror )
      pt1d(:) = zwork(:)
      !
   END SUBROUTINE mppmax_real_multiple


   SUBROUTINE mpp_minloc2d( ptab, pmask, pmin, ki,kj )
      !!------------------------------------------------------------------------
      !!             ***  routine mpp_minloc  ***
      !!
      !! ** Purpose :   Compute the global minimum of an array ptab
      !!              and also give its global position
      !!
      !! ** Method  :   Use MPI_ALLREDUCE with MPI_MINLOC
      !!
      !!--------------------------------------------------------------------------
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   ptab    ! Local 2D array
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pmask   ! Local mask
      REAL(wp)                     , INTENT(  out) ::   pmin    ! Global minimum of ptab
      INTEGER                      , INTENT(  out) ::   ki, kj  ! index of minimum in global frame
      !
      INTEGER :: ierror
      INTEGER , DIMENSION(2)   ::   ilocs
      REAL(wp) ::   zmin   ! local minimum
      REAL(wp), DIMENSION(2,1) ::   zain, zaout
      !!-----------------------------------------------------------------------
      !
      zmin  = MINVAL( ptab(:,:) , mask= pmask == 1._wp )
      ilocs = MINLOC( ptab(:,:) , mask= pmask == 1._wp )
      !
      ki = ilocs(1) + nimpp - 1
      kj = ilocs(2) + njmpp - 1
      !
      zain(1,:)=zmin
      zain(2,:)=ki+10000.*kj
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OCE,ierror)
      !
      pmin = zaout(1,1)
      kj = INT(zaout(2,1)/10000.)
      ki = INT(zaout(2,1) - 10000.*kj )
      !
   END SUBROUTINE mpp_minloc2d


   SUBROUTINE mpp_minloc3d( ptab, pmask, pmin, ki, kj ,kk)
      !!------------------------------------------------------------------------
      !!             ***  routine mpp_minloc  ***
      !!
      !! ** Purpose :   Compute the global minimum of an array ptab
      !!              and also give its global position
      !!
      !! ** Method  :   Use MPI_ALLREDUCE with MPI_MINLOC
      !!
      !!--------------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   ptab         ! Local 2D array
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pmask        ! Local mask
      REAL(wp)                  , INTENT(  out) ::   pmin         ! Global minimum of ptab
      INTEGER                   , INTENT(  out) ::   ki, kj, kk   ! index of minimum in global frame
      !
      INTEGER  ::   ierror
      REAL(wp) ::   zmin     ! local minimum
      INTEGER , DIMENSION(3)   ::   ilocs
      REAL(wp), DIMENSION(2,1) ::   zain, zaout
      !!-----------------------------------------------------------------------
      !
      zmin  = MINVAL( ptab(:,:,:) , mask= pmask == 1._wp )
      ilocs = MINLOC( ptab(:,:,:) , mask= pmask == 1._wp )
      !
      ki = ilocs(1) + nimpp - 1
      kj = ilocs(2) + njmpp - 1
      kk = ilocs(3)
      !
      zain(1,:) = zmin
      zain(2,:) = ki + 10000.*kj + 100000000.*kk
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_OCE,ierror)
      !
      pmin = zaout(1,1)
      kk   = INT( zaout(2,1) / 100000000. )
      kj   = INT( zaout(2,1) - kk * 100000000. ) / 10000
      ki   = INT( zaout(2,1) - kk * 100000000. -kj * 10000. )
      !
   END SUBROUTINE mpp_minloc3d


   SUBROUTINE mpp_maxloc2d( ptab, pmask, pmax, ki, kj )
      !!------------------------------------------------------------------------
      !!             ***  routine mpp_maxloc  ***
      !!
      !! ** Purpose :   Compute the global maximum of an array ptab
      !!              and also give its global position
      !!
      !! ** Method  :   Use MPI_ALLREDUCE with MPI_MINLOC
      !!
      !!--------------------------------------------------------------------------
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   ptab     ! Local 2D array
      REAL(wp), DIMENSION (jpi,jpj), INTENT(in   ) ::   pmask    ! Local mask
      REAL(wp)                     , INTENT(  out) ::   pmax     ! Global maximum of ptab
      INTEGER                      , INTENT(  out) ::   ki, kj   ! index of maximum in global frame
      !!
      INTEGER  :: ierror
      INTEGER, DIMENSION (2)   ::   ilocs
      REAL(wp) :: zmax   ! local maximum
      REAL(wp), DIMENSION(2,1) ::   zain, zaout
      !!-----------------------------------------------------------------------
      !
      zmax  = MAXVAL( ptab(:,:) , mask= pmask == 1._wp )
      ilocs = MAXLOC( ptab(:,:) , mask= pmask == 1._wp )
      !
      ki = ilocs(1) + nimpp - 1
      kj = ilocs(2) + njmpp - 1
      !
      zain(1,:) = zmax
      zain(2,:) = ki + 10000. * kj
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OCE,ierror)
      !
      pmax = zaout(1,1)
      kj   = INT( zaout(2,1) / 10000.     )
      ki   = INT( zaout(2,1) - 10000.* kj )
      !
   END SUBROUTINE mpp_maxloc2d


   SUBROUTINE mpp_maxloc3d( ptab, pmask, pmax, ki, kj, kk )
      !!------------------------------------------------------------------------
      !!             ***  routine mpp_maxloc  ***
      !!
      !! ** Purpose :  Compute the global maximum of an array ptab
      !!              and also give its global position
      !!
      !! ** Method : Use MPI_ALLREDUCE with MPI_MINLOC
      !!
      !!--------------------------------------------------------------------------
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   ptab         ! Local 2D array
      REAL(wp), DIMENSION (:,:,:), INTENT(in   ) ::   pmask        ! Local mask
      REAL(wp)                   , INTENT(  out) ::   pmax         ! Global maximum of ptab
      INTEGER                    , INTENT(  out) ::   ki, kj, kk   ! index of maximum in global frame
      !
      INTEGER  ::   ierror   ! local integer
      REAL(wp) ::   zmax     ! local maximum
      REAL(wp), DIMENSION(2,1) ::   zain, zaout
      INTEGER , DIMENSION(3)   ::   ilocs
      !!-----------------------------------------------------------------------
      !
      zmax  = MAXVAL( ptab(:,:,:) , mask= pmask == 1._wp )
      ilocs = MAXLOC( ptab(:,:,:) , mask= pmask == 1._wp )
      !
      ki = ilocs(1) + nimpp - 1
      kj = ilocs(2) + njmpp - 1
      kk = ilocs(3)
      !
      zain(1,:) = zmax
      zain(2,:) = ki + 10000.*kj + 100000000.*kk
      !
      CALL MPI_ALLREDUCE( zain,zaout, 1, MPI_2DOUBLE_PRECISION,MPI_MAXLOC,MPI_COMM_OCE,ierror)
      !
      pmax = zaout(1,1)
      kk   = INT( zaout(2,1) / 100000000. )
      kj   = INT( zaout(2,1) - kk * 100000000. ) / 10000
      ki   = INT( zaout(2,1) - kk * 100000000. -kj * 10000. )
      !
   END SUBROUTINE mpp_maxloc3d


   SUBROUTINE mppsync()
      !!----------------------------------------------------------------------
      !!                  ***  routine mppsync  ***
      !!
      !! ** Purpose :   Massively parallel processors, synchroneous
      !!
      !!-----------------------------------------------------------------------
      INTEGER :: ierror
      !!-----------------------------------------------------------------------
      !
      CALL mpi_barrier( mpi_comm_oce, ierror )
      !
   END SUBROUTINE mppsync


   SUBROUTINE mppstop
      !!----------------------------------------------------------------------
      !!                  ***  routine mppstop  ***
      !!
      !! ** purpose :   Stop massively parallel processors method
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   info
      !!----------------------------------------------------------------------
      !
      CALL mppsync
      CALL mpi_finalize( info )
      !
   END SUBROUTINE mppstop


   SUBROUTINE mpp_comm_free( kcom )
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kcom
      !!
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      CALL MPI_COMM_FREE(kcom, ierr)
      !
   END SUBROUTINE mpp_comm_free


   SUBROUTINE mpp_ini_ice( pindic, kumout )
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_ice  ***
      !!
      !! ** Purpose :   Initialize special communicator for ice areas
      !!      condition together with global variables needed in the ddmpp folding
      !!
      !! ** Method  : - Look for ice processors in ice routines
      !!              - Put their number in nrank_ice
      !!              - Create groups for the world processors and the ice processors
      !!              - Create a communicator for ice processors
      !!
      !! ** output
      !!      njmppmax = njmpp for northern procs
      !!      ndim_rank_ice = number of processors with ice
      !!      nrank_ice (ndim_rank_ice) = ice processors
      !!      ngrp_iworld = group ID for the world processors
      !!      ngrp_ice = group ID for the ice processors
      !!      ncomm_ice = communicator for the ice procs.
      !!      n_ice_root = number (in the world) of proc 0 in the ice comm.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   pindic
      INTEGER, INTENT(in) ::   kumout   ! ocean.output logical unit
      !!
      INTEGER :: jjproc
      INTEGER :: ii, ierr
      INTEGER, ALLOCATABLE, DIMENSION(:) ::   kice
      INTEGER, ALLOCATABLE, DIMENSION(:) ::   zwork
      !!----------------------------------------------------------------------
      !
      ALLOCATE( kice(jpnij), zwork(jpnij), STAT=ierr )
      IF( ierr /= 0 ) THEN
         WRITE(kumout, cform_err)
         WRITE(kumout,*) 'mpp_ini_ice : failed to allocate 2, 1D arrays (jpnij in length)'
         CALL mppstop
      ENDIF

      ! Look for how many procs with sea-ice
      !
      kice = 0
      DO jjproc = 1, jpnij
         IF( jjproc == narea .AND. pindic .GT. 0 )   kice(jjproc) = 1
      END DO
      !
      zwork = 0
      CALL MPI_ALLREDUCE( kice, zwork, jpnij, mpi_integer, mpi_sum, mpi_comm_oce, ierr )
      ndim_rank_ice = SUM( zwork )

      ! Allocate the right size to nrank_north
      IF( ALLOCATED ( nrank_ice ) )   DEALLOCATE( nrank_ice )
      ALLOCATE( nrank_ice(ndim_rank_ice) )
      !
      ii = 0
      nrank_ice = 0
      DO jjproc = 1, jpnij
         IF( zwork(jjproc) == 1) THEN
            ii = ii + 1
            nrank_ice(ii) = jjproc -1
         ENDIF
      END DO

      ! Create the world group
      CALL MPI_COMM_GROUP( mpi_comm_oce, ngrp_iworld, ierr )

      ! Create the ice group from the world group
      CALL MPI_GROUP_INCL( ngrp_iworld, ndim_rank_ice, nrank_ice, ngrp_ice, ierr )

      ! Create the ice communicator , ie the pool of procs with sea-ice
      CALL MPI_COMM_CREATE( mpi_comm_oce, ngrp_ice, ncomm_ice, ierr )

      ! Find proc number in the world of proc 0 in the north
      ! The following line seems to be useless, we just comment & keep it as reminder
      ! CALL MPI_GROUP_TRANSLATE_RANKS(ngrp_ice,1,0,ngrp_iworld,n_ice_root,ierr)
      !
      CALL MPI_GROUP_FREE(ngrp_ice, ierr)
      CALL MPI_GROUP_FREE(ngrp_iworld, ierr)

      DEALLOCATE(kice, zwork)
      !
   END SUBROUTINE mpp_ini_ice


   SUBROUTINE mpp_ini_znl( kumout )
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_znl  ***
      !!
      !! ** Purpose :   Initialize special communicator for computing zonal sum
      !!
      !! ** Method  : - Look for processors in the same row
      !!              - Put their number in nrank_znl
      !!              - Create group for the znl processors
      !!              - Create a communicator for znl processors
      !!              - Determine if processor should write znl files
      !!
      !! ** output
      !!      ndim_rank_znl = number of processors on the same row
      !!      ngrp_znl = group ID for the znl processors
      !!      ncomm_znl = communicator for the ice procs.
      !!      n_znl_root = number (in the world) of proc 0 in the ice comm.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kumout   ! ocean.output logical units
      !
      INTEGER :: jproc      ! dummy loop integer
      INTEGER :: ierr, ii   ! local integer
      INTEGER, ALLOCATABLE, DIMENSION(:) ::   kwork
      !!----------------------------------------------------------------------
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_world     : ', ngrp_world
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - mpi_comm_world : ', mpi_comm_world
      !-$$     WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - mpi_comm_oce   : ', mpi_comm_oce
      !
      ALLOCATE( kwork(jpnij), STAT=ierr )
      IF( ierr /= 0 ) THEN
         WRITE(kumout, cform_err)
         WRITE(kumout,*) 'mpp_ini_znl : failed to allocate 1D array of length jpnij'
         CALL mppstop
      ENDIF

      IF( jpnj == 1 ) THEN
         ngrp_znl  = ngrp_world
         ncomm_znl = mpi_comm_oce
      ELSE
         !
         CALL MPI_ALLGATHER ( njmpp, 1, mpi_integer, kwork, 1, mpi_integer, mpi_comm_oce, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - kwork pour njmpp : ', kwork
         !-$$        CALL flush(numout)
         !
         ! Count number of processors on the same row
         ndim_rank_znl = 0
         DO jproc=1,jpnij
            IF ( kwork(jproc) == njmpp ) THEN
               ndim_rank_znl = ndim_rank_znl + 1
            ENDIF
         END DO
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ndim_rank_znl : ', ndim_rank_znl
         !-$$        CALL flush(numout)
         ! Allocate the right size to nrank_znl
         IF (ALLOCATED (nrank_znl)) DEALLOCATE(nrank_znl)
         ALLOCATE(nrank_znl(ndim_rank_znl))
         ii = 0
         nrank_znl (:) = 0
         DO jproc=1,jpnij
            IF ( kwork(jproc) == njmpp) THEN
               ii = ii + 1
               nrank_znl(ii) = jproc -1
            ENDIF
         END DO
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - nrank_znl : ', nrank_znl
         !-$$        CALL flush(numout)

         ! Create the opa group
         CALL MPI_COMM_GROUP(mpi_comm_oce,ngrp_opa,ierr)
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_opa : ', ngrp_opa
         !-$$        CALL flush(numout)

         ! Create the znl group from the opa group
         CALL MPI_GROUP_INCL  ( ngrp_opa, ndim_rank_znl, nrank_znl, ngrp_znl, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ngrp_znl ', ngrp_znl
         !-$$        CALL flush(numout)

         ! Create the znl communicator from the opa communicator, ie the pool of procs in the same row
         CALL MPI_COMM_CREATE ( mpi_comm_oce, ngrp_znl, ncomm_znl, ierr )
         !-$$        WRITE (numout,*) 'mpp_ini_znl ', nproc, ' - ncomm_znl ', ncomm_znl
         !-$$        CALL flush(numout)
         !
      END IF

      ! Determines if processor if the first (starting from i=1) on the row
      IF ( jpni == 1 ) THEN
         l_znl_root = .TRUE.
      ELSE
         l_znl_root = .FALSE.
         kwork (1) = nimpp
         CALL mpp_min ( kwork(1), kcom = ncomm_znl)
         IF ( nimpp == kwork(1)) l_znl_root = .TRUE.
      END IF

      DEALLOCATE(kwork)

   END SUBROUTINE mpp_ini_znl


   SUBROUTINE mpp_ini_north
      !!----------------------------------------------------------------------
      !!               ***  routine mpp_ini_north  ***
      !!
      !! ** Purpose :   Initialize special communicator for north folding
      !!      condition together with global variables needed in the mpp folding
      !!
      !! ** Method  : - Look for northern processors
      !!              - Put their number in nrank_north
      !!              - Create groups for the world processors and the north processors
      !!              - Create a communicator for northern processors
      !!
      !! ** output
      !!      njmppmax = njmpp for northern procs
      !!      ndim_rank_north = number of processors in the northern line
      !!      nrank_north (ndim_rank_north) = number  of the northern procs.
      !!      ngrp_world = group ID for the world processors
      !!      ngrp_north = group ID for the northern processors
      !!      ncomm_north = communicator for the northern procs.
      !!      north_root = number (in the world) of proc 0 in the northern comm.
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ierr
      INTEGER ::   jjproc
      INTEGER ::   ii, ji
      !!----------------------------------------------------------------------
      !
      njmppmax = MAXVAL( njmppt )
      !
      ! Look for how many procs on the northern boundary
      ndim_rank_north = 0
      DO jjproc = 1, jpnij
         IF( njmppt(jjproc) == njmppmax )   ndim_rank_north = ndim_rank_north + 1
      END DO
      !
      ! Allocate the right size to nrank_north
      IF (ALLOCATED (nrank_north)) DEALLOCATE(nrank_north)
      ALLOCATE( nrank_north(ndim_rank_north) )

      ! Fill the nrank_north array with proc. number of northern procs.
      ! Note : the rank start at 0 in MPI
      ii = 0
      DO ji = 1, jpnij
         IF ( njmppt(ji) == njmppmax   ) THEN
            ii=ii+1
            nrank_north(ii)=ji-1
         END IF
      END DO
      !
      ! create the world group
      CALL MPI_COMM_GROUP( mpi_comm_oce, ngrp_world, ierr )
      !
      ! Create the North group from the world group
      CALL MPI_GROUP_INCL( ngrp_world, ndim_rank_north, nrank_north, ngrp_north, ierr )
      !
      ! Create the North communicator , ie the pool of procs in the north group
      CALL MPI_COMM_CREATE( mpi_comm_oce, ngrp_north, ncomm_north, ierr )
      !
   END SUBROUTINE mpp_ini_north


   SUBROUTINE mpi_init_oce( ldtxt, ksft, code )
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_init.opa  ***
      !!
      !! ** Purpose :: export and attach a MPI buffer for bsend
      !!
      !! ** Method  :: define buffer size in namelist, if 0 no buffer attachment
      !!            but classical mpi_init
      !!
      !! History :: 01/11 :: IDRIS initial version for IBM only
      !!            08/04 :: R. Benshila, generalisation
      !!---------------------------------------------------------------------
      CHARACTER(len=*),DIMENSION(:), INTENT(  out) ::   ldtxt
      INTEGER                      , INTENT(inout) ::   ksft
      INTEGER                      , INTENT(  out) ::   code
      INTEGER                                      ::   ierr, ji
      LOGICAL                                      ::   mpi_was_called
      !!---------------------------------------------------------------------
      !
      CALL mpi_initialized( mpi_was_called, code )      ! MPI initialization
      IF ( code /= MPI_SUCCESS ) THEN
         DO ji = 1, SIZE(ldtxt)
            IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
         END DO
         WRITE(*, cform_err)
         WRITE(*, *) ' lib_mpp: Error in routine mpi_initialized'
         CALL mpi_abort( mpi_comm_world, code, ierr )
      ENDIF
      !
      IF( .NOT. mpi_was_called ) THEN
         CALL mpi_init( code )
         CALL mpi_comm_dup( mpi_comm_world, mpi_comm_oce, code )
         IF ( code /= MPI_SUCCESS ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in routine mpi_comm_dup'
            CALL mpi_abort( mpi_comm_world, code, ierr )
         ENDIF
      ENDIF
      !
      IF( nn_buffer > 0 ) THEN
         WRITE(ldtxt(ksft),*) 'mpi_bsend, buffer allocation of  : ', nn_buffer   ;   ksft = ksft + 1
         ! Buffer allocation and attachment
         ALLOCATE( tampon(nn_buffer), stat = ierr )
         IF( ierr /= 0 ) THEN
            DO ji = 1, SIZE(ldtxt)
               IF( TRIM(ldtxt(ji)) /= '' )   WRITE(*,*) ldtxt(ji)      ! control print of mynode
            END DO
            WRITE(*, cform_err)
            WRITE(*, *) ' lib_mpp: Error in ALLOCATE', ierr
            CALL mpi_abort( mpi_comm_world, code, ierr )
         END IF
         CALL mpi_buffer_attach( tampon, nn_buffer, code )
      ENDIF
      !
   END SUBROUTINE mpi_init_oce


   SUBROUTINE DDPDD_MPI( ydda, yddb, ilen, itype )
      !!---------------------------------------------------------------------
      !!   Routine DDPDD_MPI: used by reduction operator MPI_SUMDD
      !!
      !!   Modification of original codes written by David H. Bailey
      !!   This subroutine computes yddb(i) = ydda(i)+yddb(i)
      !!---------------------------------------------------------------------
      INTEGER                     , INTENT(in)    ::   ilen, itype
      COMPLEX(wp), DIMENSION(ilen), INTENT(in)    ::   ydda
      COMPLEX(wp), DIMENSION(ilen), INTENT(inout) ::   yddb
      !
      REAL(wp) :: zerr, zt1, zt2    ! local work variables
      INTEGER  :: ji, ztmp           ! local scalar
      !!---------------------------------------------------------------------
      !
      ztmp = itype   ! avoid compilation warning
      !
      DO ji=1,ilen
      ! Compute ydda + yddb using Knuth's trick.
         zt1  = real(ydda(ji)) + real(yddb(ji))
         zerr = zt1 - real(ydda(ji))
         zt2  = ((real(yddb(ji)) - zerr) + (real(ydda(ji)) - (zt1 - zerr))) &
                + aimag(ydda(ji)) + aimag(yddb(ji))

         ! The result is zt1 + zt2, after normalization.
         yddb(ji) = cmplx ( zt1 + zt2, zt2 - ((zt1 + zt2) - zt1),wp )
      END DO
      !
   END SUBROUTINE DDPDD_MPI


   SUBROUTINE mpp_lbc_north_icb( pt2d, cd_type, psgn, kextj)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type  ! nature of pt3d grid-points
      !                                                     !   = T ,  U , V , F or W -points
      REAL(wp)                , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                    ! north fold, =  1. otherwise
      INTEGER                 , INTENT(in   ) ::   kextj    ! Extra halo width at north fold
      !
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ipj, ij, iproc
      !
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE  ::  ztab_e, znorthloc_e
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE  ::  znorthgloio_e
      !!----------------------------------------------------------------------
      !
      ipj=4
      ALLOCATE(        ztab_e(jpiglo, 1-kextj:ipj+kextj)       ,       &
     &            znorthloc_e(jpimax, 1-kextj:ipj+kextj)       ,       &
     &          znorthgloio_e(jpimax, 1-kextj:ipj+kextj,jpni)    )
      !
      ztab_e(:,:)      = 0._wp
      znorthloc_e(:,:) = 0._wp
      !
      ij = 1 - kextj
      ! put the last ipj+2*kextj lines of pt2d into znorthloc_e 
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         znorthloc_e(1:jpi,ij)=pt2d(1:jpi,jj)
         ij = ij + 1
      END DO
      !
      itaille = jpimax * ( ipj + 2*kextj )
      CALL MPI_ALLGATHER( znorthloc_e(1,1-kextj)    , itaille, MPI_DOUBLE_PRECISION,    &
         &                znorthgloio_e(1,1-kextj,1), itaille, MPI_DOUBLE_PRECISION,    &
         &                ncomm_north, ierr )
      !
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nrank_north(jr) + 1
         ildi = nldit (iproc)
         ilei = nleit (iproc)
         iilb = nimppt(iproc)
         DO jj = 1-kextj, ipj+kextj
            DO ji = ildi, ilei
               ztab_e(ji+iilb-1,jj) = znorthgloio_e(ji,jj,jr)
            END DO
         END DO
      END DO

      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd( ztab_e(:,1-kextj:ipj+kextj), cd_type, psgn, kextj )

      ij = 1 - kextj
      !! Scatter back to pt2d
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         DO ji= 1, jpi
            pt2d(ji,jj) = ztab_e(ji+nimpp-1,ij)
         END DO
         ij  = ij +1
      END DO
      !
      DEALLOCATE( ztab_e, znorthloc_e, znorthgloio_e )
      !
   END SUBROUTINE mpp_lbc_north_icb


   SUBROUTINE mpp_lnk_2d_icb( pt2d, cd_type, psgn, kexti, kextj )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    kexti  : number of columns for extra outer halo
      !!                    kextj  : number of rows for extra outer halo
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,1-kextj:jpj+kextj), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)                                        , INTENT(in   ) ::   cd_type  ! nature of ptab array grid-points
      REAL(wp)                                                , INTENT(in   ) ::   psgn     ! sign used across the north fold
      INTEGER                                                 , INTENT(in   ) ::   kexti    ! extra i-halo width
      INTEGER                                                 , INTENT(in   ) ::   kextj    ! extra j-halo width
      !
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! local integers
      INTEGER  ::   ipreci, iprecj             !   -       -
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!
      REAL(wp), DIMENSION(1-kexti:jpi+kexti,nn_hls+kextj,2) ::   r2dns, r2dsn
      REAL(wp), DIMENSION(1-kextj:jpj+kextj,nn_hls+kexti,2) ::   r2dwe, r2dew
      !!----------------------------------------------------------------------

      ipreci = nn_hls + kexti      ! take into account outer extra 2D overlap area
      iprecj = nn_hls + kextj


      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( l_Iperio ) THEN
         pt2d(1-kexti:     1   ,:) = pt2d(jpim1-kexti: jpim1 ,:)       ! east
         pt2d(  jpi  :jpi+kexti,:) = pt2d(     2     :2+kexti,:)       ! west
         !
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(  1-kexti   :nn_hls   ,:) = 0._wp    ! east except at F-point
                                      pt2d(jpi-nn_hls+1:jpi+kexti,:) = 0._wp    ! west
      ENDIF
      !                                      ! North-South boundaries
      IF( l_Jperio ) THEN                         !* cyclic (only with no mpp j-split)
         pt2d(:,1-kextj:     1   ) = pt2d(:,jpjm1-kextj:  jpjm1)       ! north
         pt2d(:,  jpj  :jpj+kextj) = pt2d(:,     2     :2+kextj)       ! south
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(:,  1-kextj   :nn_hls   ) = 0._wp    ! north except at F-point
                                      pt2d(:,jpj-nn_hls+1:jpj+kextj) = 0._wp    ! south
      ENDIF
      !

      ! north fold treatment
      ! -----------------------
      IF( npolj /= 0 ) THEN
         !
         SELECT CASE ( jpni )
                   CASE ( 1 )     ;   CALL lbc_nfd          ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
                   CASE DEFAULT   ;   CALL mpp_lbc_north_icb( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
         END SELECT
         !
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = jpi-nreci-kexti
         DO jl = 1, ipreci
            r2dew(:,jl,1) = pt2d(nn_hls+jl,:)
            r2dwe(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req1 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mppsend( 2, r2dwe(1-kextj,1,1), imigr, noea, ml_req2 )
         CALL mpprecv( 1, r2dew(1-kextj,1,2), imigr, noea )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, r2dew(1-kextj,1,1), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, r2dwe(1-kextj,1,2), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      !                           ! Write Dirichlet lateral conditions
      iihom = jpi - nn_hls
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
         END DO
      END SELECT


      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = jpj-nrecj-kextj
         DO jl = 1, iprecj
            r2dsn(:,jl,1) = pt2d(:,ijhom +jl)
            r2dns(:,jl,1) = pt2d(:,nn_hls+jl)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*kexti )
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req1 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mppsend( 4, r2dsn(1-kexti,1,1), imigr, nono, ml_req2 )
         CALL mpprecv( 3, r2dns(1-kexti,1,2), imigr, nono )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, r2dns(1-kexti,1,1), imigr, noso, ml_req1 )
         CALL mpprecv( 4, r2dsn(1-kexti,1,2), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = jpj - nn_hls
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      CASE ( 1 )
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
         END DO
      END SELECT
      !
   END SUBROUTINE mpp_lnk_2d_icb
   
#else
   !!----------------------------------------------------------------------
   !!   Default case:            Dummy module        share memory computing
   !!----------------------------------------------------------------------
   USE in_out_manager

   INTERFACE mpp_sum
      MODULE PROCEDURE mpp_sum_a2s, mpp_sum_as, mpp_sum_ai, mpp_sum_s, mpp_sum_i, mppsum_realdd, mppsum_a_realdd
   END INTERFACE
   INTERFACE mpp_max
      MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
   END INTERFACE
   INTERFACE mpp_min
      MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
   END INTERFACE
   INTERFACE mpp_minloc
      MODULE PROCEDURE mpp_minloc2d ,mpp_minloc3d
   END INTERFACE
   INTERFACE mpp_maxloc
      MODULE PROCEDURE mpp_maxloc2d ,mpp_maxloc3d
   END INTERFACE
   INTERFACE mpp_max_multiple
      MODULE PROCEDURE mppmax_real_multiple
   END INTERFACE

   LOGICAL, PUBLIC, PARAMETER ::   lk_mpp = .FALSE.      !: mpp flag
   LOGICAL, PUBLIC            ::   ln_nnogather          !: namelist control of northfold comms (needed here in case "key_mpp_mpi" is not used)
   INTEGER :: ncomm_ice
   INTEGER, PUBLIC            ::   mpi_comm_oce          ! opa local communicator
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION lib_mpp_alloc(kumout)          ! Dummy function
      INTEGER, INTENT(in) ::   kumout
      lib_mpp_alloc = 0
   END FUNCTION lib_mpp_alloc

   FUNCTION mynode( ldtxt, ldname, kumnam_ref, knumnam_cfg,  kumond , kstop, localComm ) RESULT (function_value)
      INTEGER, OPTIONAL            , INTENT(in   ) ::   localComm
      CHARACTER(len=*),DIMENSION(:) ::   ldtxt
      CHARACTER(len=*) ::   ldname
      INTEGER ::   kumnam_ref, knumnam_cfg , kumond , kstop
      IF( PRESENT( localComm ) ) mpi_comm_oce = localComm
      function_value = 0
      IF( .FALSE. )   ldtxt(:) = 'never done'
      CALL ctl_opn( kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE. , 1 )
   END FUNCTION mynode

   SUBROUTINE mppsync                       ! Dummy routine
   END SUBROUTINE mppsync

   SUBROUTINE mpp_sum_as( parr, kdim, kcom )      ! Dummy routine
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_as: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mpp_sum_as

   SUBROUTINE mpp_sum_a2s( parr, kdim, kcom )      ! Dummy routine
      REAL   , DIMENSION(:,:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_a2s: You should not have seen this print! error?', kdim, parr(1,1), kcom
   END SUBROUTINE mpp_sum_a2s

   SUBROUTINE mpp_sum_ai( karr, kdim, kcom )      ! Dummy routine
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_ai: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mpp_sum_ai

   SUBROUTINE mpp_sum_s( psca, kcom )            ! Dummy routine
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_s: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mpp_sum_s

   SUBROUTINE mpp_sum_i( kint, kcom )            ! Dummy routine
      integer               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mpp_sum_i: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mpp_sum_i

   SUBROUTINE mppsum_realdd( ytab, kcom )
      COMPLEX(wp), INTENT(inout)         :: ytab    ! input scalar
      INTEGER , INTENT( in  ), OPTIONAL :: kcom
      WRITE(*,*) 'mppsum_realdd: You should not have seen this print! error?', ytab
   END SUBROUTINE mppsum_realdd

   SUBROUTINE mppsum_a_realdd( ytab, kdim, kcom )
      INTEGER , INTENT( in )                        ::   kdim      ! size of ytab
      COMPLEX(wp), DIMENSION(kdim), INTENT( inout ) ::   ytab      ! input array
      INTEGER , INTENT( in  ), OPTIONAL :: kcom
      WRITE(*,*) 'mppsum_a_realdd: You should not have seen this print! error?', kdim, ytab(1), kcom
   END SUBROUTINE mppsum_a_realdd

   SUBROUTINE mppmax_a_real( parr, kdim, kcom )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mppmax_a_real

   SUBROUTINE mppmax_real( psca, kcom )
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_real: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mppmax_real

   SUBROUTINE mppmin_a_real( parr, kdim, kcom )
      REAL   , DIMENSION(:) :: parr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
   END SUBROUTINE mppmin_a_real

   SUBROUTINE mppmin_real( psca, kcom )
      REAL                  :: psca
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_real: You should not have seen this print! error?', psca, kcom
   END SUBROUTINE mppmin_real

   SUBROUTINE mppmax_a_int( karr, kdim ,kcom)
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mppmax_a_int

   SUBROUTINE mppmax_int( kint, kcom)
      INTEGER               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmax_int: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mppmax_int

   SUBROUTINE mppmin_a_int( karr, kdim, kcom )
      INTEGER, DIMENSION(:) :: karr
      INTEGER               :: kdim
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
   END SUBROUTINE mppmin_a_int

   SUBROUTINE mppmin_int( kint, kcom )
      INTEGER               :: kint
      INTEGER, OPTIONAL     :: kcom
      WRITE(*,*) 'mppmin_int: You should not have seen this print! error?', kint, kcom
   END SUBROUTINE mppmin_int

   SUBROUTINE mpp_minloc2d( ptab, pmask, pmin, ki, kj )
      REAL                   :: pmin
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mpp_minloc2d: You should not have seen this print! error?', pmin, ki, kj, ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_minloc2d

   SUBROUTINE mpp_minloc3d( ptab, pmask, pmin, ki, kj, kk )
      REAL                     :: pmin
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mpp_minloc3d: You should not have seen this print! error?', pmin, ki, kj, kk, ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_minloc3d

   SUBROUTINE mpp_maxloc2d( ptab, pmask, pmax, ki, kj )
      REAL                   :: pmax
      REAL , DIMENSION (:,:) :: ptab, pmask
      INTEGER :: ki, kj
      WRITE(*,*) 'mpp_maxloc2d: You should not have seen this print! error?', pmax, ki, kj, ptab(1,1), pmask(1,1)
   END SUBROUTINE mpp_maxloc2d

   SUBROUTINE mpp_maxloc3d( ptab, pmask, pmax, ki, kj, kk )
      REAL                     :: pmax
      REAL , DIMENSION (:,:,:) :: ptab, pmask
      INTEGER :: ki, kj, kk
      WRITE(*,*) 'mpp_maxloc3d: You should not have seen this print! error?', pmax, ki, kj, kk, ptab(1,1,1), pmask(1,1,1)
   END SUBROUTINE mpp_maxloc3d

   SUBROUTINE mppstop
      STOP      ! non MPP case, just stop the run
   END SUBROUTINE mppstop

   SUBROUTINE mpp_ini_ice( kcom, knum )
      INTEGER :: kcom, knum
      WRITE(*,*) 'mpp_ini_ice: You should not have seen this print! error?', kcom, knum
   END SUBROUTINE mpp_ini_ice

   SUBROUTINE mpp_ini_znl( knum )
      INTEGER :: knum
      WRITE(*,*) 'mpp_ini_znl: You should not have seen this print! error?', knum
   END SUBROUTINE mpp_ini_znl

   SUBROUTINE mpp_comm_free( kcom )
      INTEGER :: kcom
      WRITE(*,*) 'mpp_comm_free: You should not have seen this print! error?', kcom
   END SUBROUTINE mpp_comm_free
   
   SUBROUTINE mppmax_real_multiple( ptab, kdim , kcom  )
      REAL, DIMENSION(:) ::   ptab   ! 
      INTEGER            ::   kdim   ! 
      INTEGER, OPTIONAL  ::   kcom   ! 
      WRITE(*,*) 'mppmax_real_multiple: You should not have seen this print! error?', ptab(1), kdim
   END SUBROUTINE mppmax_real_multiple

#endif

   !!----------------------------------------------------------------------
   !!   All cases:         ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam   routines
   !!----------------------------------------------------------------------

   SUBROUTINE ctl_stop( cd1, cd2, cd3, cd4, cd5 ,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_opa  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the error number (nstop) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nstop = nstop + 1
      IF(lwp) THEN
         WRITE(numout,cform_err)
         IF( PRESENT(cd1 ) )   WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) )   WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) )   WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) )   WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) )   WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) )   WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) )   WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) )   WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) )   WRITE(numout,*) cd9
         IF( PRESENT(cd10) )   WRITE(numout,*) cd10
      ENDIF
                               CALL FLUSH(numout    )
      IF( numstp     /= -1 )   CALL FLUSH(numstp    )
      IF( numrun     /= -1 )   CALL FLUSH(numrun    )
      IF( numevo_ice /= -1 )   CALL FLUSH(numevo_ice)
      !
      IF( cd1 == 'STOP' ) THEN
         IF(lwp) WRITE(numout,*)  'huge E-R-R-O-R : immediate stop'
         CALL mppstop()
      ENDIF
      !
   END SUBROUTINE ctl_stop


   SUBROUTINE ctl_warn( cd1, cd2, cd3, cd4, cd5,   &
      &                 cd6, cd7, cd8, cd9, cd10 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE  stop_warn  ***
      !!
      !! ** Purpose :   print in ocean.outpput file a error message and
      !!                increment the warning number (nwarn) by one.
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd1, cd2, cd3, cd4, cd5
      CHARACTER(len=*), INTENT(in), OPTIONAL ::  cd6, cd7, cd8, cd9, cd10
      !!----------------------------------------------------------------------
      !
      nwarn = nwarn + 1
      IF(lwp) THEN
         WRITE(numout,cform_war)
         IF( PRESENT(cd1 ) ) WRITE(numout,*) cd1
         IF( PRESENT(cd2 ) ) WRITE(numout,*) cd2
         IF( PRESENT(cd3 ) ) WRITE(numout,*) cd3
         IF( PRESENT(cd4 ) ) WRITE(numout,*) cd4
         IF( PRESENT(cd5 ) ) WRITE(numout,*) cd5
         IF( PRESENT(cd6 ) ) WRITE(numout,*) cd6
         IF( PRESENT(cd7 ) ) WRITE(numout,*) cd7
         IF( PRESENT(cd8 ) ) WRITE(numout,*) cd8
         IF( PRESENT(cd9 ) ) WRITE(numout,*) cd9
         IF( PRESENT(cd10) ) WRITE(numout,*) cd10
      ENDIF
      CALL FLUSH(numout)
      !
   END SUBROUTINE ctl_warn


   SUBROUTINE ctl_opn( knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_opn  ***
      !!
      !! ** Purpose :   Open file and check if required file is available.
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER          , INTENT(  out) ::   knum      ! logical unit to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdfile    ! file name to open
      CHARACTER(len=*) , INTENT(in   ) ::   cdstat    ! disposition specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdform    ! formatting specifier
      CHARACTER(len=*) , INTENT(in   ) ::   cdacce    ! access specifier
      INTEGER          , INTENT(in   ) ::   klengh    ! record length
      INTEGER          , INTENT(in   ) ::   kout      ! number of logical units for write
      LOGICAL          , INTENT(in   ) ::   ldwp      ! boolean term for print
      INTEGER, OPTIONAL, INTENT(in   ) ::   karea     ! proc number
      !
      CHARACTER(len=80) ::   clfile
      INTEGER           ::   iost
      !!----------------------------------------------------------------------
      !
      ! adapt filename
      ! ----------------
      clfile = TRIM(cdfile)
      IF( PRESENT( karea ) ) THEN
         IF( karea > 1 )   WRITE(clfile, "(a,'_',i4.4)") TRIM(clfile), karea-1
      ENDIF
#if defined key_agrif
      IF( .NOT. Agrif_Root() )   clfile = TRIM(Agrif_CFixed())//'_'//TRIM(clfile)
      knum=Agrif_Get_Unit()
#else
      knum=get_unit()
#endif
      !
      iost=0
      IF( cdacce(1:6) == 'DIRECT' )  THEN
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
      ELSE
         OPEN( UNIT=knum, FILE=clfile, FORM=cdform, ACCESS=cdacce, STATUS=cdstat             , ERR=100, IOSTAT=iost )
      ENDIF
      IF( iost == 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*) '     file   : ', clfile,' open ok'
            WRITE(kout,*) '     unit   = ', knum
            WRITE(kout,*) '     status = ', cdstat
            WRITE(kout,*) '     form   = ', cdform
            WRITE(kout,*) '     access = ', cdacce
            WRITE(kout,*)
         ENDIF
      ENDIF
100   CONTINUE
      IF( iost /= 0 ) THEN
         IF(ldwp) THEN
            WRITE(kout,*)
            WRITE(kout,*) ' ===>>>> : bad opening file: ', clfile
            WRITE(kout,*) ' =======   ===  '
            WRITE(kout,*) '           unit   = ', knum
            WRITE(kout,*) '           status = ', cdstat
            WRITE(kout,*) '           form   = ', cdform
            WRITE(kout,*) '           access = ', cdacce
            WRITE(kout,*) '           iostat = ', iost
            WRITE(kout,*) '           we stop. verify the file '
            WRITE(kout,*)
         ELSE  !!! Force writing to make sure we get the information - at least once - in this violent STOP!!
            WRITE(*,*)
            WRITE(*,*) ' ===>>>> : bad opening file: ', clfile
            WRITE(*,*) ' =======   ===  '
            WRITE(*,*) '           unit   = ', knum
            WRITE(*,*) '           status = ', cdstat
            WRITE(*,*) '           form   = ', cdform
            WRITE(*,*) '           access = ', cdacce
            WRITE(*,*) '           iostat = ', iost
            WRITE(*,*) '           we stop. verify the file '
            WRITE(*,*)
         ENDIF
         CALL FLUSH( kout ) 
         STOP 'ctl_opn bad opening'
      ENDIF
      !
   END SUBROUTINE ctl_opn


   SUBROUTINE ctl_nam ( kios, cdnam, ldwp )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ctl_nam  ***
      !!
      !! ** Purpose :   Informations when error while reading a namelist
      !!
      !! ** Method  :   Fortan open
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(inout) ::   kios    ! IO status after reading the namelist
      CHARACTER(len=*), INTENT(in   ) ::   cdnam   ! group name of namelist for which error occurs
      CHARACTER(len=5)                ::   clios   ! string to convert iostat in character for print
      LOGICAL         , INTENT(in   ) ::   ldwp    ! boolean term for print
      !!----------------------------------------------------------------------
      !
      WRITE (clios, '(I5.0)')   kios
      IF( kios < 0 ) THEN         
         CALL ctl_warn( 'end of record or file while reading namelist '   &
            &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      !
      IF( kios > 0 ) THEN
         CALL ctl_stop( 'misspelled variable in namelist '   &
            &           // TRIM(cdnam) // ' iostat = ' // TRIM(clios) )
      ENDIF
      kios = 0
      RETURN
      !
   END SUBROUTINE ctl_nam


   INTEGER FUNCTION get_unit()
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  get_unit  ***
      !!
      !! ** Purpose :   return the index of an unused logical unit
      !!----------------------------------------------------------------------
      LOGICAL :: llopn
      !!----------------------------------------------------------------------
      !
      get_unit = 15   ! choose a unit that is big enough then it is not already used in NEMO
      llopn = .TRUE.
      DO WHILE( (get_unit < 998) .AND. llopn )
         get_unit = get_unit + 1
         INQUIRE( unit = get_unit, opened = llopn )
      END DO
      IF( (get_unit == 999) .AND. llopn ) THEN
         CALL ctl_stop( 'get_unit: All logical units until 999 are used...' )
         get_unit = -1
      ENDIF
      !
   END FUNCTION get_unit

   !!----------------------------------------------------------------------
END MODULE lib_mpp
