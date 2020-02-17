MODULE lib_fortran
   !!======================================================================
   !!                       ***  MODULE  lib_fortran  ***
   !! Fortran utilities:  includes some low levels fortran functionality
   !!======================================================================
   !! History :  3.2  !  2010-05  (M. Dunphy, R. Benshila)  Original code
   !!            3.4  !  2013-06  (C. Rousset)  add glob_min, glob_max 
   !!                                           + 3d dim. of input is fexible (jpk, jpl...) 
   !!            4.0  !  2016-06  (T. Lovato)  double precision global sum by default 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   glob_sum    : generic interface for global masked summation over
   !!                 the interior domain for 1 or 2 2D or 3D arrays
   !!                 it works only for T points
   !!   SIGN        : generic interface for SIGN to overwrite f95 behaviour
   !!                 of intrinsinc sign function
   !!----------------------------------------------------------------------
   USE par_oce         ! Ocean parameter
   USE dom_oce         ! ocean domain
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   glob_sum      ! used in many places (masked with tmask_i)
   PUBLIC   glob_sum_full ! used in many places (masked with tmask_h, ie only over the halos)
   PUBLIC   DDPDD         ! also used in closea module
   PUBLIC   glob_min, glob_max
#if defined key_nosignedzero
   PUBLIC SIGN
#endif

   INTERFACE glob_sum
      MODULE PROCEDURE glob_sum_1d, glob_sum_2d, glob_sum_3d, &
         &             glob_sum_2d_a, glob_sum_3d_a
   END INTERFACE
   INTERFACE glob_sum_full
      MODULE PROCEDURE glob_sum_full_2d, glob_sum_full_3d
   END INTERFACE
   INTERFACE glob_min
      MODULE PROCEDURE glob_min_2d, glob_min_3d,glob_min_2d_a, glob_min_3d_a 
   END INTERFACE
   INTERFACE glob_max
      MODULE PROCEDURE glob_max_2d, glob_max_3d,glob_max_2d_a, glob_max_3d_a 
   END INTERFACE

#if defined key_nosignedzero
   INTERFACE SIGN
      MODULE PROCEDURE SIGN_SCALAR, SIGN_ARRAY_1D, SIGN_ARRAY_2D, SIGN_ARRAY_3D,   &
         &             SIGN_ARRAY_1D_A, SIGN_ARRAY_2D_A, SIGN_ARRAY_3D_A,          &
         &             SIGN_ARRAY_1D_B, SIGN_ARRAY_2D_B, SIGN_ARRAY_3D_B
   END INTERFACE
#endif

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lib_fortran.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   ! --- SUM ---
   FUNCTION glob_sum_1d( ptab, kdim )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_1d ***
      !!
      !! ** Purpose : perform a sum in calling DDPDD routine
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in) :: kdim
      REAL(wp), INTENT(in), DIMENSION(kdim) ::   ptab
      REAL(wp)                              ::   glob_sum_1d   ! global sum
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO ji = 1, kdim
         ztmp =  ptab(ji)
         CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_1d = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_1d

   FUNCTION glob_sum_2d( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_2d ***
      !!
      !! ** Purpose : perform a sum in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab
      REAL(wp)                             ::   glob_sum_2d   ! global masked sum
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jj = 1, jpj
         DO ji =1, jpi
         ztmp =  ptab(ji,jj) * tmask_i(ji,jj)
         CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_2d = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_2d


   FUNCTION glob_sum_3d( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_3d ***
      !!
      !! ** Purpose : perform a sum on a 3D array in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab
      REAL(wp)                               ::   glob_sum_3d   ! global masked sum
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj, jk   ! dummy loop indices
      INTEGER    ::   ijpk ! local variables: size of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab,3)
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jk = 1, ijpk
         DO jj = 1, jpj
            DO ji =1, jpi
            ztmp =  ptab(ji,jj,jk) * tmask_i(ji,jj)
            CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_3d = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_3d


   FUNCTION glob_sum_2d_a( ptab1, ptab2 )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_2d_a ***
      !!
      !! ** Purpose : perform a sum on two 2D arrays in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab1, ptab2
      REAL(wp)                             ::   glob_sum_2d_a   ! global masked sum
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jj = 1, jpj
         DO ji =1, jpi
         ztmp =  ptab1(ji,jj) * tmask_i(ji,jj)
         CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         ztmp =  ptab2(ji,jj) * tmask_i(ji,jj)
         CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_2d_a = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_2d_a


   FUNCTION glob_sum_3d_a( ptab1, ptab2 )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_3d_a ***
      !!
      !! ** Purpose : perform a sum on two 3D array in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab1, ptab2
      REAL(wp)                               ::   glob_sum_3d_a   ! global masked sum
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj, jk   ! dummy loop indices
      INTEGER    ::   ijpk ! local variables: size of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab1,3)
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jk = 1, ijpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmp =  ptab1(ji,jj,jk) * tmask_i(ji,jj)
               CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
               ztmp =  ptab2(ji,jj,jk) * tmask_i(ji,jj)
               CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
            END DO
         END DO    
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_3d_a = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_3d_a   

   FUNCTION glob_sum_full_2d( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_full_2d ***
      !!
      !! ** Purpose : perform a sum in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab
      REAL(wp)                             ::   glob_sum_full_2d   ! global sum (nomask)
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jj = 1, jpj
         DO ji =1, jpi
         ztmp =  ptab(ji,jj) * tmask_h(ji,jj) 
         CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_full_2d = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_full_2d

   FUNCTION glob_sum_full_3d( ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_full_3d ***
      !!
      !! ** Purpose : perform a sum on a 3D array in calling DDPDD routine
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab
      REAL(wp)                               ::   glob_sum_full_3d   ! global sum (nomask)
      !!
      COMPLEX(wp)::   ctmp
      REAL(wp)   ::   ztmp
      INTEGER    ::   ji, jj, jk   ! dummy loop indices
      INTEGER    ::   ijpk ! local variables: size of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab,3)
      !
      ztmp = 0.e0
      ctmp = CMPLX( 0.e0, 0.e0, wp )
      DO jk = 1, ijpk
         DO jj = 1, jpj
            DO ji =1, jpi
            ztmp =  ptab(ji,jj,jk) * tmask_h(ji,jj)
            CALL DDPDD( CMPLX( ztmp, 0.e0, wp ), ctmp )
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( ctmp )   ! sum over the global domain
      glob_sum_full_3d = REAL(ctmp,wp)
      !
   END FUNCTION glob_sum_full_3d

   ! --- MIN ---
   FUNCTION glob_min_2d( ptab ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_min_2D  ***
      !!
      !! ** Purpose : perform a masked min on the inner global domain of a 2D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab          ! input 2D array
      REAL(wp)                             ::   glob_min_2d   ! global masked min
      !!-----------------------------------------------------------------------
      !
      glob_min_2d = MINVAL( ptab(:,:)*tmask_i(:,:) )
      IF( lk_mpp )   CALL mpp_min( glob_min_2d )
      !
   END FUNCTION glob_min_2d
 
   FUNCTION glob_min_3d( ptab ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_min_3D  ***
      !!
      !! ** Purpose : perform a masked min on the inner global domain of a 3D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab          ! input 3D array
      REAL(wp)                               ::   glob_min_3d   ! global masked min
      !!
      INTEGER :: jk
      INTEGER :: ijpk ! local variable: size of the 3d dimension of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab,3)
      !
      glob_min_3d = MINVAL( ptab(:,:,1)*tmask_i(:,:) )
      DO jk = 2, ijpk
         glob_min_3d = MIN( glob_min_3d, MINVAL( ptab(:,:,jk)*tmask_i(:,:) ) )
      END DO
      IF( lk_mpp )   CALL mpp_min( glob_min_3d )
      !
   END FUNCTION glob_min_3d


   FUNCTION glob_min_2d_a( ptab1, ptab2 ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_min_2D _a ***
      !!
      !! ** Purpose : perform a masked min on the inner global domain of two 2D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab1, ptab2    ! input 2D array
      REAL(wp)            , DIMENSION(2)   ::   glob_min_2d_a   ! global masked min
      !!-----------------------------------------------------------------------
      !             
      glob_min_2d_a(1) = MINVAL( ptab1(:,:)*tmask_i(:,:) )
      glob_min_2d_a(2) = MINVAL( ptab2(:,:)*tmask_i(:,:) )
      IF( lk_mpp )   CALL mpp_min( glob_min_2d_a, 2 )
      !
   END FUNCTION glob_min_2d_a
 
 
   FUNCTION glob_min_3d_a( ptab1, ptab2 ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_min_3D_a ***
      !!
      !! ** Purpose : perform a masked min on the inner global domain of two 3D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab1, ptab2    ! input 3D array
      REAL(wp)            , DIMENSION(2)     ::   glob_min_3d_a   ! global masked min
      !!
      INTEGER :: jk
      INTEGER :: ijpk ! local variable: size of the 3d dimension of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab1,3)
      !
      glob_min_3d_a(1) = MINVAL( ptab1(:,:,1)*tmask_i(:,:) )
      glob_min_3d_a(2) = MINVAL( ptab2(:,:,1)*tmask_i(:,:) )
      DO jk = 2, ijpk
         glob_min_3d_a(1) = MIN( glob_min_3d_a(1), MINVAL( ptab1(:,:,jk)*tmask_i(:,:) ) )
         glob_min_3d_a(2) = MIN( glob_min_3d_a(2), MINVAL( ptab2(:,:,jk)*tmask_i(:,:) ) )
      END DO
      IF( lk_mpp )   CALL mpp_min( glob_min_3d_a, 2 )
      !
   END FUNCTION glob_min_3d_a

   ! --- MAX ---
   FUNCTION glob_max_2d( ptab ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_max_2D  ***
      !!
      !! ** Purpose : perform a masked max on the inner global domain of a 2D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab          ! input 2D array
      REAL(wp)                             ::   glob_max_2d   ! global masked max
      !!-----------------------------------------------------------------------
      !
      glob_max_2d = MAXVAL( ptab(:,:)*tmask_i(:,:) )
      IF( lk_mpp )   CALL mpp_max( glob_max_2d )
      !
   END FUNCTION glob_max_2d
 
   FUNCTION glob_max_3d( ptab ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_max_3D  ***
      !!
      !! ** Purpose : perform a masked max on the inner global domain of a 3D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab          ! input 3D array
      REAL(wp)                               ::   glob_max_3d   ! global masked max
      !!
      INTEGER :: jk
      INTEGER :: ijpk ! local variable: size of the 3d dimension of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab,3)
      !
      glob_max_3d = MAXVAL( ptab(:,:,1)*tmask_i(:,:) )
      DO jk = 2, ijpk
         glob_max_3d = MAX( glob_max_3d, MAXVAL( ptab(:,:,jk)*tmask_i(:,:) ) )
      END DO
      IF( lk_mpp )   CALL mpp_max( glob_max_3d )
      !
   END FUNCTION glob_max_3d


   FUNCTION glob_max_2d_a( ptab1, ptab2 ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_max_2D _a ***
      !!
      !! ** Purpose : perform a masked max on the inner global domain of two 2D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptab1, ptab2    ! input 2D array
      REAL(wp)            , DIMENSION(2)   ::   glob_max_2d_a   ! global masked max
      !!-----------------------------------------------------------------------
      !             
      glob_max_2d_a(1) = MAXVAL( ptab1(:,:)*tmask_i(:,:) )
      glob_max_2d_a(2) = MAXVAL( ptab2(:,:)*tmask_i(:,:) )
      IF( lk_mpp )   CALL mpp_max( glob_max_2d_a, 2 )
      !
   END FUNCTION glob_max_2d_a
 
 
   FUNCTION glob_max_3d_a( ptab1, ptab2 ) 
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_max_3D_a ***
      !!
      !! ** Purpose : perform a masked max on the inner global domain of two 3D array
      !!-----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:,:) ::   ptab1, ptab2    ! input 3D array
      REAL(wp)            , DIMENSION(2)     ::   glob_max_3d_a   ! global masked max
      !!
      INTEGER :: jk
      INTEGER :: ijpk ! local variable: size of the 3d dimension of ptab
      !!-----------------------------------------------------------------------
      !
      ijpk = SIZE(ptab1,3)
      !
      glob_max_3d_a(1) = MAXVAL( ptab1(:,:,1)*tmask_i(:,:) )
      glob_max_3d_a(2) = MAXVAL( ptab2(:,:,1)*tmask_i(:,:) )
      DO jk = 2, ijpk
         glob_max_3d_a(1) = MAX( glob_max_3d_a(1), MAXVAL( ptab1(:,:,jk)*tmask_i(:,:) ) )
         glob_max_3d_a(2) = MAX( glob_max_3d_a(2), MAXVAL( ptab2(:,:,jk)*tmask_i(:,:) ) )
      END DO
      IF( lk_mpp )   CALL mpp_max( glob_max_3d_a, 2 )
      !
   END FUNCTION glob_max_3d_a


   SUBROUTINE DDPDD( ydda, yddb )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE DDPDD ***
      !!
      !! ** Purpose : Add a scalar element to a sum
      !!
      !!
      !! ** Method  : The code uses the compensated summation with doublet
      !!              (sum,error) emulated useing complex numbers. ydda is the
      !!               scalar to add to the summ yddb
      !!
      !! ** Action  : This does only work for MPI.
      !!
      !! References : Using Acurate Arithmetics to Improve Numerical
      !!              Reproducibility and Sability in Parallel Applications
      !!              Yun HE and Chris H. Q. DING, Journal of Supercomputing 18, 259-277, 2001
      !!----------------------------------------------------------------------
      COMPLEX(wp), INTENT(in   ) ::   ydda
      COMPLEX(wp), INTENT(inout) ::   yddb
      !
      REAL(wp) :: zerr, zt1, zt2  ! local work variables
      !!-----------------------------------------------------------------------
      !
      ! Compute ydda + yddb using Knuth's trick.
      zt1  = REAL(ydda) + REAL(yddb)
      zerr = zt1 - REAL(ydda)
      zt2  = ( (REAL(yddb) - zerr) + (REAL(ydda) - (zt1 - zerr)) )   &
         &   + AIMAG(ydda)         + AIMAG(yddb)
      !
      ! The result is t1 + t2, after normalization.
      yddb = CMPLX( zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp )
      !
   END SUBROUTINE DDPDD

#if defined key_nosignedzero
   !!----------------------------------------------------------------------
   !!   'key_nosignedzero'                                         F90 SIGN
   !!----------------------------------------------------------------------

   FUNCTION SIGN_SCALAR( pa, pb )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_SCALAR  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb          ! input
      REAL(wp) :: SIGN_SCALAR    ! result
      !!-----------------------------------------------------------------------
      IF ( pb >= 0.e0) THEN   ;   SIGN_SCALAR = ABS(pa)
      ELSE                    ;   SIGN_SCALAR =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_SCALAR


   FUNCTION SIGN_ARRAY_1D( pa, pb )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:)                   ! input
      REAL(wp) :: SIGN_ARRAY_1D(SIZE(pb,1))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_1D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_1D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_1D


   FUNCTION SIGN_ARRAY_2D(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_2D(SIZE(pb,1),SIZE(pb,2))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_2D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_2D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_2D

   FUNCTION SIGN_ARRAY_3D(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa,pb(:,:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_3D(SIZE(pb,1),SIZE(pb,2),SIZE(pb,3))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_3D = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_3D =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_3D


   FUNCTION SIGN_ARRAY_1D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:),pb(:)      ! input
      REAL(wp) :: SIGN_ARRAY_1D_A(SIZE(pb,1))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_1D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_1D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_1D_A


   FUNCTION SIGN_ARRAY_2D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:),pb(:,:)      ! input
      REAL(wp) :: SIGN_ARRAY_2D_A(SIZE(pb,1),SIZE(pb,2))  ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_2D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_2D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_2D_A


   FUNCTION SIGN_ARRAY_3D_A(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D_A  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:,:),pb(:,:,:)  ! input
      REAL(wp) :: SIGN_ARRAY_3D_A(SIZE(pb,1),SIZE(pb,2),SIZE(pb,3)) ! result
      !!-----------------------------------------------------------------------
      WHERE ( pb >= 0.e0 )   ;   SIGN_ARRAY_3D_A = ABS(pa)
      ELSEWHERE              ;   SIGN_ARRAY_3D_A =-ABS(pa)
      END WHERE
   END FUNCTION SIGN_ARRAY_3D_A


   FUNCTION SIGN_ARRAY_1D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_1D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_1D_B(SIZE(pa,1))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_1D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_1D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_1D_B


   FUNCTION SIGN_ARRAY_2D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_2D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_2D_B(SIZE(pa,1),SIZE(pa,2))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_2D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_2D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_2D_B


   FUNCTION SIGN_ARRAY_3D_B(pa,pb)
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION SIGN_ARRAY_3D_B  ***
      !!
      !! ** Purpose : overwrite f95 behaviour of intrinsinc sign function
      !!-----------------------------------------------------------------------
      REAL(wp) :: pa(:,:,:),pb      ! input
      REAL(wp) :: SIGN_ARRAY_3D_B(SIZE(pa,1),SIZE(pa,2),SIZE(pa,3))  ! result
      !!-----------------------------------------------------------------------
      IF( pb >= 0.e0 ) THEN   ;   SIGN_ARRAY_3D_B = ABS(pa)
      ELSE                    ;   SIGN_ARRAY_3D_B =-ABS(pa)
      ENDIF
   END FUNCTION SIGN_ARRAY_3D_B
#endif

   !!======================================================================
END MODULE lib_fortran
