MODULE diaptr
   !!======================================================================
   !!                       ***  MODULE  diaptr  ***
   !! Ocean physics:  Computes meridonal transports and zonal means
   !!=====================================================================
   !! History :  1.0  ! 2003-09  (C. Talandier, G. Madec)  Original code
   !!            2.0  ! 2006-01  (A. Biastoch)  Allow sub-basins computation
   !!            3.2  ! 2010-03  (O. Marti, S. Flavoni) Add fields
   !!            3.3  ! 2010-10  (G. Madec)  dynamical allocation
   !!            3.6  ! 2014-12  (C. Ethe) use of IOM
   !!            3.6  ! 2016-06  (T. Graham) Addition of diagnostics for CMIP6
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_ptr      : Poleward Transport Diagnostics module
   !!   dia_ptr_init : Initialization, namelist read
   !!   ptr_sjk      : "zonal" mean computation of a field - tracer or flux array
   !!   ptr_sj       : "zonal" and vertical sum computation of a "meridional" flux array
   !!                   (Generic interface to ptr_sj_3d, ptr_sj_2d)
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and active tracers
   USE dom_oce          ! ocean space and time domain
   USE phycst           ! physical constants
   !
   USE iom              ! IOM library
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! MPP library
   USE timing           ! preformance summary

   IMPLICIT NONE
   PRIVATE

   INTERFACE ptr_sj
      MODULE PROCEDURE ptr_sj_3d, ptr_sj_2d
   END INTERFACE

   PUBLIC   ptr_sj         ! call by tra_ldf & tra_adv routines
   PUBLIC   ptr_sjk        ! 
   PUBLIC   dia_ptr_init   ! call in memogcm
   PUBLIC   dia_ptr        ! call in step module
   PUBLIC   dia_ptr_hst    ! called from tra_ldf/tra_adv routines

   !                                  !!** namelist  namptr  **
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   htr_adv, htr_ldf, htr_eiv   !: Heat TRansports (adv, diff, Bolus.)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   str_adv, str_ldf, str_eiv   !: Salt TRansports (adv, diff, Bolus.)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   htr_ove, str_ove   !: heat Salt TRansports ( overturn.)
   REAL(wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:,:) ::   htr_btr, str_btr   !: heat Salt TRansports ( barotropic )

   LOGICAL, PUBLIC ::   ln_diaptr   !  Poleward transport flag (T) or not (F)
   LOGICAL, PUBLIC ::   ln_subbas   !  Atlantic/Pacific/Indian basins calculation
   INTEGER, PUBLIC ::   nptr        ! = 1 (l_subbas=F) or = 5 (glo, atl, pac, ind, ipc) (l_subbas=T) 

   REAL(wp) ::   rc_sv    = 1.e-6_wp   ! conversion from m3/s to Sverdrup
   REAL(wp) ::   rc_pwatt = 1.e-15_wp  ! conversion from W    to PW (further x rau0 x Cp)
   REAL(wp) ::   rc_ggram = 1.e-6_wp   ! conversion from g    to Pg

   CHARACTER(len=3), ALLOCATABLE, SAVE, DIMENSION(:)     :: clsubb
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: btmsk   ! T-point basin interior masks
   REAL(wp),         ALLOCATABLE, SAVE, DIMENSION(:,:)   :: btm30   ! mask out Southern Ocean (=0 south of 30°S)

   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:)     :: p_fval1d
   REAL(wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: p_fval2d

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diaptr.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_ptr( pvtr )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr  ***
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::   pvtr   ! j-effective transport
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zsfc,zvfc               ! local scalar
      REAL(wp), DIMENSION(jpi,jpj)     ::  z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  z3d   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zmask   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts) ::  zts   ! 3D workspace
      REAL(wp), DIMENSION(jpj)     ::  vsum   ! 1D workspace
      REAL(wp), DIMENSION(jpj,jpts)     ::  tssum   ! 1D workspace
 
      !
      !overturning calculation
      REAL(wp), DIMENSION(jpj,jpk,nptr) ::   sjk  , r1_sjk ! i-mean i-k-surface and its inverse
      REAL(wp), DIMENSION(jpj,jpk,nptr) ::   v_msf, sn_jk  , tn_jk ! i-mean T and S, j-Stream-Function
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zvn   ! 3D workspace


      CHARACTER( len = 12 )  :: cl1
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dia_ptr')

      !
      IF( PRESENT( pvtr ) ) THEN
         IF( iom_use("zomsfglo") ) THEN    ! effective MSF
            z3d(1,:,:) = ptr_sjk( pvtr(:,:,:) )  ! zonal cumulative effective transport
            DO jk = 2, jpkm1 
              z3d(1,:,jk) = z3d(1,:,jk-1) + z3d(1,:,jk)   ! effective j-Stream-Function (MSF)
            END DO
            DO ji = 1, jpi
               z3d(ji,:,:) = z3d(1,:,:)
            ENDDO
            cl1 = TRIM('zomsf'//clsubb(1) )
            CALL iom_put( cl1, z3d * rc_sv )
            DO jn = 2, nptr                                    ! by sub-basins
               z3d(1,:,:) =  ptr_sjk( pvtr(:,:,:), btmsk(:,:,jn)*btm30(:,:) ) 
               DO jk = 2, jpkm1 
                  z3d(1,:,jk) = z3d(1,:,jk-1) + z3d(1,:,jk)    ! effective j-Stream-Function (MSF)
               END DO
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zomsf'//clsubb(jn) )
               CALL iom_put( cl1, z3d * rc_sv )
            END DO
         ENDIF
         IF( iom_use("sopstove") .OR. iom_use("sophtove") .OR. iom_use("sopstbtr") .OR. iom_use("sophtbtr") ) THEN
            ! define fields multiplied by scalar
            zmask(:,:,:) = 0._wp
            zts(:,:,:,:) = 0._wp
            zvn(:,:,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1
                  DO ji = 1, jpi
                     zvfc = e1v(ji,jj) * e3v_n(ji,jj,jk)
                     zmask(ji,jj,jk)      = vmask(ji,jj,jk)      * zvfc
                     zts(ji,jj,jk,jp_tem) = (tsn(ji,jj,jk,jp_tem)+tsn(ji,jj+1,jk,jp_tem)) * 0.5 * zvfc  !Tracers averaged onto V grid
                     zts(ji,jj,jk,jp_sal) = (tsn(ji,jj,jk,jp_sal)+tsn(ji,jj+1,jk,jp_sal)) * 0.5 * zvfc
                     zvn(ji,jj,jk)        = vn(ji,jj,jk)         * zvfc
                  ENDDO
               ENDDO
             ENDDO
         ENDIF
         IF( iom_use("sopstove") .OR. iom_use("sophtove") ) THEN
             sjk(:,:,1) = ptr_sjk( zmask(:,:,:), btmsk(:,:,1) )
             r1_sjk(:,:,1) = 0._wp
             WHERE( sjk(:,:,1) /= 0._wp )   r1_sjk(:,:,1) = 1._wp / sjk(:,:,1)

             ! i-mean T and S, j-Stream-Function, global
             tn_jk(:,:,1) = ptr_sjk( zts(:,:,:,jp_tem) ) * r1_sjk(:,:,1)
             sn_jk(:,:,1) = ptr_sjk( zts(:,:,:,jp_sal) ) * r1_sjk(:,:,1)
             v_msf(:,:,1) = ptr_sjk( zvn(:,:,:) )

             htr_ove(:,1) = SUM( v_msf(:,:,1)*tn_jk(:,:,1) ,2 )
             str_ove(:,1) = SUM( v_msf(:,:,1)*sn_jk(:,:,1) ,2 )

             z2d(1,:) = htr_ove(:,1) * rc_pwatt        !  (conversion in PW)
             DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
             ENDDO
             cl1 = 'sophtove'
             CALL iom_put( TRIM(cl1), z2d )
             z2d(1,:) = str_ove(:,1) * rc_ggram        !  (conversion in Gg)
             DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
             ENDDO
             cl1 = 'sopstove'
             CALL iom_put( TRIM(cl1), z2d )
             IF( ln_subbas ) THEN
                DO jn = 2, nptr
                    sjk(:,:,jn) = ptr_sjk( zmask(:,:,:), btmsk(:,:,jn) )
                    r1_sjk(:,:,jn) = 0._wp
                    WHERE( sjk(:,:,jn) /= 0._wp )   r1_sjk(:,:,jn) = 1._wp / sjk(:,:,jn)

                    ! i-mean T and S, j-Stream-Function, basin
                    tn_jk(:,:,jn) = ptr_sjk( zts(:,:,:,jp_tem), btmsk(:,:,jn) ) * r1_sjk(:,:,jn)
                    sn_jk(:,:,jn) = ptr_sjk( zts(:,:,:,jp_sal), btmsk(:,:,jn) ) * r1_sjk(:,:,jn)
                    v_msf(:,:,jn) = ptr_sjk( zvn(:,:,:), btmsk(:,:,jn) ) 
                    htr_ove(:,jn) = SUM( v_msf(:,:,jn)*tn_jk(:,:,jn) ,2 )
                    str_ove(:,jn) = SUM( v_msf(:,:,jn)*sn_jk(:,:,jn) ,2 )

                    z2d(1,:) = htr_ove(:,jn) * rc_pwatt !  (conversion in PW)
                    DO ji = 1, jpi
                        z2d(ji,:) = z2d(1,:)
                    ENDDO
                    cl1 = TRIM('sophtove_'//clsubb(jn))
                    CALL iom_put( cl1, z2d )
                    z2d(1,:) = str_ove(:,jn) * rc_ggram        ! (conversion in Gg)
                    DO ji = 1, jpi
                        z2d(ji,:) = z2d(1,:)
                    ENDDO
                    cl1 = TRIM('sopstove_'//clsubb(jn))
                    CALL iom_put( cl1, z2d )
                END DO
             ENDIF
         ENDIF
         IF( iom_use("sopstbtr") .OR. iom_use("sophtbtr") ) THEN
         ! Calculate barotropic heat and salt transport here 
             sjk(:,1,1) = ptr_sj( zmask(:,:,:), btmsk(:,:,1) )
             r1_sjk(:,1,1) = 0._wp
             WHERE( sjk(:,1,1) /= 0._wp )   r1_sjk(:,1,1) = 1._wp / sjk(:,1,1)
            
            vsum = ptr_sj( zvn(:,:,:), btmsk(:,:,1))
            tssum(:,jp_tem) = ptr_sj( zts(:,:,:,jp_tem), btmsk(:,:,1) )
            tssum(:,jp_sal) = ptr_sj( zts(:,:,:,jp_sal), btmsk(:,:,1) )
            htr_btr(:,1) = vsum * tssum(:,jp_tem) * r1_sjk(:,1,1)
            str_btr(:,1) = vsum * tssum(:,jp_sal) * r1_sjk(:,1,1)
            z2d(1,:) = htr_btr(:,1) * rc_pwatt        !  (conversion in PW)
            DO ji = 2, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophtbtr'
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_btr(:,1) * rc_ggram        !  (conversion in Gg)
            DO ji = 2, jpi
              z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopstbtr'
            CALL iom_put( TRIM(cl1), z2d )
            IF( ln_subbas ) THEN
                DO jn = 2, nptr
                    sjk(:,1,jn) = ptr_sj( zmask(:,:,:), btmsk(:,:,jn) )
                    r1_sjk(:,1,jn) = 0._wp
                    WHERE( sjk(:,1,jn) /= 0._wp )   r1_sjk(:,1,jn) = 1._wp / sjk(:,1,jn)
                    vsum = ptr_sj( zvn(:,:,:), btmsk(:,:,jn))
                    tssum(:,jp_tem) = ptr_sj( zts(:,:,:,jp_tem), btmsk(:,:,jn) )
                    tssum(:,jp_sal) = ptr_sj( zts(:,:,:,jp_sal), btmsk(:,:,jn) )
                    htr_btr(:,jn) = vsum * tssum(:,jp_tem) * r1_sjk(:,1,jn)
                    str_btr(:,jn) = vsum * tssum(:,jp_sal) * r1_sjk(:,1,jn)
                    z2d(1,:) = htr_btr(:,jn) * rc_pwatt !  (conversion in PW)
                    DO ji = 1, jpi
                        z2d(ji,:) = z2d(1,:)
                    ENDDO
                    cl1 = TRIM('sophtbtr_'//clsubb(jn))
                    CALL iom_put( cl1, z2d )
                    z2d(1,:) = str_btr(:,jn) * rc_ggram        ! (conversion in Gg)
                    DO ji = 1, jpi
                        z2d(ji,:) = z2d(1,:)
                    ENDDO
                    cl1 = TRIM('sopstbtr_'//clsubb(jn))
                    CALL iom_put( cl1, z2d )
               ENDDO
            ENDIF !ln_subbas
         ENDIF !iom_use("sopstbtr....)
         !
      ELSE
         !
         IF( iom_use("zotemglo") ) THEN    ! i-mean i-k-surface 
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zsfc = e1t(ji,jj) * e3t_n(ji,jj,jk)
                     zmask(ji,jj,jk)      = tmask(ji,jj,jk)      * zsfc
                     zts(ji,jj,jk,jp_tem) = tsn(ji,jj,jk,jp_tem) * zsfc
                     zts(ji,jj,jk,jp_sal) = tsn(ji,jj,jk,jp_sal) * zsfc
                  END DO
               END DO
            END DO
            DO jn = 1, nptr
               zmask(1,:,:) = ptr_sjk( zmask(:,:,:), btmsk(:,:,jn) )
               cl1 = TRIM('zosrf'//clsubb(jn) )
               CALL iom_put( cl1, zmask )
               !
               z3d(1,:,:) = ptr_sjk( zts(:,:,:,jp_tem), btmsk(:,:,jn) ) &
                  &            / MAX( zmask(1,:,:), 10.e-15 )
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zotem'//clsubb(jn) )
               CALL iom_put( cl1, z3d )
               !
               z3d(1,:,:) = ptr_sjk( zts(:,:,:,jp_sal), btmsk(:,:,jn) ) &
                  &            / MAX( zmask(1,:,:), 10.e-15 )
               DO ji = 1, jpi
                  z3d(ji,:,:) = z3d(1,:,:)
               ENDDO
               cl1 = TRIM('zosal'//clsubb(jn) )
               CALL iom_put( cl1, z3d )
            END DO
         ENDIF
         !
         !                                ! Advective and diffusive heat and salt transport
         IF( iom_use("sophtadv") .OR. iom_use("sopstadv") ) THEN   
            z2d(1,:) = htr_adv(:,1) * rc_pwatt        !  (conversion in PW)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophtadv'                 
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_adv(:,1) * rc_ggram        ! (conversion in Gg)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopstadv'
            CALL iom_put( TRIM(cl1), z2d )
            IF( ln_subbas ) THEN
              DO jn=2,nptr
               z2d(1,:) = htr_adv(:,jn) * rc_pwatt        !  (conversion in PW)
               DO ji = 1, jpi
                 z2d(ji,:) = z2d(1,:)
               ENDDO
               cl1 = TRIM('sophtadv_'//clsubb(jn))                 
               CALL iom_put( cl1, z2d )
               z2d(1,:) = str_adv(:,jn) * rc_ggram        ! (conversion in Gg)
               DO ji = 1, jpi
                  z2d(ji,:) = z2d(1,:)
               ENDDO
               cl1 = TRIM('sopstadv_'//clsubb(jn))                 
               CALL iom_put( cl1, z2d )              
              ENDDO
            ENDIF
         ENDIF
         !
         IF( iom_use("sophtldf") .OR. iom_use("sopstldf") ) THEN   
            z2d(1,:) = htr_ldf(:,1) * rc_pwatt        !  (conversion in PW) 
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophtldf'
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_ldf(:,1) * rc_ggram        !  (conversion in Gg)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopstldf'
            CALL iom_put( TRIM(cl1), z2d )
            IF( ln_subbas ) THEN
              DO jn=2,nptr
               z2d(1,:) = htr_ldf(:,jn) * rc_pwatt        !  (conversion in PW)
               DO ji = 1, jpi
                 z2d(ji,:) = z2d(1,:)
               ENDDO
               cl1 = TRIM('sophtldf_'//clsubb(jn))                 
               CALL iom_put( cl1, z2d )
               z2d(1,:) = str_ldf(:,jn) * rc_ggram        ! (conversion in Gg)
               DO ji = 1, jpi
                  z2d(ji,:) = z2d(1,:)
               ENDDO
               cl1 = TRIM('sopstldf_'//clsubb(jn))                 
               CALL iom_put( cl1, z2d )              
              ENDDO
            ENDIF
         ENDIF

         IF( iom_use("sophteiv") .OR. iom_use("sopsteiv") ) THEN 
            z2d(1,:) = htr_eiv(:,1) * rc_pwatt        !  (conversion in PW) 
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sophteiv'
            CALL iom_put( TRIM(cl1), z2d )
            z2d(1,:) = str_eiv(:,1) * rc_ggram        !  (conversion in Gg)
            DO ji = 1, jpi
               z2d(ji,:) = z2d(1,:)
            ENDDO
            cl1 = 'sopsteiv'
            CALL iom_put( TRIM(cl1), z2d )
            IF( ln_subbas ) THEN
               DO jn=2,nptr
                  z2d(1,:) = htr_eiv(:,jn) * rc_pwatt        !  (conversion in PW)
                  DO ji = 1, jpi
                     z2d(ji,:) = z2d(1,:)
                  ENDDO
                  cl1 = TRIM('sophteiv_'//clsubb(jn))                 
                  CALL iom_put( cl1, z2d )
                  z2d(1,:) = str_eiv(:,jn) * rc_ggram        ! (conversion in Gg)
                  DO ji = 1, jpi
                     z2d(ji,:) = z2d(1,:)
                  ENDDO
                  cl1 = TRIM('sopsteiv_'//clsubb(jn)) 
                  CALL iom_put( cl1, z2d )              
               ENDDO
            ENDIF
         ENDIF
         !
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('dia_ptr')
      !
   END SUBROUTINE dia_ptr


   SUBROUTINE dia_ptr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr_init  ***
      !!                   
      !! ** Purpose :   Initialization, namelist read
      !!----------------------------------------------------------------------
      INTEGER ::  jn           ! local integers
      INTEGER ::  inum, ierr   ! local integers
      INTEGER ::  ios          ! Local integer output status for namelist read
      !!
      NAMELIST/namptr/ ln_diaptr, ln_subbas
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )              ! Namelist namptr in reference namelist : Poleward transport
      READ  ( numnam_ref, namptr, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namptr in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namptr in configuration namelist : Poleward transport
      READ  ( numnam_cfg, namptr, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namptr in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namptr )

      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namptr : set ptr parameters'
         WRITE(numout,*) '      Poleward heat & salt transport (T) or not (F)      ln_diaptr  = ', ln_diaptr
         WRITE(numout,*) '      Global (F) or glo/Atl/Pac/Ind/Indo-Pac basins      ln_subbas  = ', ln_subbas
      ENDIF

      IF( ln_diaptr ) THEN  
         !
         IF( ln_subbas ) THEN 
            nptr = 5            ! Global, Atlantic, Pacific, Indian, Indo-Pacific
            ALLOCATE( clsubb(nptr) )
            clsubb(1) = 'glo' ;  clsubb(2) = 'atl'  ;  clsubb(3) = 'pac'  ;  clsubb(4) = 'ind'  ;  clsubb(5) = 'ipc'
         ELSE               
            nptr = 1       ! Global only
            ALLOCATE( clsubb(nptr) )
            clsubb(1) = 'glo' 
         ENDIF

         !                                      ! allocate dia_ptr arrays
         IF( dia_ptr_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_ptr_init : unable to allocate arrays' )

         rc_pwatt = rc_pwatt * rau0_rcp          ! conversion from K.s-1 to PetaWatt

         IF( lk_mpp )   CALL mpp_ini_znl( numout )     ! Define MPI communicator for zonal sum

         IF( ln_subbas ) THEN                ! load sub-basin mask
            CALL iom_open( 'subbasins', inum,  ldstop = .FALSE.  )
            CALL iom_get( inum, jpdom_data, 'atlmsk', btmsk(:,:,2) )   ! Atlantic basin
            CALL iom_get( inum, jpdom_data, 'pacmsk', btmsk(:,:,3) )   ! Pacific  basin
            CALL iom_get( inum, jpdom_data, 'indmsk', btmsk(:,:,4) )   ! Indian   basin
            CALL iom_close( inum )
            btmsk(:,:,5) = MAX ( btmsk(:,:,3), btmsk(:,:,4) )          ! Indo-Pacific basin
            WHERE( gphit(:,:) < -30._wp)   ;   btm30(:,:) = 0._wp      ! mask out Southern Ocean
            ELSE WHERE                     ;   btm30(:,:) = ssmask(:,:)
            END WHERE
         ENDIF
   
         btmsk(:,:,1) = tmask_i(:,:)                                   ! global ocean
      
         DO jn = 1, nptr
            btmsk(:,:,jn) = btmsk(:,:,jn) * tmask_i(:,:)               ! interior domain only
         END DO

         ! Initialise arrays to zero because diatpr is called before they are first calculated
         ! Note that this means diagnostics will not be exactly correct when model run is restarted.
         htr_adv(:,:) = 0._wp  ;  str_adv(:,:) =  0._wp 
         htr_ldf(:,:) = 0._wp  ;  str_ldf(:,:) =  0._wp 
         htr_eiv(:,:) = 0._wp  ;  str_eiv(:,:) =  0._wp 
         htr_ove(:,:) = 0._wp  ;   str_ove(:,:) =  0._wp
         htr_btr(:,:) = 0._wp  ;   str_btr(:,:) =  0._wp
         !
      ENDIF 
      ! 
   END SUBROUTINE dia_ptr_init


   SUBROUTINE dia_ptr_hst( ktra, cptr, pva ) 
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ptr_hst ***
      !!----------------------------------------------------------------------
      !! Wrapper for heat and salt transport calculations to calculate them for each basin
      !! Called from all advection and/or diffusion routines
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in )  :: ktra  ! tracer index
      CHARACTER(len=3)                , INTENT(in)   :: cptr  ! transport type  'adv'/'ldf'/'eiv'
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)   :: pva   ! 3D input array of advection/diffusion
      INTEGER                                        :: jn    !

      IF( cptr == 'adv' ) THEN
         IF( ktra == jp_tem )  htr_adv(:,1) = ptr_sj( pva(:,:,:) )
         IF( ktra == jp_sal )  str_adv(:,1) = ptr_sj( pva(:,:,:) )
      ENDIF
      IF( cptr == 'ldf' ) THEN
         IF( ktra == jp_tem )  htr_ldf(:,1) = ptr_sj( pva(:,:,:) )
         IF( ktra == jp_sal )  str_ldf(:,1) = ptr_sj( pva(:,:,:) )
      ENDIF
      IF( cptr == 'eiv' ) THEN
         IF( ktra == jp_tem )  htr_eiv(:,1) = ptr_sj( pva(:,:,:) )
         IF( ktra == jp_sal )  str_eiv(:,1) = ptr_sj( pva(:,:,:) )
      ENDIF
      !
      IF( ln_subbas ) THEN
         !
         IF( cptr == 'adv' ) THEN
             IF( ktra == jp_tem ) THEN 
                DO jn = 2, nptr
                   htr_adv(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                END DO
             ENDIF
             IF( ktra == jp_sal ) THEN 
                DO jn = 2, nptr
                   str_adv(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                END DO
             ENDIF
         ENDIF
         IF( cptr == 'ldf' ) THEN
             IF( ktra == jp_tem ) THEN 
                DO jn = 2, nptr
                    htr_ldf(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                 END DO
             ENDIF
             IF( ktra == jp_sal ) THEN 
                DO jn = 2, nptr
                   str_ldf(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                END DO
             ENDIF
         ENDIF
         IF( cptr == 'eiv' ) THEN
             IF( ktra == jp_tem ) THEN 
                DO jn = 2, nptr
                    htr_eiv(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                 END DO
             ENDIF
             IF( ktra == jp_sal ) THEN 
                DO jn = 2, nptr
                   str_eiv(:,jn) = ptr_sj( pva(:,:,:), btmsk(:,:,jn) )
                END DO
             ENDIF
         ENDIF
         !
      ENDIF
   END SUBROUTINE dia_ptr_hst


   FUNCTION dia_ptr_alloc()
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_ptr_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER               ::   dia_ptr_alloc   ! return value
      INTEGER, DIMENSION(3) ::   ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( btmsk(jpi,jpj,nptr) ,              &
         &      htr_adv(jpj,nptr) , str_adv(jpj,nptr) ,   &
         &      htr_eiv(jpj,nptr) , str_eiv(jpj,nptr) ,   &
         &      htr_ove(jpj,nptr) , str_ove(jpj,nptr) ,   &
         &      htr_btr(jpj,nptr) , str_btr(jpj,nptr) ,   &
         &      htr_ldf(jpj,nptr) , str_ldf(jpj,nptr) , STAT=ierr(1)  )
         !
      ALLOCATE( p_fval1d(jpj), p_fval2d(jpj,jpk), Stat=ierr(2))
      !
      ALLOCATE( btm30(jpi,jpj), STAT=ierr(3)  )

         !
      dia_ptr_alloc = MAXVAL( ierr )
      IF(lk_mpp)   CALL mpp_sum( dia_ptr_alloc )
      !
   END FUNCTION dia_ptr_alloc


   FUNCTION ptr_sj_3d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_3d  ***
      !!
      !! ** Purpose :   i-k sum computation of a j-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!              pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk)       ::   pva   ! mask flux array at V-point
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj), OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji, jj, jk   ! dummy loop arguments
      INTEGER                  ::   ijpj         ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval  ! function value
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      IF( PRESENT( pmsk ) ) THEN 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! Vector opt.
                  p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj) * pmsk(ji,jj)
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! Vector opt.
                  p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj) 
               END DO
            END DO
         END DO
      ENDIF
#if defined key_mpp_mpi
      IF(lk_mpp)   CALL mpp_sum( p_fval, ijpj, ncomm_znl)
#endif
      !
   END FUNCTION ptr_sj_3d


   FUNCTION ptr_sj_2d( pva, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sj_2d  ***
      !!
      !! ** Purpose :   "zonal" and vertical sum computation of a i-flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!----------------------------------------------------------------------
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)           ::   pva   ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj), OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !
      INTEGER                  ::   ji,jj       ! dummy loop arguments
      INTEGER                  ::   ijpj        ! ???
      REAL(wp), POINTER, DIMENSION(:) :: p_fval ! function value
      !!--------------------------------------------------------------------
      ! 
      p_fval => p_fval1d

      ijpj = jpj
      p_fval(:) = 0._wp
      IF( PRESENT( pmsk ) ) THEN 
         DO jj = 2, jpjm1
            DO ji = nldi, nlei   ! No vector optimisation here. Better use a mask ?
               p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj) * pmsk(ji,jj)
            END DO
         END DO
      ELSE
         DO jj = 2, jpjm1
            DO ji = nldi, nlei   ! No vector optimisation here. Better use a mask ?
               p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
      ENDIF
#if defined key_mpp_mpi
      CALL mpp_sum( p_fval, ijpj, ncomm_znl )
#endif
      ! 
   END FUNCTION ptr_sj_2d


   FUNCTION ptr_sjk( pta, pmsk )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_sjk  ***
      !!
      !! ** Purpose :   i-sum computation of an array
      !!
      !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
      !!
      !! ** Action  : - p_fval: i-mean poleward flux of pva
      !!----------------------------------------------------------------------
      !!
      IMPLICIT none
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk)           ::   pta    ! mask flux array at V-point
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj)    , OPTIONAL ::   pmsk   ! Optional 2D basin mask
      !!
      INTEGER                           :: ji, jj, jk ! dummy loop arguments
      REAL(wp), POINTER, DIMENSION(:,:) :: p_fval     ! return function value
#if defined key_mpp_mpi
      INTEGER, DIMENSION(1) ::   ish
      INTEGER, DIMENSION(2) ::   ish2
      INTEGER               ::   ijpjjpk
      REAL(wp), DIMENSION(jpj*jpk) ::   zwork    ! mask flux array at V-point
#endif
      !!--------------------------------------------------------------------
      !
      p_fval => p_fval2d

      p_fval(:,:) = 0._wp
      !
      IF( PRESENT( pmsk ) ) THEN 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
!!gm here, use of tmask_i  ==> no need of loop over nldi, nlei....
               DO ji =  nldi, nlei   ! No vector optimisation here. Better use a mask ?
                  p_fval(jj,jk) = p_fval(jj,jk) + pta(ji,jj,jk) * pmsk(ji,jj)
               END DO
            END DO
         END DO
      ELSE 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji =  nldi, nlei   ! No vector optimisation here. Better use a mask ?
                  p_fval(jj,jk) = p_fval(jj,jk) + pta(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
      END IF
      !
#if defined key_mpp_mpi
      ijpjjpk = jpj*jpk
      ish(1) = ijpjjpk  ;   ish2(1) = jpj   ;   ish2(2) = jpk
      zwork(1:ijpjjpk) = RESHAPE( p_fval, ish )
      CALL mpp_sum( zwork, ijpjjpk, ncomm_znl )
      p_fval(:,:) = RESHAPE( zwork, ish2 )
#endif
      !
   END FUNCTION ptr_sjk


   !!======================================================================
END MODULE diaptr
