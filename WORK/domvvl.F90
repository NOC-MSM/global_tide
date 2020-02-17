MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean : 
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!            3.3  !  2011-10  (M. Leclair) totally rewrote domvvl: vvl option includes z_star and z_tilde coordinates
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!            4.0  !  2018-09  (J. Chanut) improve z_tilde robustness
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_vvl_init     : define initial vertical scale factors, depths and column thickness
   !!   dom_vvl_sf_nxt   : Compute next vertical scale factors
   !!   dom_vvl_sf_swp   : Swap vertical scale factors and update the vertical grid
   !!   dom_vvl_interpol : Interpolate vertical scale factors from one grid point to another
   !!   dom_vvl_rst      : read/write restart file
   !!   dom_vvl_ctl      : Check the vvl options
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE phycst          ! physical constant
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! ocean surface boundary condition
   USE wet_dry         ! wetting and drying
   USE usrdef_istate   ! user defined initial state (wad only)
   USE restart         ! ocean restart
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! Timing
   USE bdy_oce         ! ocean open boundary conditions
   USE sbcrnf          ! river runoff 
   USE dynspg_ts, ONLY: un_adv, vn_adv

   IMPLICIT NONE
   PRIVATE

   PUBLIC  dom_vvl_init       ! called by domain.F90
   PUBLIC  dom_vvl_sf_nxt     ! called by step.F90
   PUBLIC  dom_vvl_sf_swp     ! called by step.F90
   PUBLIC  dom_vvl_interpol   ! called by dynnxt.F90

   !                                                      !!* Namelist nam_vvl
   LOGICAL , PUBLIC :: ln_vvl_zstar           = .FALSE.    ! zstar  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde          = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_layer           = .FALSE.    ! level  vertical coordinate
   LOGICAL          :: ln_vvl_ztilde_as_zstar = .FALSE.    ! ztilde vertical coordinate
   LOGICAL          :: ln_vvl_zstar_at_eqtor  = .FALSE.    ! ztilde vertical coordinate
   LOGICAL          :: ln_vvl_zstar_on_shelf  = .FALSE.    ! revert to zstar on shelves
   LOGICAL          :: ln_vvl_adv_fct         = .FALSE.    ! Centred thickness advection
   LOGICAL          :: ln_vvl_adv_cn2         = .TRUE.     ! FCT thickness advection
   LOGICAL          :: ln_vvl_dbg             = .FALSE.    ! debug control prints
   LOGICAL          :: ln_vvl_ramp            = .FALSE.    ! Ramp on interfaces displacement
   LOGICAL          :: ln_vvl_lap             = .FALSE.    ! Laplacian thickness diffusion
   LOGICAL          :: ln_vvl_blp             = .FALSE.    ! Bilaplacian thickness diffusion 
   LOGICAL          :: ln_vvl_regrid          = .FALSE.    ! ensure layer separation 
   LOGICAL          :: ll_shorizd             = .FALSE.    ! Use "shelf horizon depths" 
   LOGICAL          :: ln_vvl_kepe            = .FALSE.    ! kinetic/potential energy transfer
   !                                                       ! conservation: not used yet
   INTEGER          :: nn_filt_order=1
   REAL(wp)         :: rn_ahe3_lap               ! thickness diffusion coefficient (Laplacian)
   REAL(wp)         :: rn_ahe3_blp               ! thickness diffusion coefficient (Bilaplacian)
   REAL(wp)         :: rn_rst_e3t                ! ztilde to zstar restoration timescale [days]
   REAL(wp)         :: rn_lf_cutoff              ! cutoff frequency for low-pass filter  [days]
   REAL(wp)         :: rn_day_ramp               ! Duration of linear ramp  [days]
   REAL(wp)         :: hsmall=0.01_wp            ! small thickness [m]

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: un_td, vn_td                ! thickness diffusion transport
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: un_lf, vn_lf, hdivn_lf    ! low frequency fluxes and divergence
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_b, tilde_e3t_n    ! baroclinic scale factors
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_a                 ! baroclinic scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tildemask                   ! mask tilde tendency
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_e3t                 ! restoring period for scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_hdv                 ! restoring period for low freq. divergence
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: hsm, dsm                    ! 
   INTEGER         , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: i_int_bot

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domvvl.F90 10167 2018-10-03 13:52:45Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dom_vvl_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION dom_vvl_alloc  ***
      !!----------------------------------------------------------------------
      IF( ln_vvl_zstar )   dom_vvl_alloc = 0
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         ALLOCATE( tilde_e3t_b(jpi,jpj,jpk)  , tilde_e3t_n(jpi,jpj,jpk) , tilde_e3t_a(jpi,jpj,jpk) ,   &
            &      un_td  (jpi,jpj,jpk)      , vn_td  (jpi,jpj,jpk)     ,                              &
            &      tildemask(jpi,jpj) , hsm(jpi,jpj) , dsm(jpi,jpj) , i_int_bot(jpi,jpj), STAT = dom_vvl_alloc )
         IF( lk_mpp             )   CALL mpp_sum ( dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_warn('dom_vvl_alloc: failed to allocate arrays')
         un_td = 0._wp
         vn_td = 0._wp
      ENDIF
      IF( ln_vvl_ztilde ) THEN
         ALLOCATE( frq_rst_e3t(jpi,jpj) , frq_rst_hdv(jpi,jpj) , hdivn_lf(jpi,jpj,jpk,nn_filt_order),  &
            &      un_lf(jpi,jpj,jpk,nn_filt_order), vn_lf(jpi,jpj,jpk,nn_filt_order), STAT= dom_vvl_alloc )
         IF( lk_mpp             )   CALL mpp_sum ( dom_vvl_alloc )
         IF( dom_vvl_alloc /= 0 )   CALL ctl_warn('dom_vvl_alloc: failed to allocate arrays')
      ENDIF
      !
   END FUNCTION dom_vvl_alloc


   SUBROUTINE dom_vvl_init
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_init  ***
      !!                   
      !! ** Purpose :  Initialization of all scale factors, depths
      !!               and water column heights
      !!
      !! ** Method  :  - use restart file and/or initialize
      !!               - interpolate scale factors
      !!
      !! ** Action  : - e3t_(n/b) and tilde_e3t_(n/b)
      !!              - Regrid: e3(u/v)_n
      !!                        e3(u/v)_b       
      !!                        e3w_n           
      !!                        e3(u/v)w_b      
      !!                        e3(u/v)w_n      
      !!                        gdept_n, gdepw_n and gde3w_n
      !!              - h(t/u/v)_0
      !!              - frq_rst_e3t and frq_rst_hdv
      !!
      !! Reference  : Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk
      INTEGER ::   ii0, ii1, ij0, ij1
      REAL(wp)::   zcoef, zwgt, ztmp, zhmin, zhmax
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dom_vvl_init : Variable volume activated'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      !
      CALL dom_vvl_ctl     ! choose vertical coordinate (z_star, z_tilde or layer)
      !
      !                    ! Allocate module arrays
      IF( dom_vvl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dom_vvl_init : unable to allocate arrays' )
      !
      !                    ! Read or initialize e3t_(b/n), tilde_e3t_(b/n) and hdivn_lf
      CALL dom_vvl_rst( nit000, 'READ' )
      e3t_a(:,:,jpk) = e3t_0(:,:,jpk)  ! last level always inside the sea floor set one for all
      !
      !                    !== Set of all other vertical scale factors  ==!  (now and before)
      !                                ! Horizontal interpolation of e3t
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3u_b(:,:,:), 'U' )    ! from T to U
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3u_n(:,:,:), 'U' )
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3v_b(:,:,:), 'V' )    ! from T to V 
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3v_n(:,:,:), 'V' )
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3f_n(:,:,:), 'F' )    ! from U to F
      !                                ! Vertical interpolation of e3t,u,v 
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3w_n (:,:,:), 'W'  )  ! from T to W
      CALL dom_vvl_interpol( e3t_b(:,:,:), e3w_b (:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3uw_n(:,:,:), 'UW' )  ! from U to UW
      CALL dom_vvl_interpol( e3u_b(:,:,:), e3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_n(:,:,:), e3vw_n(:,:,:), 'VW' )  ! from V to UW
      CALL dom_vvl_interpol( e3v_b(:,:,:), e3vw_b(:,:,:), 'VW' )

      ! We need to define e3[tuv]_a for AGRIF initialisation (should not be a problem for the restartability...)
      e3t_a(:,:,:) = e3t_n(:,:,:)
      e3u_a(:,:,:) = e3u_n(:,:,:)
      e3v_a(:,:,:) = e3v_n(:,:,:)
      !
      !                    !==  depth of t and w-point  ==!   (set the isf depth as it is in the initial timestep)
      gdept_n(:,:,1) = 0.5_wp * e3w_n(:,:,1)       ! reference to the ocean surface (used for MLD and light penetration)
      gdepw_n(:,:,1) = 0.0_wp
      gde3w_n(:,:,1) = gdept_n(:,:,1) - sshn(:,:)  ! reference to a common level z=0 for hpg
      gdept_b(:,:,1) = 0.5_wp * e3w_b(:,:,1)
      gdepw_b(:,:,1) = 0.0_wp
      DO jk = 2, jpk                               ! vertical sum
         DO jj = 1,jpj
            DO ji = 1,jpi
               !    zcoef = tmask - wmask    ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
               !                             ! 1 everywhere from mbkt to mikt + 1 or 1 (if no isf)
               !                             ! 0.5 where jk = mikt     
!!gm ???????   BUG ?  gdept_n as well as gde3w_n  does not include the thickness of ISF ??
               zcoef = ( tmask(ji,jj,jk) - wmask(ji,jj,jk) )
               gdepw_n(ji,jj,jk) = gdepw_n(ji,jj,jk-1) + e3t_n(ji,jj,jk-1)
               gdept_n(ji,jj,jk) =      zcoef  * ( gdepw_n(ji,jj,jk  ) + 0.5 * e3w_n(ji,jj,jk))  &
                  &                + (1-zcoef) * ( gdept_n(ji,jj,jk-1) +       e3w_n(ji,jj,jk)) 
               gde3w_n(ji,jj,jk) = gdept_n(ji,jj,jk) - sshn(ji,jj)
               gdepw_b(ji,jj,jk) = gdepw_b(ji,jj,jk-1) + e3t_b(ji,jj,jk-1)
               gdept_b(ji,jj,jk) =      zcoef  * ( gdepw_b(ji,jj,jk  ) + 0.5 * e3w_b(ji,jj,jk))  &
                  &                + (1-zcoef) * ( gdept_b(ji,jj,jk-1) +       e3w_b(ji,jj,jk)) 
            END DO
         END DO
      END DO
      !
      !                    !==  thickness of the water column  !!   (ocean portion only)
      ht_n(:,:) = e3t_n(:,:,1) * tmask(:,:,1)   !!gm  BUG  :  this should be 1/2 * e3w(k=1) ....
      hu_b(:,:) = e3u_b(:,:,1) * umask(:,:,1)
      hu_n(:,:) = e3u_n(:,:,1) * umask(:,:,1)
      hv_b(:,:) = e3v_b(:,:,1) * vmask(:,:,1)
      hv_n(:,:) = e3v_n(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         ht_n(:,:) = ht_n(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
         hu_b(:,:) = hu_b(:,:) + e3u_b(:,:,jk) * umask(:,:,jk)
         hu_n(:,:) = hu_n(:,:) + e3u_n(:,:,jk) * umask(:,:,jk)
         hv_b(:,:) = hv_b(:,:) + e3v_b(:,:,jk) * vmask(:,:,jk)
         hv_n(:,:) = hv_n(:,:) + e3v_n(:,:,jk) * vmask(:,:,jk)
      END DO
      !
      !                    !==  inverse of water column thickness   ==!   (u- and v- points)
      r1_hu_b(:,:) = ssumask(:,:) / ( hu_b(:,:) + 1._wp - ssumask(:,:) )    ! _i mask due to ISF
      r1_hu_n(:,:) = ssumask(:,:) / ( hu_n(:,:) + 1._wp - ssumask(:,:) )
      r1_hv_b(:,:) = ssvmask(:,:) / ( hv_b(:,:) + 1._wp - ssvmask(:,:) )
      r1_hv_n(:,:) = ssvmask(:,:) / ( hv_n(:,:) + 1._wp - ssvmask(:,:) )

      !                    !==   z_tilde coordinate case  ==!   (Restoring frequencies)
      tildemask(:,:) = 1._wp

      IF( ln_vvl_ztilde ) THEN
!!gm : idea: add here a READ in a file of custumized restoring frequency
         ! Values in days provided via the namelist; use rsmall to avoid possible division by zero errors with faulty settings
         frq_rst_e3t(:,:) = 2.0_wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.0_wp )
         frq_rst_hdv(:,:) = 2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.0_wp )
         !
         IF( ln_vvl_ztilde_as_zstar ) THEN
            ! Ignore namelist settings and use these next two to emulate z-star using z-tilde
            frq_rst_e3t(:,:) = 0.0_wp 
            frq_rst_hdv(:,:) = 1.0_wp / rdt
            rn_lf_cutoff     = 2.0_wp * rpi * rdt / 86400._wp
            tildemask(:,:) = 0._wp
         ENDIF
        
         IF ( ln_vvl_zstar_at_eqtor ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
!!gm  case |gphi| >= 6 degrees is useless   initialized just above by default  
                  IF( ABS(gphit(ji,jj)) >= 6.) THEN
                     ! values outside the equatorial band and transition zone (ztilde)
                     frq_rst_e3t(ji,jj) =  2.0_wp * rpi / ( MAX( rn_rst_e3t  , rsmall ) * 86400.e0_wp )
!                     frq_rst_hdv(ji,jj) =  2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.e0_wp )
              
                  ELSEIF( ABS(gphit(ji,jj)) <= 2.5) THEN
                     ! values inside the equatorial band (ztilde as zstar)
                     frq_rst_e3t(ji,jj) =  0.0_wp
!                     frq_rst_hdv(ji,jj) =  1.0_wp / rdt
                     tildemask(ji,jj) = 0._wp
                  ELSE
                     ! values in the transition band (linearly vary from ztilde to ztilde as zstar values)
                     frq_rst_e3t(ji,jj) = 0.0_wp + (frq_rst_e3t(ji,jj)-0.0_wp)*0.5_wp   &
                        &            * (  1.0_wp - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                          * 180._wp / 3.5_wp ) )
!                     frq_rst_hdv(ji,jj) = (1.0_wp / rdt)                                &
!                        &            + (  frq_rst_hdv(ji,jj)-(1.e0_wp / rdt) )*0.5_wp   &
!                        &            * (  1._wp  - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
!                        &                                          * 180._wp / 3.5_wp ) )
                     tildemask(ji,jj) = 0.5_wp * (  1._wp  - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                                 * 180._wp / 3.5_wp ) )
                  ENDIF
               END DO
            END DO
         ENDIF
         !
         IF ( ln_vvl_zstar_on_shelf ) THEN
            zhmin = 50._wp
            zhmax = 100._wp
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zwgt = 1._wp
                  IF(( ht_0(ji,jj)>zhmin).AND.(ht_0(ji,jj) <=zhmax)) THEN
                     zwgt = (ht_0(ji,jj)-zhmin)/(zhmax-zhmin)
                  ELSEIF ( ht_0(ji,jj)<=zhmin) THEN
                     zwgt = 0._wp
                  ENDIF
                  frq_rst_e3t(ji,jj) = MIN(frq_rst_e3t(ji,jj), frq_rst_e3t(ji,jj)*zwgt)
                  tildemask(ji,jj)   = MIN(tildemask(ji,jj), zwgt)
               END DO
            END DO
         ENDIF
         !
         ztmp = MAXVAL( frq_rst_hdv(:,:) )
         IF( lk_mpp )   CALL mpp_max( ztmp )                 ! max over the global domain
         !
         IF ( (ztmp*rdt) > 1._wp) CALL ctl_stop( 'dom_vvl_init: rn_lf_cuttoff is too small' )
         !
      ENDIF

      IF( ln_vvl_layer ) THEN
         IF ( ln_vvl_zstar_on_shelf ) THEN
            zhmin = 50._wp
            zhmax = 100._wp
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zwgt = 1._wp
                  IF(( ht_0(ji,jj)>zhmin).AND.(ht_0(ji,jj) <=zhmax)) THEN
                     zwgt = (ht_0(ji,jj)-zhmin)/(zhmax-zhmin)
                  ELSEIF ( ht_0(ji,jj)<=zhmin) THEN
                     zwgt = 0._wp
                  ENDIF
                  tildemask(ji,jj)   = MIN(tildemask(ji,jj), zwgt)
               END DO
            END DO
         ENDIF
         IF ( ln_vvl_zstar_at_eqtor ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
!!gm  case |gphi| >= 6 degrees is useless   initialized just above by default  
                  IF( ABS(gphit(ji,jj)) >= 6.) THEN
                     ! values outside the equatorial band and transition zone (ztilde)
              
                  ELSEIF( ABS(gphit(ji,jj)) <= 2.5) THEN
                     ! values inside the equatorial band (ztilde as zstar)
                     tildemask(ji,jj) = 0._wp
                  ELSE
                     tildemask(ji,jj) = 0.5_wp * (  1._wp  - COS( rad*(ABS(gphit(ji,jj))-2.5_wp)  &
                        &                                                 * 180._wp / 3.5_wp ) )
                  ENDIF
               END DO
            END DO
         ENDIF
      ENDIF
      !
      IF(lwxios) THEN
! define variables in restart file when writing with XIOS
         CALL iom_set_rstw_var_active('e3t_b')
         CALL iom_set_rstw_var_active('e3t_n')
         !                                           ! ----------------------- !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN  ! z_tilde and layer cases !
            !                                        ! ----------------------- !
            CALL iom_set_rstw_var_active('tilde_e3t_b')
            CALL iom_set_rstw_var_active('tilde_e3t_n')
         END IF
         !                                           ! -------------!    
         IF( ln_vvl_ztilde ) THEN                    ! z_tilde case !
            !                                        ! ------------ !
            CALL iom_set_rstw_var_active('un_lf')
            CALL iom_set_rstw_var_active('vn_lf')
            CALL iom_set_rstw_var_active('hdivn_lf')
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE dom_vvl_init


   SUBROUTINE dom_vvl_sf_nxt( kt, kcall ) 
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_nxt  ***
      !!                   
      !! ** Purpose :  - compute the after scale factors used in tra_zdf, dynnxt,
      !!                 tranxt and dynspg routines
      !!
      !! ** Method  :  - z_star case:  Repartition of ssh INCREMENT proportionnaly to the level thickness.
      !!               - z_tilde_case: after scale factor increment = 
      !!                                    high frequency part of horizontal divergence
      !!                                  + retsoring towards the background grid
      !!                                  + thickness difusion
      !!                               Then repartition of ssh INCREMENT proportionnaly
      !!                               to the "baroclinic" level thickness.
      !!
      !! ** Action  :  - hdivn_lf   : restoring towards full baroclinic divergence in z_tilde case
      !!               - tilde_e3t_a: after increment of vertical scale factor 
      !!                              in z_tilde case
      !!               - e3(t/u/v)_a
      !!
      !! Reference  : Leclair, M., and Madec, G. 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )           ::   kt      ! time step
      INTEGER, INTENT( in ), OPTIONAL ::   kcall   ! optional argument indicating call sequence
      !
      LOGICAL                ::   ll_do_bclinic         ! local logical
      INTEGER                ::   ji, jj, jk, jo        ! dummy loop indices
      INTEGER                ::   ib, ib_bdy, ip, jp    !   "     "     "
      INTEGER , DIMENSION(3) ::   ijk_max, ijk_min      ! temporary integers
      INTEGER                ::   ncall
      REAL(wp)               ::   z2dt  , z_tmin, z_tmax! local scalars                
      REAL(wp)               ::   zalpha, zwgt          ! temporary scalars
      REAL(wp)               ::   zdu, zdv, zramp, zmet
      REAL(wp)               ::   ztra, zbtr, ztout, ztin, zfac, zmsku, zmskv
      REAL(wp), DIMENSION(jpi,jpj)     ::   zht, z_scale, zwu, zwv, zhdiv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze3t, ztu, ztv
      !!----------------------------------------------------------------------
      !
      IF( ln_linssh )   RETURN      ! No calculation in linear free surface
      !
      IF( ln_timing )   CALL timing_start('dom_vvl_sf_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_nxt : compute after scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      zmet = 1._wp

      ll_do_bclinic = .TRUE.
      ncall = 1

      IF( PRESENT(kcall) ) THEN
         ! comment line below if tilda coordinate has to be computed at each call
         IF (kcall == 2 .AND. (ln_vvl_ztilde.OR.ln_vvl_layer) ) ll_do_bclinic = .FALSE.
         ncall = kcall  
      ENDIF

      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dt =  rdt
      ELSE
         z2dt = 2.0_wp * rdt
      ENDIF

      ! ******************************* !
      ! After scale factors at t-points !
      ! ******************************* !
      !                                                ! --------------------------------------------- !
      !                                                ! z_star coordinate and barotropic z-tilde part !
      !                                                ! --------------------------------------------- !
      !
      z_scale(:,:) = ( ssha(:,:) - sshb(:,:) ) * ssmask(:,:) / ( ht_0(:,:) + sshn(:,:) + 1. - ssmask(:,:) )
      DO jk = 1, jpkm1
         ! formally this is the same as e3t_a = e3t_0*(1+ssha/ht_0)
         e3t_a(:,:,jk) = e3t_b(:,:,jk) + e3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
      END DO
      !
      IF ((ln_vvl_ztilde.OR.ln_vvl_layer).AND.(zmet==1._wp)) THEN
         DO jk = 1, jpkm1
            e3t_a(:,:,jk) = e3t_a(:,:,jk) - tilde_e3t_n(:,:,jk) * z_scale(:,:) * tmask(:,:,jk)
         END DO  
      ENDIF

      IF( (ln_vvl_ztilde.OR.ln_vvl_layer) .AND. ll_do_bclinic ) THEN   ! z_tilde or layer coordinate !
         !                                                             ! ------baroclinic part------ !
         tilde_e3t_a(:,:,:) = 0.0_wp  ! tilde_e3t_a used to store tendency terms
         un_td(:,:,:) = 0.0_wp        ! Transport corrections
         vn_td(:,:,:) = 0.0_wp

         zhdiv(:,:) = 0.
         DO jk = 1, jpkm1
            zhdiv(:,:) = zhdiv(:,:) + e3t_n(:,:,jk) * hdivn(:,:,jk)
         END DO
         zhdiv(:,:) = zhdiv(:,:) / ( ht_0(:,:) + sshn(:,:) + 1._wp - tmask_i(:,:) ) 

         ze3t(:,:,:) = 0._wp
         IF( ln_rnf ) THEN
            CALL sbc_rnf_div( ze3t )          ! runoffs 
            DO jk=1,jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - ze3t(:,:,jk) * e3t_n(:,:,jk)
            END DO   
         ENDIF 

         ! Thickness advection:
         ! --------------------
         ! Set advection velocities and source term
         IF ( ln_vvl_ztilde ) THEN
            IF ( ncall==1 ) THEN
               zalpha  = rdt * 2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.0_wp )
               DO jk = 1, jpkm1
                  ztu(:,:,jk) = un(:,:,jk) * e3u_n(:,:,jk) * e2u(:,:)
                  ztv(:,:,jk) = vn(:,:,jk) * e3v_n(:,:,jk) * e1v(:,:)
                  ze3t(:,:,jk) = -tilde_e3t_a(:,:,jk) - (e3t_n(:,:,jk)-zmet*tilde_e3t_n(:,:,jk)) * zhdiv(:,:)
               END DO
               !
                  un_lf(:,:,:,nn_filt_order)  =    un_lf(:,:,:,nn_filt_order) * (1._wp - zalpha) + zalpha * ztu(:,:,:)
                  vn_lf(:,:,:,nn_filt_order)  =    vn_lf(:,:,:,nn_filt_order) * (1._wp - zalpha) + zalpha * ztv(:,:,:)
               hdivn_lf(:,:,:,nn_filt_order)  = hdivn_lf(:,:,:,nn_filt_order) * (1._wp - zalpha) + zalpha * ze3t(:,:,:)

               DO jo = nn_filt_order-1,1,-1
                     un_lf(:,:,:,jo)  =    un_lf(:,:,:,jo) * (1._wp - zalpha) + zalpha *    un_lf(:,:,:,jo+1)
                     vn_lf(:,:,:,jo)  =    vn_lf(:,:,:,jo) * (1._wp - zalpha) + zalpha *    vn_lf(:,:,:,jo+1)
                  hdivn_lf(:,:,:,jo)  = hdivn_lf(:,:,:,jo) * (1._wp - zalpha) + zalpha * hdivn_lf(:,:,:,jo+1)
               END DO
            ENDIF
            !
            DO jk = 1, jpkm1
               ztu(:,:,jk) = (un(:,:,jk)-un_lf(:,:,jk,1)/e3u_n(:,:,jk)*r1_e2u(:,:))*umask(:,:,jk)
               ztv(:,:,jk) = (vn(:,:,jk)-vn_lf(:,:,jk,1)/e3v_n(:,:,jk)*r1_e1v(:,:))*vmask(:,:,jk)
               tilde_e3t_a(:,:,jk) =  tilde_e3t_a(:,:,jk) + hdivn_lf(:,:,jk,1) - frq_rst_e3t(:,:) * tilde_e3t_b(:,:,jk)
            END DO
            !
         ELSEIF ( ln_vvl_layer ) THEN
            !
            DO jk = 1, jpkm1
               ztu(:,:,jk) = un(:,:,jk)
               ztv(:,:,jk) = vn(:,:,jk)
            END DO
            !
         ENDIF
         !
         ! Block fluxes through small layers:
!         DO jk=1,jpkm1
!            DO ji = 1, jpi
!               DO jj= 1, jpj
!                  zmsku = 0.5_wp * ( 1._wp + SIGN(1._wp, e3u_n(ji,jj,jk) - hsmall) )
!                  un_td(ji,jj,jk) = un_td(ji,jj,jk) - (1. - zmsku) * un(ji,jj,jk) * e3u_n(ji,jj,jk) * e2u(ji,jj)
!                  ztu(ji,jj,jk) = zmsku * ztu(ji,jj,jk)
!                  IF ( ln_vvl_ztilde ) un_lf(ji,jj,jk) = zmsku * un_lf(ji,jj,jk)
!                  !
!                  zmskv = 0.5_wp * ( 1._wp + SIGN(1._wp, e3v_n(ji,jj,jk) - hsmall) )
!                  vn_td(ji,jj,jk) = vn_td(ji,jj,jk) - (1. - zmskv) * vn(ji,jj,jk) * e3v_n(ji,jj,jk) * e1v(ji,jj)
!                  ztv(ji,jj,jk) = zmskv * ztv(ji,jj,jk)
!                  IF ( ln_vvl_ztilde ) vn_lf(ji,jj,jk) = zmskv * vn_lf(ji,jj,jk)
!               END DO
!            END DO
!         END DO      
         !
         ! Do advection
         IF     (ln_vvl_adv_fct) THEN
            CALL dom_vvl_adv_fct( kt, tilde_e3t_a, ztu, ztv )
            !
         ELSEIF (ln_vvl_adv_cn2) THEN
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     tilde_e3t_a(ji,jj,jk) =   tilde_e3t_a(ji,jj,jk) &
                     & -(  e2u(ji,jj)*e3u_n(ji,jj,jk) * ztu(ji,jj,jk) - e2u(ji-1,jj  )*e3u_n(ji-1,jj  ,jk) * ztu(ji-1,jj  ,jk)       &
                     &   + e1v(ji,jj)*e3v_n(ji,jj,jk) * ztv(ji,jj,jk) - e1v(ji  ,jj-1)*e3v_n(ji  ,jj-1,jk) * ztv(ji  ,jj-1,jk)  )    &
                     & / ( e1t(ji,jj) * e2t(ji,jj) )
                  END DO
               END DO
            END DO
         ENDIF
         ! 
         ! Thickness anomaly diffusion:
         ! ----------------------------
         ztu(:,:,:) = 0.0_wp
         ztv(:,:,:) = 0.0_wp

         ze3t(:,:,1) = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 2, jpk
                  ze3t(ji,jj,jk) =  ze3t(ji,jj,jk-1) + tilde_e3t_b(ji,jj,jk-1) * tmask(ji,jj,jk-1) 
               END DO
            END DO
         END DO

         IF ( ln_vvl_blp ) THEN  ! Bilaplacian
            DO jk = 1, jpkm1
               DO jj = 1, jpjm1                 ! First derivative (gradient)
                  DO ji = 1, fs_jpim1   ! vector opt.                  
                     ztu(ji,jj,jk) =  umask(ji,jj,jk) * e2_e1u(ji,jj) &
                                     &  * ( ze3t(ji,jj,jk) - ze3t(ji+1,jj  ,jk) )
                     ztv(ji,jj,jk) =  vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                                     &  * ( ze3t(ji,jj,jk) - ze3t(ji  ,jj+1,jk) )
                  END DO
               END DO

               DO jj = 2, jpjm1                 ! Second derivative (divergence) time the eddy diffusivity coefficient
                  DO ji = fs_2, fs_jpim1 ! vector opt.
                     zht(ji,jj) =  rn_ahe3_blp * r1_e1e2t(ji,jj) * (   ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                        &                                            + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)   )
                  END DO
               END DO

               ! Open boundary conditions:
               IF ( ln_bdy ) THEN
                  DO ib_bdy=1, nb_bdy
                     DO ib = 1, idx_bdy(ib_bdy)%nblenrim(1)
                        ji = idx_bdy(ib_bdy)%nbi(ib,1)
                        jj = idx_bdy(ib_bdy)%nbj(ib,1)
                        zht(ji,jj) = 0._wp
                     END DO 
                  END DO
               END IF

               CALL lbc_lnk( zht, 'T', 1. )     ! Lateral boundary conditions (unchanged sgn)

               DO jj = 1, jpjm1                 ! third derivative (gradient)
                  DO ji = 1, fs_jpim1   ! vector opt.
                     ztu(ji,jj,jk) = umask(ji,jj,jk) * e2_e1u(ji,jj) * ( zht(ji+1,jj  ) - zht(ji,jj) )
                     ztv(ji,jj,jk) = vmask(ji,jj,jk) * e1_e2v(ji,jj) * ( zht(ji  ,jj+1) - zht(ji,jj) )
                  END DO
               END DO
            END DO
         ENDIF

         IF ( ln_vvl_lap ) THEN    ! Laplacian
            DO jk = 1, jpkm1                    ! First derivative (gradient)
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.                  
                     zdu =  rn_ahe3_lap * umask(ji,jj,jk) * e2_e1u(ji,jj) &
                         &  * ( ze3t(ji,jj,jk) - ze3t(ji+1,jj  ,jk) )
                     zdv =  rn_ahe3_lap * vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                         &  * ( ze3t(ji,jj,jk) - ze3t(ji  ,jj+1,jk) )
                     ztu(ji,jj,jk) = ztu(ji,jj,jk) + zdu
                     ztv(ji,jj,jk) = ztv(ji,jj,jk) + zdv
                  END DO
               END DO
            END DO
         ENDIF 

         ! divergence of diffusive fluxes
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  un_td(ji,jj,jk) = un_td(ji,jj,jk) + ztu(ji,jj,jk+1) - ztu(ji,jj,jk  ) 
                  vn_td(ji,jj,jk) = vn_td(ji,jj,jk) + ztv(ji,jj,jk+1) - ztv(ji,jj,jk  ) 
                  tilde_e3t_a(ji,jj,jk) = tilde_e3t_a(ji,jj,jk) + (   ztu(ji-1,jj  ,jk+1) - ztu(ji,jj,jk+1)    &
                     &                                               +ztv(ji  ,jj-1,jk+1) - ztv(ji,jj,jk+1)    &
                     &                                               -ztu(ji-1,jj  ,jk  ) + ztu(ji,jj,jk  )    &
                     &                                               -ztv(ji  ,jj-1,jk  ) + ztv(ji,jj,jk  )    &
                     &                                            ) * r1_e1e2t(ji,jj)
               END DO
            END DO
         END DO

         CALL lbc_lnk_multi( un_td, 'U', -1., vn_td, 'V', -1. )     !* local domain boundaries
         !
         CALL dom_vvl_ups_cor( kt, tilde_e3t_a, un_td, vn_td )

!         IF ( ln_vvl_ztilde ) THEN
!            ztu(:,:,:) = -un_lf(:,:,:)
!            ztv(:,:,:) = -vn_lf(:,:,:)
!            CALL dom_vvl_ups_cor( kt, tilde_e3t_a, ztu, ztv )
!         ENDIF
         !
         ! Remove "external thickness" tendency:
         DO jk = 1, jpkm1
            tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) +   (e3t_n(:,:,jk)-zmet*tilde_e3t_n(:,:,jk)) * zhdiv(:,:) * tmask(:,:,jk)
         END DO 
                   
         ! Leapfrog time stepping
         ! ~~~~~~~~~~~~~~~~~~~~~~
         zramp = 1._wp
         IF ((.NOT.ln_rstart).AND.ln_vvl_ramp) zramp = MIN(MAX( ((kt-nit000)*rdt)/(rn_day_ramp*rday),0._wp),1._wp)

         DO jk=1,jpkm1
            tilde_e3t_a(:,:,jk) = tilde_e3t_b(:,:,jk) + z2dt * tmask(:,:,jk) * tilde_e3t_a(:,:,jk) &
                               & * tildemask(:,:) * zramp
         END DO

         ! Ensure layer separation:
         ! ~~~~~~~~~~~~~~~~~~~~~~~~
         IF ( ln_vvl_regrid ) CALL dom_vvl_regrid( kt )
  
         ! Boundary conditions:
         ! ~~~~~~~~~~~~~~~~~~~~
         IF ( ln_bdy ) THEN
            DO ib_bdy = 1, nb_bdy
               DO ib = 1, idx_bdy(ib_bdy)%nblenrim(1)
!!               DO ib = 1, idx_bdy(ib_bdy)%nblen(1)
                  ji = idx_bdy(ib_bdy)%nbi(ib,1)
                  jj = idx_bdy(ib_bdy)%nbj(ib,1)
                  zwgt = idx_bdy(ib_bdy)%nbw(ib,1)
                  ip = bdytmask(ji+1,jj  ) - bdytmask(ji-1,jj  )
                  jp = bdytmask(ji  ,jj+1) - bdytmask(ji  ,jj-1)
                  DO jk = 1, jpkm1  
                     tilde_e3t_a(ji,jj,jk) = 0.e0
!!                     tilde_e3t_a(ji,jj,jk) = tilde_e3t_a(ji,jj,jk) * (1._wp - zwgt)
!!                     tilde_e3t_a(ji,jj,jk) = tilde_e3t_a(ji+ip,jj+jp,jk) * tmask(ji+ip,jj+jp,jk)
                  END DO
               END DO 
            END DO
         ENDIF

         CALL lbc_lnk( tilde_e3t_a(:,:,:), 'T', 1. )        
        
      ENDIF

      IF( ln_vvl_ztilde.AND.( ncall==1))  THEN
         zalpha  = rdt * 2.0_wp * rpi / ( MAX( rn_lf_cutoff, rsmall ) * 86400.0_wp )
         !
         ! divergence of diffusive fluxes
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ze3t(ji,jj,jk) = (  un_td(ji,jj,jk) - un_td(ji-1,jj  ,jk) &
                                 &  + vn_td(ji,jj,jk) - vn_td(ji  ,jj-1,jk) &
                                 & ) / ( e1t(ji,jj) * e2t(ji,jj) )
               END DO
            END DO
         END DO
         CALL lbc_lnk( ze3t(:,:,:), 'T', 1. )

         DO jo = nn_filt_order,1,-1
            hdivn_lf(:,:,:,jo) =  hdivn_lf(:,:,:,jo)  + zalpha**(nn_filt_order-jo+1) * ze3t(:,:,:)
         END DO
      ENDIF

      IF( ln_vvl_ztilde .OR. ln_vvl_layer )  THEN   ! z_tilde or layer coordinate !
      !                                             ! ---baroclinic part--------- !

         IF ( (ncall==2).AND.(.NOT.ll_do_bclinic) ) THEN

            DO jk = 1, jpkm1
               ztu(:,:,jk) = (un_adv(:,:)*r1_hu_n(:,:) - un_b(:,:) ) * e3u_n(:,:,jk) * e2u(:,:) * umask(:,:,jk)
               ztv(:,:,jk) = (vn_adv(:,:)*r1_hv_n(:,:) - vn_b(:,:) ) * e3v_n(:,:,jk) * e1v(:,:) * vmask(:,:,jk)
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     ze3t(ji,jj,jk) = (  ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk) &
                                    &  + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk) &
                                    & ) / ( e1t(ji,jj) * e2t(ji,jj) )
                  END DO
               END DO
            END DO
            !
            zhdiv(:,:) = 0.
            DO jk = 1, jpkm1
               zhdiv(:,:) = zhdiv(:,:) + ze3t(:,:,jk) * tmask(:,:,jk)
            END DO
            zhdiv(:,:) = zhdiv(:,:) / ( ht_0(:,:) + sshn(:,:) + 1._wp - tmask_i(:,:) ) 
            !
            DO jk = 1, jpkm1
               tilde_e3t_a(:,:,jk) = tilde_e3t_a(:,:,jk) - z2dt * (ze3t(:,:,jk) & 
                                   & - zhdiv(:,:)*(e3t_n(:,:,jk)-zmet*tilde_e3t_n(:,:,jk))*tmask(:,:,jk))
            END DO
            CALL lbc_lnk( tilde_e3t_a(:,:,:), 'T', 1. )
         ENDIF
         !
         DO jk = 1, jpkm1
            e3t_a(:,:,jk) = e3t_a(:,:,jk) + tilde_e3t_a(:,:,jk) - tilde_e3t_b(:,:,jk)
         END DO
!         IF( ln_vvl_ztilde ) THEN !  Relax barotropic component:
!            DO jk = 1, jpkm1
!               e3t_a(:,:,jk) = e3t_a(:,:,jk)  &
!                             & - z2dt * frq_rst_e3t(:,:) * (e3t_b(:,:,jk) - tilde_e3t_b(:,:,jk)   & 
!                             & - e3t_0(:,:,jk) * (ht_0(:,:) + sshb(:,:))/ (ht_0(:,:)*tmask(:,:,1) + 1._wp - tmask(:,:,1)))                 
!            END DO
!         ENDIF
      ENDIF

      IF( ln_vvl_dbg.AND.(ln_vvl_ztilde .OR. ln_vvl_layer) ) THEN   ! - ML - test: control prints for debuging
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + tilde_e3t_a(:,:,jk) * tmask(:,:,jk)
         END DO
         IF( lwp ) WRITE(numout, *) 'kt = ', kt
         IF( lwp ) WRITE(numout, *) 'ncall = ', ncall
         IF( lwp ) WRITE(numout, *) 'll_do_bclinic', ll_do_bclinic 
         IF ( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
            z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( zht(:,:) ), mask = tmask(:,:,1) == 1.e0 )
            IF( lk_mpp ) CALL mpp_max( z_tmax )                             ! max over the global domain
            IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(SUM(tilde_e3t_a))) =', z_tmax
         END IF
         !
         z_tmin = MINVAL( e3t_n(:,:,:)  )
         IF( lk_mpp )   CALL mpp_min( z_tmin )                 ! min over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MINVAL(e3t_n) =', z_tmin
         IF ( z_tmin .LE. 0._wp ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_minloc(e3t_n(:,:,:), tmask, z_tmin, ijk_min(1), ijk_min(2), ijk_min(3) )
            ELSE
               ijk_min = MINLOC( e3t_n(:,:,:)  )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               ji = ijk_min(1) ; jj = ijk_min(2) ; jk = ijk_min(3)
               WRITE(numout, *) 'Negative scale factor, e3t_n =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_stop('dom_vvl_sf_nxt: Negative scale factor')
            ENDIF
         ENDIF         
         !
         z_tmin = MINVAL( e3u_n(:,:,:))
         IF( lk_mpp )   CALL mpp_min( z_tmin )                 ! min over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MINVAL(e3u_n) =', z_tmin
         IF ( z_tmin .LE. 0._wp ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_minloc(e3u_n(:,:,:), umask, z_tmin, ijk_min(1), ijk_min(2), ijk_min(3) )
            ELSE
               ijk_min = MINLOC( e3u_n(:,:,:)  )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               WRITE(numout, *) 'Negative scale factor, e3u_n =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_stop('dom_vvl_sf_nxt: Negative scale factor')
            ENDIF
         ENDIF 
         !
         z_tmin = MINVAL( e3v_n(:,:,:) )
         IF( lk_mpp )   CALL mpp_min( z_tmin )                 ! min over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MINVAL(e3v_n) =', z_tmin
         IF ( z_tmin .LE. 0._wp ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_minloc(e3v_n(:,:,:), vmask, z_tmin, ijk_min(1), ijk_min(2), ijk_min(3) )
            ELSE
               ijk_min = MINLOC( e3v_n(:,:,:), mask = vmask(:,:,:) == 1.e0   )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               WRITE(numout, *) 'Negative scale factor, e3v_n =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_stop('dom_vvl_sf_nxt: Negative scale factor')
            ENDIF
         ENDIF 
         !
         z_tmin = MINVAL( e3f_n(:,:,:))
         IF( lk_mpp )   CALL mpp_min( z_tmin )                 ! min over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MINVAL(e3f_n) =', z_tmin
         IF ( z_tmin .LE. 0._wp ) THEN
            IF( lk_mpp ) THEN
               CALL mpp_minloc(e3f_n(:,:,:), fmask, z_tmin, ijk_min(1), ijk_min(2), ijk_min(3) )
            ELSE
               ijk_min = MINLOC( e3f_n(:,:,:) )
               ijk_min(1) = ijk_min(1) + nimpp - 1
               ijk_min(2) = ijk_min(2) + njmpp - 1
            ENDIF
            IF (lwp) THEN
               WRITE(numout, *) 'Negative scale factor, e3f_n =', z_tmin
               WRITE(numout, *) 'at i, j, k=', ijk_min            
               CALL ctl_stop('dom_vvl_sf_nxt: Negative scale factor')
            ENDIF
         ENDIF
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
         END DO
         z_tmax = MAXVAL( tmask(:,:,1) * tmask_i(:,:) * ABS( ht_0(:,:) + sshn(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(ht_0+sshn  -SUM(e3t_n))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3u_n(:,:,jk) * umask(:,:,jk)
         END DO
         zwu(:,:) = 0._wp
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zwu(ji,jj) = 0.5_wp * umask(ji,jj,1) * r1_e1e2u(ji,jj)                             &
                        &          * ( e1e2t(ji,jj) * sshn(ji,jj) + e1e2t(ji+1,jj) * sshn(ji+1,jj) )
            END DO
         END DO
         CALL lbc_lnk( zwu(:,:), 'U', 1._wp )
         z_tmax = MAXVAL( ssumask(:,:) * ssumask(:,:) * ABS( hu_0(:,:) + zwu(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(hu_0+sshu_n-SUM(e3u_n))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            zht(:,:) = zht(:,:) + e3v_n(:,:,jk) * vmask(:,:,jk)
         END DO
         zwv(:,:) = 0._wp
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zwv(ji,jj) = 0.5_wp * vmask(ji,jj,1) * r1_e1e2v(ji,jj)                             &
                        &          * ( e1e2t(ji,jj) * sshn(ji,jj) + e1e2t(ji,jj+1) * sshn(ji,jj+1) )
            END DO
         END DO
         CALL lbc_lnk( zwv(:,:), 'V', 1._wp )
         z_tmax = MAXVAL( ssvmask(:,:) * ssvmask(:,:) * ABS( hv_0(:,:) + zwv(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(hv_0+sshv_n-SUM(e3v_n))) =', z_tmax
         !
         zht(:,:) = 0.0_wp
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               zht(:,jj) = zht(:,jj) + e3f_n(:,jj,jk) * umask(:,jj,jk)*umask(:,jj+1,jk)
            END DO
         END DO
         zwu(:,:) = 0._wp
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zwu(ji,jj) = 0.25_wp * umask(ji,jj,1) * umask(ji,jj+1,1) * r1_e1e2f(ji,jj)  &
                        &          * (  e1e2t(ji  ,jj) * sshn(ji  ,jj) + e1e2t(ji  ,jj+1) * sshn(ji  ,jj+1) &
                        &             + e1e2t(ji+1,jj) * sshn(ji+1,jj) + e1e2t(ji+1,jj+1) * sshn(ji+1,jj+1) )
            END DO
         END DO
         CALL lbc_lnk( zht(:,:), 'F', 1._wp )
         CALL lbc_lnk( zwu(:,:), 'F', 1._wp )
         z_tmax = MAXVAL( fmask(:,:,1) * fmask(:,:,1) * ABS( hf_0(:,:) + zwu(:,:) - zht(:,:) ) )
         IF( lk_mpp ) CALL mpp_max( z_tmax )                                ! max over the global domain
         IF( lwp    ) WRITE(numout, *) kt,' MAXVAL(abs(hf_0+sshf_n-SUM(e3f_n))) =', z_tmax
         !
      END IF

      ! *********************************** !
      ! After scale factors at u- v- points !
      ! *********************************** !

      CALL dom_vvl_interpol( e3t_a(:,:,:), e3u_a(:,:,:), 'U' )
      CALL dom_vvl_interpol( e3t_a(:,:,:), e3v_a(:,:,:), 'V' )

      ! *********************************** !
      ! After depths at u- v points         !
      ! *********************************** !

      hu_a(:,:) = e3u_a(:,:,1) * umask(:,:,1)
      hv_a(:,:) = e3v_a(:,:,1) * vmask(:,:,1)
      DO jk = 2, jpkm1
         hu_a(:,:) = hu_a(:,:) + e3u_a(:,:,jk) * umask(:,:,jk)
         hv_a(:,:) = hv_a(:,:) + e3v_a(:,:,jk) * vmask(:,:,jk)
      END DO
      !                                        ! Inverse of the local depth
!!gm BUG ?  don't understand the use of umask_i here .....
      r1_hu_a(:,:) = ssumask(:,:) / ( hu_a(:,:) + 1._wp - ssumask(:,:) )
      r1_hv_a(:,:) = ssvmask(:,:) / ( hv_a(:,:) + 1._wp - ssvmask(:,:) )
      !
      IF( ln_timing )   CALL timing_stop('dom_vvl_sf_nxt')
      !
   END SUBROUTINE dom_vvl_sf_nxt


   SUBROUTINE dom_vvl_sf_swp( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_sf_swp  ***
      !!                   
      !! ** Purpose :  compute time filter and swap of scale factors 
      !!               compute all depths and related variables for next time step
      !!               write outputs and restart file
      !!
      !! ** Method  :  - swap of e3t with trick for volume/tracer conservation
      !!               - reconstruct scale factor at other grid points (interpolate)
      !!               - recompute depths and water height fields
      !!
      !! ** Action  :  - e3t_(b/n), tilde_e3t_(b/n) and e3(u/v)_n ready for next time step
      !!               - Recompute:
      !!                    e3(u/v)_b       
      !!                    e3w_n           
      !!                    e3(u/v)w_b      
      !!                    e3(u/v)w_n      
      !!                    gdept_n, gdepw_n  and gde3w_n
      !!                    h(u/v) and h(u/v)r
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!              Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! time step
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zcoef        ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_linssh )   RETURN      ! No calculation in linear free surface
      !
      IF( ln_timing )   CALL timing_start('dom_vvl_sf_swp')
      !
      IF( kt == nit000 )   THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dom_vvl_sf_swp : - time filter and swap of scale factors'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   - interpolate scale factors and compute depths for next time step'
      ENDIF
      !
      ! Time filter and swap of scale factors
      ! =====================================
      ! - ML - e3(t/u/v)_b are allready computed in dynnxt.
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:)
         ELSE
            tilde_e3t_b(:,:,:) = tilde_e3t_n(:,:,:) & 
            &         + atfp * ( tilde_e3t_b(:,:,:) - 2.0_wp * tilde_e3t_n(:,:,:) + tilde_e3t_a(:,:,:) )
         ENDIF
         tilde_e3t_n(:,:,:) = tilde_e3t_a(:,:,:)
      ENDIF
      gdept_b(:,:,:) = gdept_n(:,:,:)
      gdepw_b(:,:,:) = gdepw_n(:,:,:)

      e3t_n(:,:,:) = e3t_a(:,:,:)
      e3u_n(:,:,:) = e3u_a(:,:,:)
      e3v_n(:,:,:) = e3v_a(:,:,:)

      ! Compute all missing vertical scale factor and depths
      ! ====================================================
      ! Horizontal scale factor interpolations
      ! --------------------------------------
      ! - ML - e3u_b and e3v_b are allready computed in dynnxt
      ! - JC - hu_b, hv_b, hur_b, hvr_b also
      
      CALL dom_vvl_interpol( e3t_n(:,:,:), e3f_n(:,:,:), 'F'  )
      
      ! Vertical scale factor interpolations
      CALL dom_vvl_interpol( e3t_n(:,:,:),  e3w_n(:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_n(:,:,:), e3uw_n(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_n(:,:,:), e3vw_n(:,:,:), 'VW' )
      CALL dom_vvl_interpol( e3t_b(:,:,:),  e3w_b(:,:,:), 'W'  )
      CALL dom_vvl_interpol( e3u_b(:,:,:), e3uw_b(:,:,:), 'UW' )
      CALL dom_vvl_interpol( e3v_b(:,:,:), e3vw_b(:,:,:), 'VW' )

      ! t- and w- points depth (set the isf depth as it is in the initial step)
      gdept_n(:,:,1) = 0.5_wp * e3w_n(:,:,1)
      gdepw_n(:,:,1) = 0.0_wp
      gde3w_n(:,:,1) = gdept_n(:,:,1) - sshn(:,:)
      DO jk = 2, jpk
         DO jj = 1,jpj
            DO ji = 1,jpi
              !    zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))   ! 0 everywhere tmask = wmask, ie everywhere expect at jk = mikt
                                                                 ! 1 for jk = mikt
               zcoef = (tmask(ji,jj,jk) - wmask(ji,jj,jk))
               gdepw_n(ji,jj,jk) = gdepw_n(ji,jj,jk-1) + e3t_n(ji,jj,jk-1)
               gdept_n(ji,jj,jk) =    zcoef  * ( gdepw_n(ji,jj,jk  ) + 0.5 * e3w_n(ji,jj,jk) )  &
                   &             + (1-zcoef) * ( gdept_n(ji,jj,jk-1) +       e3w_n(ji,jj,jk) ) 
               gde3w_n(ji,jj,jk) = gdept_n(ji,jj,jk) - sshn(ji,jj)
            END DO
         END DO
      END DO

      ! Local depth and Inverse of the local depth of the water
      ! -------------------------------------------------------
      hu_n(:,:) = hu_a(:,:)   ;   r1_hu_n(:,:) = r1_hu_a(:,:)
      hv_n(:,:) = hv_a(:,:)   ;   r1_hv_n(:,:) = r1_hv_a(:,:)
      !
      ht_n(:,:) = e3t_n(:,:,1) * tmask(:,:,1)
      DO jk = 2, jpkm1
         ht_n(:,:) = ht_n(:,:) + e3t_n(:,:,jk) * tmask(:,:,jk)
      END DO

      ! Output some diagnostics:
      ! ------------------------
      IF (ln_vvl_ztilde .OR. ln_vvl_layer)  CALL dom_vvl_dia( kt )

      ! write restart file
      ! ==================
      IF( lrst_oce  )   CALL dom_vvl_rst( kt, 'WRITE' )
      !
      IF( ln_timing )   CALL timing_stop('dom_vvl_sf_swp')
      !
   END SUBROUTINE dom_vvl_sf_swp


   SUBROUTINE dom_vvl_interpol( pe3_in, pe3_out, pout )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl__interpol  ***
      !!
      !! ** Purpose :   interpolate scale factors from one grid point to another
      !!
      !! ** Method  :   e3_out = e3_0 + interpolation(e3_in - e3_0)
      !!                - horizontal interpolation: grid cell surface averaging
      !!                - vertical interpolation: simple averaging
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::  pe3_in    ! input e3 to be interpolated
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::  pe3_out   ! output interpolated e3
      CHARACTER(LEN=*)                , INTENT(in   ) ::  pout      ! grid point of out scale factors
      !                                                             !   =  'U', 'V', 'W, 'F', 'UW' or 'VW'
      !
      INTEGER ::   ji, jj, jk, jkbot                                ! dummy loop indices
      INTEGER ::   nmet                                             ! horizontal interpolation method
      REAL(wp) ::  zlnwd                                            ! =1./0. when ln_wd_il = T/F
      REAL(wp) ::  ztap, zsmall                                     ! Parameters defining minimum thicknesses UVF-points
      REAL(wp) ::  zmin
      REAL(wp) ::  zdo, zup                                         ! Lower and upper interfaces depths anomalies
      REAL(wp), DIMENSION(jpi,jpj) :: zs                            ! Surface interface depth anomaly
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zw                        ! Interface depth anomaly
      !!----------------------------------------------------------------------
      !
!     nmet = 0  ! Original method (Surely wrong)
     nmet = 2  ! Interface interpolation
!      nmet = 2  ! Internal interfaces interpolation only, spread barotropic increment
!     Note that we kept surface weighted interpolation for barotropic increment to be compliant
!     with what is done in surface pressure module.
      !
      IF(ln_wd_il) THEN
        zlnwd = 1.0_wp
      ELSE
        zlnwd = 0.0_wp
      END IF
      !
      ztap   = 0._wp   ! Minimum fraction of T-point thickness at cell interfaces
      zsmall = 1.e-8_wp ! Minimum thickness at U or V points (m)
      !
      IF ( (nmet==1).OR.(nmet==2) ) THEN
         SELECT CASE ( pout )
            !
         CASE( 'U', 'V', 'F' )
            ! Compute interface depth anomaly at T-points
            !
            zw(:,:,:) =  0._wp    
            !
            DO jk=2,jpk
               zw(:,:,jk) = zw(:,:,jk-1) + pe3_in(:,:,jk-1)*tmask(:,:,jk-1)
            END DO 
            ! Interface depth anomalies:
            DO jk=1,jpkm1
               zw(:,:,jk) = zw(:,:,jk) - zw(:,:,jpk) + ht_0(:,:)
            END DO
            zw(:,:,jpk) = ht_0(:,:)
            !
            IF (nmet==2) THEN        ! Consider "internal" interfaces only
               zs(:,:) = - zw(:,:,1) ! Save surface anomaly (ssh)
               !
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     DO jk=1,jpk
                        zw(ji,jj,jk) = (zw(ji,jj,jk) + zs(ji,jj))                                          &
                                     & * ht_0(ji,jj) / (ht_0(ji,jj) + zs(ji,jj) + 1._wp - tmask(ji,jj,1))  &
                                     & * tmask(ji,jj,jk)
                     END DO
                  END DO
               END DO 
            ENDIF 
            zw(:,:,:) = (zw(:,:,:) - gdepw_0(:,:,:))*tmask(:,:,:)
            !
         END SELECT
      END IF
      !
      pe3_out(:,:,:) = 0.0_wp
      !
      SELECT CASE ( pout )    !==  type of interpolation  ==!
         !
      CASE( 'U' )                   !* from T- to U-point : hor. surface weighted mean
         IF (nmet==0) THEN
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     pe3_out(ji,jj,jk) = 0.5_wp * (  umask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd ) * r1_e1e2u(ji,jj)   &
                        &                       * (   e1e2t(ji  ,jj) * ( pe3_in(ji  ,jj,jk) - e3t_0(ji  ,jj,jk) )     &
                        &                           + e1e2t(ji+1,jj) * ( pe3_in(ji+1,jj,jk) - e3t_0(ji+1,jj,jk) ) )
                  END DO
               END DO
            END DO
         ELSE
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! Correction at last level:
                  jkbot = mbku(ji,jj)
                  zdo = 0._wp
                  DO jk=jkbot,1,-1
                     zup = 0.5_wp * ( e1e2t(ji  ,jj)*zw(ji  ,jj,jk) &
                         &          + e1e2t(ji+1,jj)*zw(ji+1,jj,jk) ) * r1_e1e2u(ji,jj)
                     !
                     ! If there is a step, taper bottom interface:
!                     IF ((hu_0(ji,jj) < 0.5_wp * ( ht_0(ji,jj) + ht_0(ji+1,jj) ) ).AND.(jk==jkbot)) THEN
!                        IF ( ht_0(ji+1,jj) < ht_0(ji,jj) ) THEN
!                           zmin = ztap * (zw(ji+1,jj,jk+1)-zw(ji+1,jj,jk))
!                        ELSE
!                           zmin = ztap * (zw(ji  ,jj,jk+1)-zw(ji  ,jj,jk))
!                        ENDIF
!                        zup = MIN(zup, zdo-zmin)
!                     ENDIF
                     zup = MIN(zup, zdo+e3u_0(ji,jj,jk)-zsmall)
                     pe3_out(ji,jj,jk) = (zdo - zup) * ( umask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd )
                     zdo = zup
                  END DO
               END DO
            END DO
         END IF
         ! 
         IF (nmet==2) THEN           ! Spread sea level anomaly
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  DO jk=1,jpk
                     pe3_out(ji,jj,jk) =       pe3_out(ji,jj,jk)                                   &
                           &               + ( pe3_out(ji,jj,jk) + e3u_0(ji,jj,jk) )               & 
                           &               / ( hu_0(ji,jj) + 1._wp - ssumask(ji,jj) )              &
                           &               * 0.5_wp * r1_e1e2u(ji,jj)                              &
                           &               * ( umask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd )        &
                           &               * ( e1e2t(ji,jj)*zs(ji,jj) + e1e2t(ji+1,jj)*zs(ji+1,jj) )       
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk( pe3_out(:,:,:), 'U', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3u_0(:,:,:)
         !
      CASE( 'V' )                   !* from T- to V-point : hor. surface weighted mean
         IF (nmet==0) THEN
            DO jk = 1, jpk
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     pe3_out(ji,jj,jk) = 0.5_wp * ( vmask(ji,jj,jk)  * (1.0_wp - zlnwd) + zlnwd ) * r1_e1e2v(ji,jj)   &
                        &                       * (   e1e2t(ji,jj  ) * ( pe3_in(ji,jj  ,jk) - e3t_0(ji,jj  ,jk) )     &
                        &                           + e1e2t(ji,jj+1) * ( pe3_in(ji,jj+1,jk) - e3t_0(ji,jj+1,jk) ) )
                  END DO
               END DO
            END DO
         ELSE
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! Correction at last level:
                  jkbot = mbkv(ji,jj)
                  zdo = 0._wp
                  DO jk=jkbot,1,-1
                     zup = 0.5_wp * ( e1e2t(ji,jj  ) * zw(ji,jj  ,jk) & 
                         &          + e1e2t(ji,jj+1) * zw(ji,jj+1,jk) ) * r1_e1e2v(ji,jj)
                     !
                     ! If there is a step, taper bottom interface:
!                     IF ((hv_0(ji,jj) < 0.5_wp * ( ht_0(ji,jj) + ht_0(ji,jj+1) ) ).AND.(jk==jkbot)) THEN
!                        IF ( ht_0(ji,jj+1) < ht_0(ji,jj) ) THEN
!                           zmin = ztap * (zw(ji,jj+1,jk+1)-zw(ji,jj+1,jk))
!                        ELSE
!                           zmin = ztap * (zw(ji  ,jj,jk+1)-zw(ji  ,jj,jk))
!                        ENDIF
!                        zup = MIN(zup, zdo-zmin)
!                     ENDIF
                     zup = MIN(zup, zdo + e3v_0(ji,jj,jk) - zsmall)
                     pe3_out(ji,jj,jk) = (zdo - zup) * ( vmask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd )
                     zdo = zup
                  END DO
               END DO
            END DO
         END IF
         !
         IF (nmet==2) THEN           ! Spread sea level anomaly
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  DO jk=1,jpk
                     pe3_out(ji,jj,jk) =       pe3_out(ji,jj,jk)                                                       &
                           &               + ( pe3_out(ji,jj,jk) + e3v_0(ji,jj,jk) )                                   & 
                           &               / ( hv_0(ji,jj) + 1._wp - ssvmask(ji,jj) )                                  &
                           &               * 0.5_wp * r1_e1e2v(ji,jj)                                                  &
                                           * ( vmask(ji,jj,jk) * (1.0_wp - zlnwd) + zlnwd )                            &
                           &               * ( e1e2t(ji,jj)*zs(ji,jj) + e1e2t(ji,jj+1)*zs(ji,jj+1) )       
                  END DO
               END DO
            END DO
            !
         ENDIF
         !
         CALL lbc_lnk( pe3_out(:,:,:), 'V', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3v_0(:,:,:)
         !
      CASE( 'F' )                   !* from T-point to F-point : hor. surface weighted mean
         IF (nmet==0) THEN
            DO jk=1,jpk
               DO jj = 1, jpjm1
                  DO ji = 1, fs_jpim1   ! vector opt.
                     pe3_out(ji,jj,jk) = 0.25_wp * (  umask(ji,jj,jk) * umask(ji,jj+1,jk) * (1.0_wp - zlnwd) + zlnwd ) &
                           &                     *    r1_e1e2f(ji,jj)                                                  &
                           &                     * (  e1e2t(ji  ,jj  ) * ( pe3_in(ji  ,jj  ,jk)-e3t_0(ji  ,jj  ,jk) )  & 
                           &                        + e1e2t(ji  ,jj+1) * ( pe3_in(ji  ,jj+1,jk)-e3t_0(ji  ,jj+1,jk) )  &
                           &                        + e1e2t(ji+1,jj  ) * ( pe3_in(ji+1,jj  ,jk)-e3t_0(ji+1,jj  ,jk) )  & 
                           &                        + e1e2t(ji+1,jj+1) * ( pe3_in(ji+1,jj+1,jk)-e3t_0(ji+1,jj+1,jk) ) )                                                  
                  END DO
               END DO
            END DO
         ELSE
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! bottom correction:
                  jkbot = MIN(mbku(ji,jj), mbku(ji,jj+1)) 
                  zdo = 0._wp
                  DO jk=jkbot,1,-1
                     zup =  0.25_wp * (  e1e2t(ji  ,jj  ) * zw(ji  ,jj  ,jk)  & 
                           &           + e1e2t(ji+1,jj  ) * zw(ji+1,jj  ,jk)  & 
                           &           + e1e2t(ji  ,jj+1) * zw(ji  ,jj+1,jk)  & 
                           &           + e1e2t(ji+1,jj+1) * zw(ji+1,jj+1,jk)  ) *    r1_e1e2f(ji,jj)
                     !
                     ! If there is a step, taper bottom interface:
!                     IF ((hf_0(ji,jj) < 0.5_wp * ( hu_0(ji,jj  ) + hu_0(ji,jj+1) ) ).AND.(jk==jkbot)) THEN
!                        IF ( hu_0(ji,jj+1) < hu_0(ji,jj) ) THEN
!                           IF ( ht_0(ji+1,jj+1) < ht_0(ji  ,jj+1) ) THEN
!                              zmin = ztap * (zw(ji+1,jj+1,jk+1)-zw(ji+1,jj+1,jk))
!                           ELSE
!                              zmin = ztap * (zw(ji  ,jj+1,jk+1)-zw(ji  ,jj+1,jk))
!                           ENDIF
!                        ELSE
!                           IF ( ht_0(ji+1,jj  ) < ht_0(ji  ,jj  ) ) THEN
!                              zmin = ztap * (zw(ji+1,jj  ,jk+1)-zw(ji+1,jj  ,jk))
!                           ELSE
!                              zmin = ztap * (zw(ji  ,jj  ,jk+1)-zw(ji  ,jj  ,jk))
!                           ENDIF
!                        ENDIF
!                        zup = MIN(zup, zdo-zmin)
!                     ENDIF
                     zup = MIN(zup, zdo+e3f_0(ji,jj,jk)-zsmall)
                     !
                     pe3_out(ji,jj,jk) = ( zdo - zup ) & 
                                      & *( umask(ji,jj,jk) * umask(ji,jj+1,jk) * (1.0_wp - zlnwd) + zlnwd )
                     zdo = zup
                  END DO
               END DO
            END DO
         END IF
         !
         IF (nmet==2) THEN           ! Spread sea level anomaly
            !
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  DO jk=1,jpk
                     pe3_out(ji,jj,jk) =       pe3_out(ji,jj,jk)                                           &
                           &               + ( pe3_out(ji,jj,jk) + e3f_0(ji,jj,jk) )                       & 
                           &               / ( hf_0(ji,jj) + 1._wp - umask(ji,jj,1)*umask(ji,jj+1,1) )     &
                           &               * 0.25_wp * r1_e1e2f(ji,jj)                                     & 
                           &               * ( umask(ji,jj,jk)*umask(ji,jj+1,jk)*(1.0_wp - zlnwd) + zlnwd )&
                           &               * ( e1e2t(ji  ,jj)*zs(ji  ,jj) + e1e2t(ji  ,jj+1)*zs(ji  ,jj+1) &
                           &                  +e1e2t(ji+1,jj)*zs(ji+1,jj) + e1e2t(ji+1,jj+1)*zs(ji+1,jj+1) )               
                  END DO
               END DO
            END DO
         END IF
         !
         CALL lbc_lnk( pe3_out(:,:,:), 'F', 1._wp )
         pe3_out(:,:,:) = pe3_out(:,:,:) + e3f_0(:,:,:)
         !
      CASE( 'W' )                   !* from T- to W-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3w_0(:,:,1) + pe3_in(:,:,1) - e3t_0(:,:,1)
         ! - ML - The use of mask in this formulea enables the special treatment of the last w-point without indirect adressing
!!gm BUG? use here wmask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3w_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( tmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )   &
               &                            * ( pe3_in(:,:,jk-1) - e3t_0(:,:,jk-1) )                               &
               &                            +            0.5_wp * ( tmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )     &
               &                            * ( pe3_in(:,:,jk  ) - e3t_0(:,:,jk  ) )
         END DO
         !
      CASE( 'UW' )                  !* from U- to UW-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3uw_0(:,:,1) + pe3_in(:,:,1) - e3u_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
!!gm BUG? use here wumask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3uw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( umask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )  &
               &                             * ( pe3_in(:,:,jk-1) - e3u_0(:,:,jk-1) )                              &
               &                             +            0.5_wp * ( umask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )    &
               &                             * ( pe3_in(:,:,jk  ) - e3u_0(:,:,jk  ) )
         END DO
         !
      CASE( 'VW' )                  !* from V- to VW-point : vertical simple mean
         !
         pe3_out(:,:,1) = e3vw_0(:,:,1) + pe3_in(:,:,1) - e3v_0(:,:,1)
         ! - ML - The use of mask in this formaula enables the special treatment of the last w- point without indirect adressing
!!gm BUG? use here wvmask in case of ISF ?  to be checked
         DO jk = 2, jpk
            pe3_out(:,:,jk) = e3vw_0(:,:,jk) + ( 1.0_wp - 0.5_wp * ( vmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd ) )  &
               &                             * ( pe3_in(:,:,jk-1) - e3v_0(:,:,jk-1) )                              &
               &                             +            0.5_wp * ( vmask(:,:,jk) * (1.0_wp - zlnwd) + zlnwd )    &
               &                             * ( pe3_in(:,:,jk  ) - e3v_0(:,:,jk  ) )
         END DO
      END SELECT
      !
   END SUBROUTINE dom_vvl_interpol


   SUBROUTINE dom_vvl_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dom_vvl_rst  ***
      !!                     
      !! ** Purpose :   Read or write VVL file in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain vertical scale factors,
      !!                they are set to the _0 values
      !!                if the restart does not contain vertical scale factors increments (z_tilde),
      !!                they are set to 0.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !
      INTEGER ::   ji, jj, jk
      INTEGER ::   id1, id2, id3, id4, id5, id6, id7     ! local integers
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ===============
         IF( ln_rstart ) THEN                   !* Read the restart file
            CALL rst_read_open                  !  open the restart file if necessary
            CALL iom_get( numror, jpdom_autoglo, 'sshn'   , sshn, ldxios = lrxios    )
            !
            id1 = iom_varid( numror, 'e3t_b', ldstop = .FALSE. )
            id2 = iom_varid( numror, 'e3t_n', ldstop = .FALSE. )
            id3 = iom_varid( numror, 'tilde_e3t_b', ldstop = .FALSE. )
            id4 = iom_varid( numror, 'tilde_e3t_n', ldstop = .FALSE. )
            id5 = iom_varid( numror, 'hdivn_lf', ldstop = .FALSE. )
            id6 = iom_varid( numror, 'un_lf', ldstop = .FALSE. )
            id7 = iom_varid( numror, 'vn_lf', ldstop = .FALSE. )
            !                             ! --------- !
            !                             ! all cases !
            !                             ! --------- !
            IF( MIN( id1, id2 ) > 0 ) THEN       ! all required arrays exist
               CALL iom_get( numror, jpdom_autoglo, 'e3t_b', e3t_b(:,:,:), ldxios = lrxios )
               CALL iom_get( numror, jpdom_autoglo, 'e3t_n', e3t_n(:,:,:), ldxios = lrxios )
               ! needed to restart if land processor not computed 
               IF(lwp) write(numout,*) 'dom_vvl_rst : e3t_b and e3t_n found in restart files'
               WHERE ( tmask(:,:,:) == 0.0_wp ) 
                  e3t_n(:,:,:) = e3t_0(:,:,:)
                  e3t_b(:,:,:) = e3t_0(:,:,:)
               END WHERE
               IF( neuler == 0 ) THEN
                  e3t_b(:,:,:) = e3t_n(:,:,:)
               ENDIF
            ELSE IF( id1 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_n not found in restart files'
               IF(lwp) write(numout,*) 'e3t_n set equal to e3t_b.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'e3t_b', e3t_b(:,:,:), ldxios = lrxios )
               e3t_n(:,:,:) = e3t_b(:,:,:)
               neuler = 0
            ELSE IF( id2 > 0 ) THEN
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_b not found in restart files'
               IF(lwp) write(numout,*) 'e3t_b set equal to e3t_n.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'e3t_n', e3t_n(:,:,:), ldxios = lrxios )
               e3t_b(:,:,:) = e3t_n(:,:,:)
               neuler = 0
            ELSE
               IF(lwp) write(numout,*) 'dom_vvl_rst WARNING : e3t_n not found in restart file'
               IF(lwp) write(numout,*) 'Compute scale factor from sshn'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               DO jk = 1, jpk
                  e3t_n(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshn(:,:) ) &
                      &                          / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                      &          + e3t_0(:,:,jk)                               * (1._wp -tmask(:,:,jk))
               END DO
               e3t_b(:,:,:) = e3t_n(:,:,:)
               neuler = 0
            ENDIF
            !                             ! ----------- !
            IF( ln_vvl_zstar ) THEN       ! z_star case !
               !                          ! ----------- !
               IF( MIN( id3, id4 ) > 0 ) THEN
                  CALL ctl_stop( 'dom_vvl_rst: z_star cannot restart from a z_tilde or layer run' )
               ENDIF
               !                          ! ----------------------- !
            ELSE                          ! z_tilde and layer cases !
               !                          ! ----------------------- !
               IF( MIN( id3, id4 ) > 0 ) THEN  ! all required arrays exist
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_b', tilde_e3t_b(:,:,:), ldxios = lrxios )
                  CALL iom_get( numror, jpdom_autoglo, 'tilde_e3t_n', tilde_e3t_n(:,:,:), ldxios = lrxios )
               ELSE                            ! one at least array is missing
                  tilde_e3t_b(:,:,:) = 0.0_wp
                  tilde_e3t_n(:,:,:) = 0.0_wp
               ENDIF
               !                          ! ------------ !
               IF( ln_vvl_ztilde ) THEN   ! z_tilde case !
                  !                       ! ------------ !
                  IF( MIN(id5, id6, id7) > 0 ) THEN  ! required arrays exist
                     CALL iom_get( numror, jpdom_autoglo, 'hdivn_lf', hdivn_lf(:,:,:,1), ldxios = lrxios )
                     CALL iom_get( numror, jpdom_autoglo, 'un_lf', un_lf(:,:,:,1), ldxios = lrxios )
                     CALL iom_get( numror, jpdom_autoglo, 'vn_lf', vn_lf(:,:,:,1), ldxios = lrxios )
                  ELSE                ! array is missing
                     hdivn_lf(:,:,:,:) = 0.0_wp
                     un_lf(:,:,:,:) = 0.0_wp
                     vn_lf(:,:,:,:) = 0.0_wp
                  ENDIF
               ENDIF
            ENDIF
            !
         ELSE                                   !* Initialize at "rest"
            !

            IF( ll_wd ) THEN   ! MJB ll_wd edits start here - these are essential 
               !
               IF( cn_cfg == 'wad' ) THEN
                  ! Wetting and drying test case
                  CALL usr_def_istate( gdept_b, tmask, tsb, ub, vb, sshb  )
                  tsn  (:,:,:,:) = tsb (:,:,:,:)       ! set now values from to before ones
                  sshn (:,:)     = sshb(:,:)
                  un   (:,:,:)   = ub  (:,:,:)
                  vn   (:,:,:)   = vb  (:,:,:)
               ELSE
                  ! if not test case
                  sshn(:,:) = -ssh_ref
                  sshb(:,:) = -ssh_ref

                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        IF( ht_0(ji,jj)-ssh_ref <  rn_wdmin1 ) THEN ! if total depth is less than min depth

                           sshb(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                           sshn(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                           ssha(ji,jj) = rn_wdmin1 - (ht_0(ji,jj) )
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF !If test case else

               ! Adjust vertical metrics for all wad
               DO jk = 1, jpk
                  e3t_n(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshn(:,:)  ) &
                    &                            / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                    &            + e3t_0(:,:,jk) * ( 1._wp - tmask(:,:,jk) )
               END DO
               e3t_b(:,:,:) = e3t_n(:,:,:)

               DO ji = 1, jpi
                  DO jj = 1, jpj
                     IF ( ht_0(ji,jj) .LE. 0.0 .AND. NINT( ssmask(ji,jj) ) .EQ. 1) THEN
                       CALL ctl_stop( 'dom_vvl_rst: ht_0 must be positive at potentially wet points' )
                     ENDIF
                  END DO 
               END DO 
               !
            ELSE
               !
               ! Just to read set ssh in fact, called latter once vertical grid
               ! is set up:
!               CALL usr_def_istate( gdept_0, tmask, tsb, ub, vb, sshb  )
!               !
!               DO jk=1,jpk
!                  e3t_b(:,:,jk) =  e3t_0(:,:,jk) * ( ht_0(:,:) + sshb(:,:) ) &
!                     &            / ( ht_0(:,:) + 1._wp -ssmask(:,:) ) * tmask(:,:,jk)
!               END DO
!               e3t_n(:,:,:) = e3t_b(:,:,:)
                sshn(:,:)=0._wp
                e3t_n(:,:,:)=e3t_0(:,:,:)
                e3t_b(:,:,:)=e3t_0(:,:,:)
               !
            END IF           ! end of ll_wd edits

            IF( ln_vvl_ztilde .OR. ln_vvl_layer) THEN
               tilde_e3t_b(:,:,:) = 0._wp
               tilde_e3t_n(:,:,:) = 0._wp
               IF( ln_vvl_ztilde ) THEN 
                  hdivn_lf(:,:,:,:) = 0._wp 
                  un_lf(:,:,:,:) = 0._wp 
                  vn_lf(:,:,:,:) = 0._wp 
               ENDIF
            END IF
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! ===================
         IF(lwp) WRITE(numout,*) '---- dom_vvl_rst ----'
         IF( lwxios ) CALL iom_swap(      cwxios_context          )
         !                                           ! --------- !
         !                                           ! all cases !
         !                                           ! --------- !
         CALL iom_rstput( kt, nitrst, numrow, 'e3t_b', e3t_b(:,:,:), ldxios = lwxios )
         CALL iom_rstput( kt, nitrst, numrow, 'e3t_n', e3t_n(:,:,:), ldxios = lwxios )
         !                                           ! ----------------------- !
         IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN  ! z_tilde and layer cases !
            !                                        ! ----------------------- !
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_b', tilde_e3t_b(:,:,:), ldxios = lwxios)
            CALL iom_rstput( kt, nitrst, numrow, 'tilde_e3t_n', tilde_e3t_n(:,:,:), ldxios = lwxios)
         END IF
         !                                           ! -------------!    
         IF( ln_vvl_ztilde ) THEN                    ! z_tilde case !
            !                                        ! ------------ !
            CALL iom_rstput( kt, nitrst, numrow, 'hdivn_lf', hdivn_lf(:,:,:,1), ldxios = lwxios)
            CALL iom_rstput( kt, nitrst, numrow, 'un_lf', un_lf(:,:,:,1), ldxios = lwxios)
            CALL iom_rstput( kt, nitrst, numrow, 'vn_lf', vn_lf(:,:,:,1), ldxios = lwxios)
         ENDIF
         !
         IF( lwxios ) CALL iom_swap(      cxios_context          )
      ENDIF
      !
   END SUBROUTINE dom_vvl_rst


   SUBROUTINE dom_vvl_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between namelist options
      !!                for vertical coordinate
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios

      NAMELIST/nam_vvl/ ln_vvl_zstar               , ln_vvl_ztilde                       , &
                      & ln_vvl_layer               , ln_vvl_ztilde_as_zstar              , &
                      & ln_vvl_zstar_at_eqtor      , ln_vvl_zstar_on_shelf               , &
                      & ln_vvl_adv_cn2             , ln_vvl_adv_fct                      , &
                      & ln_vvl_lap                 , ln_vvl_blp                          , &
                      & rn_ahe3_lap                , rn_ahe3_blp                         , &
                      & rn_rst_e3t                 , rn_lf_cutoff                        , &
                      & ln_vvl_regrid                                                    , &
                      & ln_vvl_ramp                , rn_day_ramp                         , &
                      & ln_vvl_dbg   ! not yet implemented: ln_vvl_kepe
      !!----------------------------------------------------------------------  
      !
      REWIND( numnam_ref )              ! Namelist nam_vvl in reference namelist : 
      READ  ( numnam_ref, nam_vvl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_vvl in reference namelist', lwp )
      REWIND( numnam_cfg )              ! Namelist nam_vvl in configuration namelist : Parameters of the run
      READ  ( numnam_cfg, nam_vvl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nam_vvl in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_vvl )
      !
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_vvl_ctl : choice/control of the variable vertical coordinate'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '           Namelist nam_vvl : chose a vertical coordinate'
         WRITE(numout,*) '              zstar                      ln_vvl_zstar   = ', ln_vvl_zstar
         WRITE(numout,*) '              ztilde                     ln_vvl_ztilde  = ', ln_vvl_ztilde
         WRITE(numout,*) '              layer                      ln_vvl_layer   = ', ln_vvl_layer
         WRITE(numout,*) '              ztilde as zstar   ln_vvl_ztilde_as_zstar  = ', ln_vvl_ztilde_as_zstar
         ! WRITE(numout,*) '           Namelist nam_vvl : chose kinetic-to-potential energy conservation'
         ! WRITE(numout,*) '                                         ln_vvl_kepe    = ', ln_vvl_kepe
         WRITE(numout,*) '      ztilde near the equator    ln_vvl_zstar_at_eqtor  = ', ln_vvl_zstar_at_eqtor
         WRITE(numout,*) '      ztilde on shelves          ln_vvl_zstar_on_shelf  = ', ln_vvl_zstar_on_shelf
         WRITE(numout,*) '           Namelist nam_vvl : thickness advection scheme'
         WRITE(numout,*) '              2nd order                  ln_vvl_adv_cn2 = ', ln_vvl_adv_cn2
         WRITE(numout,*) '              2nd order FCT              ln_vvl_adv_fct = ', ln_vvl_adv_fct                    
         WRITE(numout,*) '           Namelist nam_vvl : thickness diffusion scheme'
         WRITE(numout,*) '              Laplacian                  ln_vvl_lap     = ', ln_vvl_lap
         WRITE(numout,*) '              Bilaplacian                ln_vvl_blp     = ', ln_vvl_blp
         WRITE(numout,*) '              Laplacian   coefficient    rn_ahe3_lap    = ', rn_ahe3_lap
         WRITE(numout,*) '              Bilaplacian coefficient    rn_ahe3_blp    = ', rn_ahe3_blp
         WRITE(numout,*) '           Namelist nam_vvl : layers regriding'
         WRITE(numout,*) '              ln_vvl_regrid                             = ', ln_vvl_regrid
         WRITE(numout,*) '           Namelist nam_vvl : linear ramp at startup'
         WRITE(numout,*) '              ln_vvl_ramp                               = ', ln_vvl_ramp
         WRITE(numout,*) '              rn_day_ramp                               = ', rn_day_ramp
         IF( ln_vvl_ztilde_as_zstar ) THEN
            WRITE(numout,*) '           ztilde running in zstar emulation mode; '
            WRITE(numout,*) '           ignoring namelist timescale parameters and using:'
            WRITE(numout,*) '                 hard-wired : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                                         rn_rst_e3t     =    0.0'
            WRITE(numout,*) '                 hard-wired : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                                         rn_lf_cutoff   =    1.0/rdt'
         ELSE
            WRITE(numout,*) '           Namelist nam_vvl : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                                         rn_rst_e3t     = ', rn_rst_e3t
            WRITE(numout,*) '           Namelist nam_vvl : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                                         rn_lf_cutoff   = ', rn_lf_cutoff
         ENDIF
         WRITE(numout,*) '           Namelist nam_vvl : debug prints'
         WRITE(numout,*) '                                         ln_vvl_dbg     = ', ln_vvl_dbg
      ENDIF
      !
      IF ( ln_vvl_ztilde.OR.ln_vvl_layer ) THEN 
         ioptio = 0                      ! Choose one advection scheme at most
         IF( ln_vvl_adv_cn2         )        ioptio = ioptio + 1
         IF( ln_vvl_adv_fct         )        ioptio = ioptio + 1
         IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE thickness advection scheme in namelist nam_vvl' )
      ENDIF
      !
      ioptio = 0                      ! Parameter control
      IF( ln_vvl_ztilde_as_zstar )   ln_vvl_ztilde = .true.
      IF( ln_vvl_zstar           )   ioptio = ioptio + 1
      IF( ln_vvl_ztilde          )   ioptio = ioptio + 1
      IF( ln_vvl_layer           )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE vertical coordinate in namelist nam_vvl' )
      IF( .NOT. ln_vvl_zstar .AND. ln_isf ) CALL ctl_stop( 'Only vvl_zstar has been tested with ice shelf cavity' )
      !
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( ln_vvl_zstar           ) WRITE(numout,*) '      ==>>>   zstar vertical coordinate is used'
         IF( ln_vvl_ztilde          ) WRITE(numout,*) '      ==>>>   ztilde vertical coordinate is used'
         IF( ln_vvl_layer           ) WRITE(numout,*) '      ==>>>   layer vertical coordinate is used'
         IF( ln_vvl_ztilde_as_zstar ) WRITE(numout,*) '      ==>>>   to emulate a zstar coordinate'
      ENDIF
      !
      ! Use of "shelf horizon depths" should be allowed with s-z coordinates, but we restrict it to zco and zps 
      ! for the time being
      IF ( ln_sco ) THEN 
        ll_shorizd=.FALSE.
      ELSE
        ll_shorizd=.TRUE.
      ENDIF
      !
#if defined key_agrif
      IF( (.NOT.Agrif_Root()).AND.(.NOT.ln_vvl_zstar) )   CALL ctl_stop( 'AGRIF is implemented with zstar coordinate only' )
#endif
      !
   END SUBROUTINE dom_vvl_ctl

   SUBROUTINE dom_vvl_regrid( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_regrid  ***
      !!                   
      !! ** Purpose :  Ensure "well-behaved" vertical grid
      !!
      !! ** Method  :  More or less adapted from references below.
      !!regrid
      !! ** Action  :  Ensure that thickness are above a given value, spaced enough
      !!               and revert to Eulerian coordinates near the bottom.      
      !!
      !! References :  Bleck, R. and S. Benjamin, 1993: Regional Weather Prediction
      !!               with a Model Combining Terrain-following and Isentropic
      !!               coordinates. Part I: Model Description. Monthly Weather Rev.,
      !!               121, 1770-1785.
      !!               Toy, M., 2011: Incorporating Condensational Heating into a
      !!               Nonhydrostatic Atmospheric Model Based on a Hybrid Isentropic-
      !!               Sigma Vertical Coordinate. Monthly Weather Rev., 139, 2940-2954.
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in )               :: kt         ! time step

      !! * Local declarations
      INTEGER                             :: ji, jj, jk ! dummy loop indices
      LOGICAL                             :: ll_chk_bot2top, ll_chk_top2bot, ll_lapdiff_cond
      LOGICAL                             :: ll_zdiff_cond, ll_blpdiff_cond
      INTEGER                             :: jkbot
      REAL(wp)                            :: zh_min, zh_0, zh2, zdiff, zh_max, ztmph, ztmpd
      REAL(wp)                            :: zufim1, zufi, zvfjm1, zvfj, dzmin_int, dzmin_surf
      REAL(wp)                            :: zh_new, zh_old, zh_bef, ztmp, ztmp1, z2dt, zh_up, zh_dwn
      REAL(wp)                            :: zeu2, zev2, zfrch_stp, zfrch_rel, zfrac_bot, zscal_bot
      REAL(wp)                            :: zhdiff, zhdiff2, zvdiff, zhlim, zhlim2, zvlim
      REAL(wp), DIMENSION(jpi,jpj)        :: zdw, zwu, zwv
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: zwdw, zwdw_b
      !!----------------------------------------------------------------------

      IF( ln_timing )  CALL timing_start('dom_vvl_regrid')
      !
      !
      ! Some user defined parameters below:
      ll_chk_bot2top   = .TRUE.
      ll_chk_top2bot   = .TRUE.
      dzmin_int  = 1.0_wp   ! Absolute minimum depth in the interior (in meters)
      dzmin_surf = 1.0_wp   ! Absolute minimum depth at the surface (in meters)
      zfrch_stp  = 5._wp   ! Maximum fractionnal thickness change in one time step (<= 1.)
      zfrch_rel  = 0.4_wp   ! Maximum relative thickness change in the vertical (<= 1.)
      zfrac_bot  = 0.05_wp  ! Fraction of bottom level allowed to change
      zscal_bot  = 2.0_wp   ! Depth lengthscale
      ll_zdiff_cond    = .TRUE.  ! Conditionnal vertical diffusion of interfaces
         zvdiff  = 0.2_wp   ! m
         zvlim   = 0.5_wp   ! max d2h/dh
      ll_lapdiff_cond  = .TRUE.  ! Conditionnal Laplacian diffusion of interfaces
         zhdiff  = 0.01_wp  ! ad.
         zhlim   = 0.03_wp  ! ad. max lap(z)*e1
      ll_blpdiff_cond  = .TRUE.  ! Conditionnal Bilaplacian diffusion of interfaces
         zhdiff2 = 0.2_wp   ! ad.
         zhlim2  = 0.01_wp  ! ad. max bilap(z)*e1**3
      ! --------------------------------------------------------------------------------------- 
      !
      ! Set arrays determining maximum vertical displacement at the bottom:
      !--------------------------------------------------------------------
      IF ( kt==nit000 ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               jk = MIN(mbkt(ji,jj), mbkt(ji+1,jj), mbkt(ji-1,jj), mbkt(ji,jj+1), mbkt(ji,jj-1))
               jk = MIN(jk,mbkt(ji-1,jj-1), mbkt(ji-1,jj+1), mbkt(ji+1,jj+1), mbkt(ji+1,jj-1))
               i_int_bot(ji,jj) = jk
            END DO
         END DO
         dsm(:,:) = REAL( i_int_bot(:,:), wp )   ;   CALL lbc_lnk(dsm(:,:),'T',1.)
         i_int_bot(:,:) = MAX( INT( dsm(:,:) ), 1 )

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zdw(ji,jj) = MAX(ABS(ht_0(ji,jj)-ht_0(ji+1,jj))*umask(ji  ,jj,1), &
                          &     ABS(ht_0(ji,jj)-ht_0(ji-1,jj))*umask(ji-1,jj,1), &
                          &     ABS(ht_0(ji,jj)-ht_0(ji,jj+1))*vmask(ji,jj  ,1), &
                          &     ABS(ht_0(ji,jj)-ht_0(ji,jj-1))*vmask(ji,jj-1,1)  )
                zdw(ji,jj) = MAX(zscal_bot * zdw(ji,jj), rsmall )
            END DO
         END DO
         CALL lbc_lnk( zdw(:,:), 'T', 1. )

         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               dsm(ji,jj) = 1._wp/16._wp * (            zdw(ji-1,jj-1) + zdw(ji+1,jj-1)      &
                  &                           +         zdw(ji-1,jj+1) + zdw(ji+1,jj+1)      &
                  &                           + 2._wp*( zdw(ji  ,jj-1) + zdw(ji-1,jj  )      &
                  &                           +         zdw(ji+1,jj  ) + zdw(ji  ,jj+1) )    &
                  &                           + 4._wp*  zdw(ji  ,jj  )                       )
            END DO
         END DO         

         CALL lbc_lnk( dsm(:,:), 'T', 1. )   
     
         IF (ln_zps) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  jk = i_int_bot(ji,jj)
                  hsm(ji,jj) = zfrac_bot * e3w_1d(jk)
                  dsm(ji,jj) = MAX(dsm(ji,jj), 0.05_wp*ht_0(ji,jj))
               END DO
            END DO
         ELSE
            DO jj = 1, jpj
               DO ji = 1, jpi
                  jk = i_int_bot(ji,jj)
                  hsm(ji,jj) = zfrac_bot * e3w_0(ji,jj,jk)
                  dsm(ji,jj) = MAX(dsm(ji,jj), 0.05_wp*ht_0(ji,jj))
               END DO
            END DO
         ENDIF
      END IF

      ! Provisionnal interface depths:
      !-------------------------------
      zwdw(:,:,1) = 0.e0
      DO jj = 1, jpj
         DO ji = 1, jpi
            DO jk = 2, jpk
               zwdw(ji,jj,jk) =  zwdw(ji,jj,jk-1) + & 
                              & (tilde_e3t_a(ji,jj,jk-1)+e3t_0(ji,jj,jk-1)) * tmask(ji,jj,jk-1)
            END DO
         END DO
      END DO
      !
      ! Conditionnal horizontal Laplacian diffusion:
      !---------------------------------------------
      IF ( ll_lapdiff_cond ) THEN
         !
         zwdw_b(:,:,1) = 0._wp
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk=2,jpk
                  zwdw_b(ji,jj,jk) =  zwdw_b(ji,jj,jk-1) + & 
                                   & (tilde_e3t_b(ji,jj,jk-1)+e3t_0(ji,jj,jk-1)) * tmask(ji,jj,jk-1)
               END DO
            END DO
         END DO
         !
         DO jk = 2, jpkm1
            zwu(:,:) = 0._wp
            zwv(:,:) = 0._wp
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.                  
                  zwu(ji,jj) =  umask(ji,jj,jk) * e2_e1u(ji,jj) &
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji+1,jj  ,jk) )
                  zwv(ji,jj) =  vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji  ,jj+1,jk) )
               END DO
            END DO
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztmp1 = ( zwu(ji-1,jj  ) - zwu(ji,jj) &
                     &  +   zwv(ji  ,jj-1) - zwv(ji,jj) ) * r1_e1e2t(ji,jj)
                  zh2 = MAX(abs(ztmp1)-zhlim*SQRT(r1_e1e2t(ji,jj)), 0._wp)
                  ztmp = SIGN(zh2, ztmp1)
                  zeu2 = zhdiff * e1e2t(ji,jj)*e1e2t(ji,jj)/(e1t(ji,jj)*e1t(ji,jj) + e2t(ji,jj)*e2t(ji,jj))
                  zwdw(ji,jj,jk) = zwdw(ji,jj,jk) + zeu2 * ztmp *  tmask(ji,jj,jk)
               END DO
            END DO
         END DO         
         !
      ENDIF

      ! Conditionnal horizontal Bilaplacian diffusion:
      !-----------------------------------------------
      IF ( ll_blpdiff_cond ) THEN
         !
         zwdw_b(:,:,1) = 0._wp
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 2,jpkm1
                  zwdw_b(ji,jj,jk) =  zwdw_b(ji,jj,jk-1) + & 
                                   & (tilde_e3t_b(ji,jj,jk-1)+e3t_0(ji,jj,jk-1)) * tmask(ji,jj,jk-1)
               END DO
            END DO
         END DO
         !
         DO jk = 2, jpkm1
            zwu(:,:) = 0._wp
            zwv(:,:) = 0._wp
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.                  
                  zwu(ji,jj) =  umask(ji,jj,jk) * e2_e1u(ji,jj) &
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji+1,jj  ,jk) )
                  zwv(ji,jj) =  vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji  ,jj+1,jk) )
               END DO
            END DO
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwdw_b(ji,jj,jk) = -( (zwu(ji-1,jj  ) - zwu(ji,jj)) &
                                &  +    (zwv(ji  ,jj-1) - zwv(ji,jj)) ) * r1_e1e2t(ji,jj)
               END DO
            END DO
         END DO   
         !
         CALL lbc_lnk( zwdw_b(:,:,:), 'T', 1. )      
         !
         DO jk = 2, jpkm1
            zwu(:,:) = 0._wp
            zwv(:,:) = 0._wp
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.                  
                  zwu(ji,jj) =  umask(ji,jj,jk) * e2_e1u(ji,jj) &
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji+1,jj  ,jk) )
                  zwv(ji,jj) =  vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                             &  * ( zwdw_b(ji,jj,jk) - zwdw_b(ji  ,jj+1,jk) )
               END DO
            END DO
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztmp1 = ( (zwu(ji-1,jj  ) - zwv(ji,jj)) &
                     &  +   (zwv(ji  ,jj-1) - zwv(ji,jj)) ) * r1_e1e2t(ji,jj)
                  zh2 = MAX(abs(ztmp1)-zhlim2*SQRT(r1_e1e2t(ji,jj))*r1_e1e2t(ji,jj), 0._wp)
                  ztmp = SIGN(zh2, ztmp1)
                  zeu2 = zhdiff2 * e1e2t(ji,jj)*e1e2t(ji,jj) / 16._wp
                  zwdw(ji,jj,jk) = zwdw(ji,jj,jk) + zeu2 * ztmp *  tmask(ji,jj,jk)
               END DO
            END DO
         END DO   
         !
      ENDIF

      ! Conditionnal vertical diffusion:
      !---------------------------------
      IF ( ll_zdiff_cond ) THEN
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.    
                  ztmp  = -( (tilde_e3t_b(ji,jj,jk-1)+e3t_0(ji,jj,jk-1))*tmask(ji,jj,jk-1) & 
                            -(tilde_e3t_b(ji,jj,jk  )+e3t_0(ji,jj,jk  ))*tmask(ji,jj,jk  ) ) 
                  ztmp1 = 0.5_wp * ( tilde_e3t_b(ji,jj,jk-1) + e3t_0(ji,jj,jk-1) & 
                        &           +tilde_e3t_b(ji,jj,jk  ) + e3t_0(ji,jj,jk  ) )
                  zh2   = MAX(abs(ztmp)-zvlim*ztmp1, 0._wp)
                  ztmp  = SIGN(zh2, ztmp)
                  IF ((jk==mbkt(ji,jj)).AND.(ln_zps)) ztmp=0.e0
                  zwdw(ji,jj,jk) = zwdw(ji,jj,jk) + zvdiff * ztmp * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF
      !
      ! Check grid from the bottom to the surface
      !------------------------------------------
      IF ( ll_chk_bot2top ) THEN
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               jkbot = mbkt(ji,jj)                   
               DO jk = jkbot,2,-1
                  !
                  zh_0   = e3t_0(ji,jj,jk)
                  zh_bef = MIN(tilde_e3t_b(ji,jj,jk) + zh_0, tilde_e3t_b(ji,jj,jk-1) + e3t_0(ji,jj,jk-1))
                  zh_old = zwdw(ji,jj,jk+1) - zwdw(ji,jj,jk)
                  zh_min = MIN(zh_0/3._wp, dzmin_int)
!                  zh_min = MAX(zh_min, zh_min-e3t_a(ji,jj,jk)+e3t_0(ji,jj,jk))
                  ! 
                  ! Set maximum and minimum vertical excursions   
                  ztmph = hsm(ji,jj)
                  ztmpd = dsm(ji,jj)
                  zh2   = ztmph * exp(-(gdepw_0(ji,jj,jk)-gdepw_0(ji,jj,i_int_bot(ji,jj)))/ztmpd)
!                  zh2   = ztmph * exp(-(gdepw_0(ji,jj,jk)-gdepw_0(ji,jj,i_int_bot(ji,jj)+1))/ztmpd)
                  zh2 = MAX(zh2,0.001_wp) ! Extend tolerance a bit for stability reasons (to be explored)
                  zdiff = cush_max(gdepw_0(ji,jj,jk)-zwdw(ji,jj,jk), zh2 )
                  zwdw(ji,jj,jk) = MAX(zwdw(ji,jj,jk), gdepw_0(ji,jj,jk) - zdiff)
                  zdiff = cush_max(zwdw(ji,jj,jk)-gdepw_0(ji,jj,jk), zh2 )
                  zwdw(ji,jj,jk) = MIN(zwdw(ji,jj,jk), gdepw_0(ji,jj,jk) + zdiff)
                  !
                  ! New layer thickness:
                  zh_new =  zwdw(ji,jj,jk+1) - zwdw(ji,jj,jk)
                  !
                  ! Ensure minimum layer thickness:                  
!                  zh_new = MAX((1._wp-zfrch_stp)*zh_bef, zh_new)
!                  zh_new = cush(zh_new, zh_min)
                  zh_new = MAX(zh_new, zh_min)           
                  !
                  ! Final flux:
                  zdiff = (zh_new - zh_old) * tmask(ji,jj,jk)
                  !
                  ! Limit thickness change in 1 time step: 
!                  ztmp = MIN( ABS(zdiff), zfrch_stp*zh_bef )
!                  zdiff = SIGN(ztmp, zh_new - zh_old)
                  zh_new = zdiff + zh_old
                  !
                  zwdw(ji,jj,jk) = zwdw(ji,jj,jk+1) - zh_new       
               END DO    
            END DO 
         END DO
      END IF    
      !
      ! Check grid from the surface to the bottom 
      !------------------------------------------ 
      IF ( ll_chk_top2bot ) THEN      
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               jkbot = mbkt(ji,jj)   
               DO jk = 1, jkbot-1
                  !
                  zh_0   = e3t_0(ji,jj,jk)
                  zh_bef = MIN(tilde_e3t_b(ji,jj,jk) + zh_0, tilde_e3t_b(ji,jj,jk+1) + e3t_0(ji,jj,jk+1))
                  zh_old = zwdw(ji,jj,jk+1) - zwdw(ji,jj,jk)
                  zh_min = MIN(zh_0/3._wp, dzmin_int)
!                  zh_min = MAX(zh_min, zh_min-e3t_a(ji,jj,jk)+e3t_0(ji,jj,jk))
                  !
!                  zwdw(ji,jj,jk+1) = MAX(zwdw(ji,jj,jk+1), REAL(jk)*dzmin_surf)
                  !
                  ! New layer thickness:
                  zh_new =  zwdw(ji,jj,jk+1) - zwdw(ji,jj,jk)
                  !
                  ! Ensure minimum layer thickness:
!                  zh_new = MAX((1._wp-zfrch_stp)*zh_bef, zh_new)
!                  zh_new = cush(zh_new, zh_min)
                  zh_new = MAX(zh_new, zh_min)   
                  !
                  ! Final flux:
                  zdiff = (zh_new -zh_old) * tmask(ji,jj,jk)
                  ! 
                  ! Limit flux:                 
!                  ztmp = MIN( ABS(zdiff), zfrch_stp*zh_bef )
!                  zdiff = SIGN(ztmp, zh_new - zh_old)
                  zh_new = zdiff + zh_old
                  !
                  zwdw(ji,jj,jk+1) = zwdw(ji,jj,jk) + zh_new
               END DO
               !               
            END DO
         END DO
      ENDIF
      !
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            DO jk = 1, jpkm1
               tilde_e3t_a(ji,jj,jk) =  (zwdw(ji,jj,jk+1)-zwdw(ji,jj,jk)-e3t_0(ji,jj,jk)) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      ! 
      !
      IF( ln_timing )  CALL timing_stop('dom_vvl_regrid')
      !
   END SUBROUTINE dom_vvl_regrid

   FUNCTION cush(hin, hmin)  RESULT(hout)
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION cush  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(wp), INTENT(in) ::  hin, hmin
      REAL(wp)             ::  hout, zx, zh_cri
      !!----------------------------------------------------------------------
      zh_cri = 3._wp * hmin
      !
      IF ( hin<=0._wp ) THEN
         hout = hmin
      !
      ELSEIF ( (hin>0._wp).AND.(hin<=zh_cri) ) THEN
         zx = hin/zh_cri
         hout = hmin * (1._wp + zx + zx*zx)      
      !
      ELSEIF ( hin>zh_cri ) THEN
         hout = hin
      !
      ENDIF
      !
   END FUNCTION cush

   FUNCTION cush_max(hin, hmax)  RESULT(hout)
      !!----------------------------------------------------------------------
      !!                 ***  FUNCTION cush  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(wp), INTENT(in) ::  hin, hmax
      REAL(wp)             ::  hout, hmin, zx, zh_cri
      !!----------------------------------------------------------------------
      hmin = 0.1_wp * hmax
      zh_cri = 3._wp * hmin
      !
      IF ( (hin>=(hmax-zh_cri)).AND.(hin<=(hmax-hmin))) THEN
         zx = (hmax-hin)/zh_cri
         hout = hmax - hmin * (1._wp + zx + zx*zx)
         !
      ELSEIF ( hin>(hmax-zh_cri) ) THEN
         hout = hmax - hmin
         !
      ELSE
         hout = hin
         !
      ENDIF
      !
   END FUNCTION cush_max

   SUBROUTINE dom_vvl_adv_fct( kt, pta, uin, vin )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl_adv_fct  ***
      !! 
      !! **  Purpose :  Do thickness advection
      !!
      !! **  Method  :  FCT scheme to ensure positivity 
      !!
      !! **  Action  : - Update pta thickness tendency and diffusive fluxes
      !!               - this is the total trend, hence it does include sea level motions
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta    ! thickness baroclinic trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   uin, vin  ! input velocities
      !
      INTEGER  ::   ji, jj, jk, ib, ib_bdy               ! dummy loop indices  
      INTEGER  ::   ikbu, ikbv, ibot
      REAL(wp) ::   z2dtt, zbtr, ztra        ! local scalar
      REAL(wp) ::   zdi, zdj, zmin           !   -      -
      REAL(wp) ::   zfp_ui, zfp_vj           !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj           !   -      -
      REAL(wp) ::   zfp_hi, zfp_hj           !   -      -
      REAL(wp) ::   zfm_hi, zfm_hj           !   -      -
      REAL(wp) ::   ztout,  ztin, zfac       !   -      -
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwx, zwy, zwi
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('dom_vvl_adv_fct')
      !
      !
      ! 1. Initializations
      ! ------------------
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtt =  rdt
      ELSE
         z2dtt = 2.0_wp * rdt
      ENDIF
      !
      zwi(:,:,:) = 0.e0
      zwx(:,:,:) = 0.e0 
      zwy(:,:,:) = 0.e0
      !
      !
      ! 2. upstream advection with initial mass fluxes & intermediate update
      ! --------------------------------------------------------------------
      IF ( ll_shorizd ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  !
                  zfp_hi = MAX(hu_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hi = MIN(zfp_hi, e3t_b(ji  ,jj  ,jk))
                  zfp_hi = 0.5_wp *(zfp_hi + SIGN(zfp_hi, zfp_hi-hsmall) ) 
                  !
                  zfm_hi = MAX(hu_b(ji,jj) - gdepw_b(ji+1,jj  ,jk), 0._wp)
                  zfm_hi = MIN(zfm_hi, e3t_b(ji+1,jj  ,jk))
                  zfm_hi = 0.5_wp *(zfm_hi + SIGN(zfm_hi, zfm_hi-hsmall) ) 
                  ! 
                  zfp_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hj = MIN(zfp_hj, e3t_b(ji  ,jj  ,jk))
                  zfp_hj = 0.5_wp *(zfp_hj + SIGN(zfp_hj, zfp_hj-hsmall) ) 
                  !
                  zfm_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj+1,jk), 0._wp)
                  zfm_hj = MIN(zfm_hj, e3t_b(ji  ,jj+1,jk))
                  zfm_hj = 0.5_wp *(zfm_hj + SIGN(zfm_hj, zfm_hj-hsmall) )
                  !
                  zfp_ui = uin(ji,jj,jk) + ABS( uin(ji,jj,jk) )
                  zfm_ui = uin(ji,jj,jk) - ABS( uin(ji,jj,jk) )
                  zfp_vj = vin(ji,jj,jk) + ABS( vin(ji,jj,jk) )
                  zfm_vj = vin(ji,jj,jk) - ABS( vin(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5 * e2u(ji,jj) * ( zfp_ui * zfp_hi + zfm_ui * zfm_hi ) * umask(ji,jj,jk)
                  zwy(ji,jj,jk) = 0.5 * e1v(ji,jj) * ( zfp_vj * zfp_hj + zfm_vj * zfm_hj ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  !
                  zfp_hi = e3t_b(ji  ,jj  ,jk)
                  zfm_hi = e3t_b(ji+1,jj  ,jk)             
                  zfp_hj = e3t_b(ji  ,jj  ,jk)
                  zfm_hj = e3t_b(ji  ,jj+1,jk)
                  !
                  zfp_ui = uin(ji,jj,jk) + ABS( uin(ji,jj,jk) )
                  zfm_ui = uin(ji,jj,jk) - ABS( uin(ji,jj,jk) )
                  zfp_vj = vin(ji,jj,jk) + ABS( vin(ji,jj,jk) )
                  zfm_vj = vin(ji,jj,jk) - ABS( vin(ji,jj,jk) )
                  zwx(ji,jj,jk) = 0.5 * e2u(ji,jj) * ( zfp_ui * zfp_hi + zfm_ui * zfm_hi ) * umask(ji,jj,jk)
                  zwy(ji,jj,jk) = 0.5 * e1v(ji,jj) * ( zfp_vj * zfp_hj + zfm_vj * zfm_hj ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! total advective trend
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr = r1_e1e2t(ji,jj)
               ! total intermediate advective trends
               ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )  &
                  &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  )  )
               !
               ! update and guess with monotonic sheme
               pta(ji,jj,jk) =   pta(ji,jj,jk) + ztra
               zwi(ji,jj,jk) = (e3t_b(ji,jj,jk) + z2dtt * ztra ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL lbc_lnk( zwi, 'T', 1. )  

      IF ( ln_bdy ) THEN
         DO ib_bdy=1, nb_bdy
            DO ib = 1, idx_bdy(ib_bdy)%nblenrim(1)
               ji = idx_bdy(ib_bdy)%nbi(ib,1)
               jj = idx_bdy(ib_bdy)%nbj(ib,1)
               DO jk = 1, jpkm1  
                  zwi(ji,jj,jk) = e3t_a(ji,jj,jk)
               END DO
            END DO 
         END DO
      ENDIF

!      IF ( ln_vvl_dbg ) THEN
!         zmin = MINVAL( zwi(:,:,:), mask = tmask(:,:,:) == 1.e0 ) 
!         IF( lk_mpp )   CALL mpp_min( zmin )
!         IF( zmin < 0._wp) THEN
!            IF(lwp) CALL ctl_warn('vvl_adv: CFL issue here')
!            IF(lwp) WRITE(numout,*) zmin
!         ENDIF
!      ENDIF

      ! 3. antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      ! antidiffusive flux on i and j
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zwx(ji,jj,jk) = (e2u(ji,jj) * uin(ji,jj,jk) * e3u_n(ji,jj,jk) &
                                & - zwx(ji,jj,jk)) * umask(ji,jj,jk)
               zwy(ji,jj,jk) = (e1v(ji,jj) * vin(ji,jj,jk) * e3v_n(ji,jj,jk) & 
                                & - zwy(ji,jj,jk)) * vmask(ji,jj,jk)
               !
               ! Update advective fluxes
               un_td(ji,jj,jk) = un_td(ji,jj,jk) - zwx(ji,jj,jk)
               vn_td(ji,jj,jk) = vn_td(ji,jj,jk) - zwy(ji,jj,jk)
            END DO
         END DO
      END DO
      
      CALL lbc_lnk_multi( zwx, 'U', -1., zwy, 'V', -1. )     !* local domain boundaries

      ! 4. monotonicity algorithm
      ! -------------------------
      CALL nonosc_2d( e3t_b(:,:,:), zwx, zwy, zwi, z2dtt )

      ! 5. final trend with corrected fluxes
      ! ------------------------------------
      !
      ! Update advective fluxes
      un_td(:,:,:) = (un_td(:,:,:) + zwx(:,:,:))*umask(:,:,:)
      vn_td(:,:,:) = (vn_td(:,:,:) + zwy(:,:,:))*vmask(:,:,:)
      !
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.  
               !
               zbtr = r1_e1e2t(ji,jj)
               ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
               ! add them to the general tracer trends
               pta(ji,jj,jk) = pta(ji,jj,jk) + ztra
               zwi(ji,jj,jk) = (e3t_b(ji,jj,jk) + z2dtt * pta(ji,jj,jk)* bdytmask(ji,jj) ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      IF ( ln_vvl_dbg ) THEN
         zmin = MINVAL( zwi(:,:,:), mask = tmask(:,:,:) == 1.e0 ) 
         IF( lk_mpp )   CALL mpp_min( zmin )
         IF( zmin < 0._wp) THEN
            IF(lwp) CALL ctl_warn('vvl_adv: CFL issue here')
            IF(lwp) WRITE(numout,*) zmin
         ENDIF
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('dom_vvl_adv_fct')
      !
   END SUBROUTINE dom_vvl_adv_fct

   SUBROUTINE dom_vvl_ups_cor( kt, pta, uin, vin ) 
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_vvl_adv_fct  ***
      !! 
      !! **  Purpose :  Correct for addionnal barotropic fluxes 
      !!                in the upstream direction
      !!
      !! **  Method  :   
      !!
      !! **  Action  : - Update diffusive fluxes uin, vin
      !!               - Remove divergence from thickness tendency
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta       ! thickness baroclinic trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   uin, vin  ! input fluxes
      INTEGER  ::   ji, jj, jk               ! dummy loop indices  
      INTEGER  ::   ikbu, ikbv, ibot
      REAL(wp) ::   zbtr, ztra               ! local scalar
      REAL(wp) ::   zdi, zdj                 !   -      -
      REAL(wp) ::   zfp_hi, zfp_hj           !   -      -
      REAL(wp) ::   zfm_hi, zfm_hj           !   -      -
      REAL(wp) ::   zfp_ui, zfp_vj           !   -      -
      REAL(wp) ::   zfm_ui, zfm_vj           !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zbu, zbv, zhu_b, zhv_b
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwx, zwy
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('dom_vvl_ups_cor')
      !
      ! Compute barotropic flux difference:
      zbu(:,:) = 0.e0
      zbv(:,:) = 0.e0
      DO jj = 1, jpj
         DO ji = 1, jpi   ! vector opt.
            DO jk = 1, jpkm1
               zbu(ji,jj) = zbu(ji,jj) - uin(ji,jj,jk) * umask(ji,jj,jk)
               zbv(ji,jj) = zbv(ji,jj) - vin(ji,jj,jk) * vmask(ji,jj,jk)
            END DO
         END DO
      ENDDO
 
      ! Compute upstream depths:
      zhu_b(:,:) = 0.e0
      zhv_b(:,:) = 0.e0

      IF ( ll_shorizd ) THEN
         ! Correct bottom value
         ! considering "shelf horizon depth"
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zdi = 0.5_wp + 0.5_wp * SIGN(1._wp, zbu(ji,jj))
               zdj = 0.5_wp + 0.5_wp * SIGN(1._wp, zbv(ji,jj))
               DO jk=1, jpkm1
                  zfp_hi = MAX(hu_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hi = MIN(zfp_hi, e3t_b(ji  ,jj  ,jk))
                  zfp_hi = 0.5_wp *(zfp_hi + SIGN(zfp_hi, zfp_hi-hsmall) ) 
                  !
                  zfm_hi = MAX(hu_b(ji,jj) - gdepw_b(ji+1,jj  ,jk), 0._wp)
                  zfm_hi = MIN(zfm_hi, e3t_b(ji+1,jj  ,jk))
                  zfm_hi = 0.5_wp *(zfm_hi + SIGN(zfm_hi, zfm_hi-hsmall) )
                  ! 
                  zfp_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hj = MIN(zfp_hj, e3t_b(ji  ,jj  ,jk))
                  zfp_hj = 0.5_wp *(zfp_hj + SIGN(zfp_hj, zfp_hj-hsmall) ) 
                  !
                  zfm_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj+1,jk), 0._wp)
                  zfm_hj = MIN(zfm_hj, e3t_b(ji  ,jj+1,jk))
                  zfm_hj = 0.5_wp *(zfm_hj + SIGN(zfm_hj, zfm_hj-hsmall) )
                  !
                  zhu_b(ji,jj) = zhu_b(ji,jj) + (         zdi  * zfp_hi &
                               &                 + (1._wp-zdi) * zfm_hi &
                               &                ) * umask(ji,jj,jk)
                  zhv_b(ji,jj) = zhv_b(ji,jj) + (         zdj  * zfp_hj &
                               &                 + (1._wp-zdj) * zfm_hj &
                               &                ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zdi = 0.5_wp + 0.5_wp * SIGN(1._wp, zbu(ji,jj)) 
               zdj = 0.5_wp + 0.5_wp * SIGN(1._wp, zbv(ji,jj))
               DO jk = 1, jpkm1
                  zfp_hi = e3t_b(ji  ,jj  ,jk)
                  zfm_hi = e3t_b(ji+1,jj  ,jk)
                  zfp_hj = e3t_b(ji  ,jj  ,jk)
                  zfm_hj = e3t_b(ji  ,jj+1,jk)
                  !
                  zhu_b(ji,jj) = zhu_b(ji,jj) + (          zdi  * zfp_hi &
                               &                  + (1._wp-zdi) * zfm_hi &
                               &                ) * umask(ji,jj,jk)
                  !
                  zhv_b(ji,jj) = zhv_b(ji,jj) + (          zdj  * zfp_hj &
                               &                  + (1._wp-zdj) * zfm_hj &
                               &                ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      CALL lbc_lnk_multi( zhu_b(:,:), 'U', 1., zhv_b(:,:), 'V', 1. )     !* local domain boundaries

      ! Corrective barotropic velocity (times hor. scale factor)
      zbu(:,:) = zbu(:,:)/ (zhu_b(:,:)*umask(:,:,1)+1._wp-umask(:,:,1))
      zbv(:,:) = zbv(:,:)/ (zhv_b(:,:)*vmask(:,:,1)+1._wp-vmask(:,:,1))
      
      ! Set corrective fluxes in upstream direction:
      !
      zwx(:,:,:) = 0.e0
      zwy(:,:,:) = 0.e0

      IF ( ll_shorizd ) THEN
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ! upstream scheme
               zfp_ui = zbu(ji,jj) + ABS( zbu(ji,jj) )
               zfm_ui = zbu(ji,jj) - ABS( zbu(ji,jj) )
               zfp_vj = zbv(ji,jj) + ABS( zbv(ji,jj) )
               zfm_vj = zbv(ji,jj) - ABS( zbv(ji,jj) )
               DO jk = 1, jpkm1
                  zfp_hi = MAX(hu_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hi = MIN(e3t_b(ji  ,jj  ,jk), zfp_hi)
                  zfp_hi = 0.5_wp *(zfp_hi + SIGN(zfp_hi, zfp_hi-hsmall) ) 
                  !
                  zfm_hi = MAX(hu_b(ji,jj) - gdepw_b(ji+1,jj  ,jk), 0._wp)
                  zfm_hi = MIN(e3t_b(ji+1,jj  ,jk), zfm_hi)
                  zfm_hi = 0.5_wp *(zfm_hi + SIGN(zfm_hi, zfm_hi-hsmall) ) 
                  !
                  zfp_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj  ,jk), 0._wp)
                  zfp_hj = MIN(e3t_b(ji  ,jj  ,jk), zfp_hj) 
                  zfp_hj = 0.5_wp *(zfp_hj + SIGN(zfp_hj, zfp_hj-hsmall) ) 
                  !
                  zfm_hj = MAX(hv_b(ji,jj) - gdepw_b(ji  ,jj+1,jk), 0._wp)
                  zfm_hj = MIN(e3t_b(ji  ,jj+1,jk), zfm_hj)
                  zfm_hj = 0.5_wp *(zfm_hj + SIGN(zfm_hj, zfm_hj-hsmall) ) 
                  !
                  zwx(ji,jj,jk) = 0.5 * ( zfp_ui * zfp_hi + zfm_ui * zfm_hi ) * umask(ji,jj,jk)
                  zwy(ji,jj,jk) = 0.5 * ( zfp_vj * zfp_hj + zfm_vj * zfm_hj ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ! upstream scheme
               zfp_ui = zbu(ji,jj) + ABS( zbu(ji,jj) )
               zfm_ui = zbu(ji,jj) - ABS( zbu(ji,jj) )
               zfp_vj = zbv(ji,jj) + ABS( zbv(ji,jj) )
               zfm_vj = zbv(ji,jj) - ABS( zbv(ji,jj) )
               DO jk = 1, jpkm1
                  zfp_hi = e3t_b(ji  ,jj  ,jk)
                  zfm_hi = e3t_b(ji+1,jj  ,jk)
                  zfp_hj = e3t_b(ji  ,jj  ,jk)
                  zfm_hj = e3t_b(ji  ,jj+1,jk)
                  !
                  zwx(ji,jj,jk) = 0.5 * ( zfp_ui * zfp_hi + zfm_ui * zfm_hi ) * umask(ji,jj,jk)
                  zwy(ji,jj,jk) = 0.5 * ( zfp_vj * zfp_hj + zfm_vj * zfm_hj ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      CALL lbc_lnk_multi( zwx, 'U', -1., zwy, 'V', -1. )     !* local domain boundaries

      uin(:,:,:) = uin(:,:,:) + zwx(:,:,:)
      vin(:,:,:) = vin(:,:,:) + zwy(:,:,:)
      !
      ! Update trend with corrective fluxes:
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.  
               !
               zbtr = r1_e1e2t(ji,jj)
               ! total advective trends
               ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk  )   &
                  &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk  ) )
               ! add them to the general tracer trends
               pta(ji,jj,jk) = pta(ji,jj,jk) + ztra
            END DO
         END DO
      END DO
      !
      IF( ln_timing )  CALL timing_stop('dom_vvl_ups_cor')
      !
   END SUBROUTINE dom_vvl_ups_cor

   SUBROUTINE nonosc_2d( pbef, paa, pbb, paft, p2dt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc_2d  ***
      !!     
      !! **  Purpose :   compute monotonic thickness fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!----------------------------------------------------------------------
      !
      !!----------------------------------------------------------------------
      REAL(wp)                         , INTENT(in   ) ::   p2dt            ! vertical profile of tracer time-step
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(in   ) ::   pbef, paft      ! before & after field
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT(inout) ::   paa, pbb        ! monotonic fluxes in the 3 directions
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt   ! local scalars
      REAL(wp) ::   zau, zbu, zcu, zav, zbv, zcv, zup, zdo            !   -      -
      REAL(wp) ::   zupip1, zupim1, zupjp1, zupjm1, zupb, zupa
      REAL(wp) ::   zdoip1, zdoim1, zdojp1, zdojm1, zdob, zdoa
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zbetup, zbetdo, zbup, zbdo
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('nonosc2')
      !
      zbig  = 1.e+40_wp
      zrtrn = 1.e-15_wp
      zbetup(:,:,jpk) = 0._wp   ;   zbetdo(:,:,jpk) = 0._wp


      ! Search local extrema
      ! --------------------
      ! max/min of pbef & paft with large negative/positive value (-/+zbig) inside land
      zbup = MAX( pbef * tmask - zbig * ( 1.e0 - tmask ),   &
         &        paft * tmask - zbig * ( 1.e0 - tmask )  )
      zbdo = MIN( pbef * tmask + zbig * ( 1.e0 - tmask ),   &
         &        paft * tmask + zbig * ( 1.e0 - tmask )  )

      DO jk = 1, jpkm1
         z2dtt = p2dt
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               ! search maximum in neighbourhood
               zup = MAX(  zbup(ji  ,jj  ,jk  ),   &
                  &        zbup(ji-1,jj  ,jk  ), zbup(ji+1,jj  ,jk  ),   &
                  &        zbup(ji  ,jj-1,jk  ), zbup(ji  ,jj+1,jk  ))

               ! search minimum in neighbourhood
               zdo = MIN(  zbdo(ji  ,jj  ,jk  ),   &
                  &        zbdo(ji-1,jj  ,jk  ), zbdo(ji+1,jj  ,jk  ),   &
                  &        zbdo(ji  ,jj-1,jk  ), zbdo(ji  ,jj+1,jk  ))

               ! positive part of the flux
               zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   

               ! negative part of the flux
               zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )

               ! up & down beta terms
               zbt = e1t(ji,jj) * e2t(ji,jj) / z2dtt
               zbetup(ji,jj,jk) = ( zup            - paft(ji,jj,jk) ) / ( zpos + zrtrn ) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zdo            ) / ( zneg + zrtrn ) * zbt
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( zbetup, 'T', 1. , zbetdo, 'T', 1. )   ! lateral boundary cond. (unchanged sign)

      ! 3. monotonic flux in the i & j direction (paa & pbb)
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zau = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
               zbu = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
               zcu =       ( 0.5  + SIGN( 0.5 , paa(ji,jj,jk) ) )
               paa(ji,jj,jk) = paa(ji,jj,jk) * ( zcu * zau + ( 1.e0 - zcu) * zbu )

               zav = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
               zbv = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
               zcv =       ( 0.5  + SIGN( 0.5 , pbb(ji,jj,jk) ) )
               pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zcv * zav + ( 1.e0 - zcv) * zbv )
            END DO
         END DO
      END DO
      CALL lbc_lnk_multi( paa, 'U', -1., pbb, 'V', -1. )     !* local domain boundaries
      !
      IF( ln_timing )  CALL timing_stop('nonosc2')
      !
   END SUBROUTINE nonosc_2d

   SUBROUTINE dom_vvl_dia( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_vvl_dia  ***
      !!                   
      !! ** Purpose :  Output some diagnostics in ztilde/zlayer cases
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in )               :: kt       ! time step
      !! * Local declarations
      INTEGER                             :: ji,jj,jk,jkbot ! dummy loop indices
      REAL(wp)                            :: ztmp1
      REAL(wp), DIMENSION(4)              :: zr1
      REAL(wp), DIMENSION(jpi,jpj       ) :: zout2d
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: zwdw, zout
      !!----------------------------------------------------------------------
      IF( ln_timing )  CALL timing_start('dom_vvl_dia')
      !
      ! Compute internal interfaces depths:
      !------------------------------------
      IF ( iom_use("dh_tilde").OR.iom_use("depw_tilde").OR.iom_use("stiff_tilde")) THEN
         zwdw(:,:,1) = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 2, jpkm1
                  zwdw(ji,jj,jk) =  zwdw(ji,jj,jk-1) + & 
                                 & (tilde_e3t_n(ji,jj,jk-1)+e3t_0(ji,jj,jk-1)) * tmask(ji,jj,jk-1)
               END DO
            END DO
         END DO
      ENDIF
      !
      ! Output interface depth anomaly:
      ! -------------------------------
      IF ( iom_use("depw_tilde") ) CALL iom_put( "depw_tilde", (zwdw(:,:,:)-gdepw_0(:,:,:))*tmask(:,:,:) )
      !
      ! Output grid stiffness (T-points):
      ! ---------------------------------
      IF ( iom_use("stiff_tilde"  ) ) THEN
         zr1(:)   = 0.e0
         zout(:,:,:) = 0.e0   
         zout2d(:,:) = 0.e0   
         DO ji = 2, jpim1
            DO jj = 2, jpjm1
               ! Exclude last level because of partial bottom cells
               jkbot = MIN(mbkt(ji,jj)-1,mbkt(ji-1,jj)-1,mbkt(ji+1,jj)-1,mbkt(ji,jj-1)-1,mbkt(ji,jj+1)-1)
               DO jk = 1, jkbot
                  zr1(1) = umask(ji-1,jj  ,jk) *abs( (zwdw(ji  ,jj  ,jk  )-zwdw(ji-1,jj  ,jk  )  & 
                       &                             +zwdw(ji  ,jj  ,jk+1)-zwdw(ji-1,jj  ,jk+1)) &
                       &                            /(zwdw(ji  ,jj  ,jk  )+zwdw(ji-1,jj  ,jk  )  &
                       &                             -zwdw(ji  ,jj  ,jk+1)-zwdw(ji-1,jj  ,jk+1) + rsmall) )
                  zr1(2) = umask(ji  ,jj  ,jk) *abs( (zwdw(ji+1,jj  ,jk  )-zwdw(ji  ,jj  ,jk  )  &
                       &                             +zwdw(ji+1,jj  ,jk+1)-zwdw(ji  ,jj  ,jk+1)) &
                       &                            /(zwdw(ji+1,jj  ,jk  )+zwdw(ji  ,jj  ,jk  )  &
                       &                             -zwdw(ji+1,jj  ,jk+1)-zwdw(ji  ,jj  ,jk+1) + rsmall) )
                  zr1(3) = vmask(ji  ,jj  ,jk) *abs( (zwdw(ji  ,jj+1,jk  )-zwdw(ji  ,jj  ,jk  )  &
                       &                             +zwdw(ji  ,jj+1,jk+1)-zwdw(ji  ,jj  ,jk+1)) &
                       &                            /(zwdw(ji  ,jj+1,jk  )+zwdw(ji  ,jj  ,jk  )  &
                       &                             -zwdw(ji  ,jj+1,jk+1)-zwdw(ji  ,jj  ,jk+1) + rsmall) )
                  zr1(4) = vmask(ji  ,jj-1,jk) *abs( (zwdw(ji  ,jj  ,jk  )-zwdw(ji  ,jj-1,jk  )  &
                       &                             +zwdw(ji  ,jj  ,jk+1)-zwdw(ji  ,jj-1,jk+1)) &
                       &                            /(zwdw(ji  ,jj  ,jk  )+zwdw(ji  ,jj-1,jk  )  &
                       &                             -zwdw(ji,  jj  ,jk+1)-zwdw(ji  ,jj-1,jk+1) + rsmall) )
                  ztmp1 = MAXVAL( zr1(1:4) )
                  zout(ji,jj,jk) = ztmp1
                  zout2d(ji,jj)  = MAX( zout2d(ji,jj), ztmp1 )
               END DO
            END DO
         END DO
         CALL lbc_lnk( zout2d(:,:), 'T', 1. )
         CALL iom_put( "stiff_tilde", zout2d(:,:) ) 
      END IF
      ! Output Horizontal Laplacian of interfaces depths (W-points):
      ! ------------------------------------------------------------
      IF ( iom_use("dh_tilde")   ) THEN
         !
         zout(:,:,1  )=0._wp
         zout(:,:,:)=0._wp
         DO jk = 2, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.                  
                  ua(ji,jj,jk) =  umask(ji,jj,jk) * e2_e1u(ji,jj) &
                                  &  * ( zwdw(ji,jj,jk) - zwdw(ji+1,jj  ,jk) )
                  va(ji,jj,jk) =  vmask(ji,jj,jk) * e1_e2v(ji,jj) & 
                                  &  * ( zwdw(ji,jj,jk) - zwdw(ji  ,jj+1,jk) )
               END DO
            END DO
         END DO
            
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztmp1 = ( (ua(ji-1,jj  ,jk) - ua(ji,jj,jk))    &
                     &  +   (va(ji  ,jj-1,jk) - va(ji,jj,jk)) ) * SQRT(r1_e1e2t(ji,jj))
                  zout(ji,jj,jk) = ABS(ztmp1)*tmask(ji,jj,jk)           
               END DO
            END DO
         END DO 
         ! Mask open boundaries:
#if defined key_bdy
         IF (lk_bdy) THEN
            DO jk = 1, jpkm1
               zout(:,:,jk) = zout(:,:,jk) * bdytmask(:,:)
            END DO
         ENDIF
#endif
         zout2d(:,:) = 0.e0 
         DO jk=1,jpkm1
            zout2d(:,:) = max( zout2d(:,:), zout(:,:,jk))
         END DO
         CALL lbc_lnk( zout2d(:,:), 'T', 1. )
         !
         CALL iom_put( "dh_tilde", zout2d(:,:) )
      ENDIF
      !
      ! Output vertical Laplacian of interfaces depths (W-points):
      ! ----------------------------------------------------------
      IF ( iom_use("dz_tilde"  ) ) THEN
         zout(:,:,1  ) = 0._wp
         zout(:,:,:) = 0._wp
         DO ji = 2, jpim1
            DO jj = 2, jpjm1
               DO jk=2,mbkt(ji,jj)-1
                  zout(ji,jj,jk) = 2._wp*ABS(tilde_e3t_n(ji,jj,jk)+e3t_0(ji,jj,jk)-tilde_e3t_n(ji,jj,jk-1)-e3t_0(ji,jj,jk-1)) & 
                                  &   /(tilde_e3t_n(ji,jj,jk)+e3t_0(ji,jj,jk)+tilde_e3t_n(ji,jj,jk-1)+e3t_0(ji,jj,jk-1)) &
                                  &   * tmask(ji,jj,jk)
               END DO 
            END DO
         END DO
         zout2d(:,:) = 0.e0 
         DO jk=1,jpkm1
            zout2d(:,:) = max( zout2d(:,:), zout(:,:,jk))
         END DO
         CALL lbc_lnk( zout2d(:,:), 'T', 1. )
         CALL iom_put( "dz_tilde", zout2d(:,:) ) 

      END IF
      !
      !
      ! Output low pass U-velocity:
      ! ---------------------------
      IF ( iom_use("un_lf_tilde"  ).AND.ln_vvl_ztilde ) THEN
         zout(:,:,jpk) = 0.e0  
         DO jk=1,jpkm1
            zout(:,:,jk) = un_lf(:,:,jk,1)/e3u_n(:,:,jk)*r1_e2u(:,:)
         END DO
         CALL iom_put( "un_lf_tilde", zout(:,:,:) )
      END IF
      !
      ! Output low pass V-velocity:
      ! ---------------------------
      IF ( iom_use("vn_lf_tilde"  ).AND.ln_vvl_ztilde ) THEN
         zout(:,:,jpk) = 0.e0  
         DO jk=1,jpkm1
            zout(:,:,jk) = vn_lf(:,:,jk,1)/e3v_n(:,:,jk)*r1_e1v(:,:)
         END DO
         CALL iom_put( "vn_lf_tilde", zout(:,:,:) )
      END IF   
      !
      ! Barotropic cell thickness anomaly:
      ! ---------------------------------- 
      IF( iom_use("e3t_star") ) THEN
         zout(:,:,:) = (e3t_n(:,:,:)-tilde_e3t_n(:,:,:)-e3t_0(:,:,:))*tmask(:,:,:) 
         CALL iom_put( "e3t_star" , zout(:,:,:) )
      ENDIF
      !
      ! Baroclinic cell thickness anomaly:
      ! ---------------------------------- 
      IF( iom_use("e3t_tilde") )  THEN
         CALL iom_put( "e3t_tilde" , tilde_e3t_n(:,:,:) )
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('dom_vvl_dia')
      !
   END SUBROUTINE dom_vvl_dia

   !!======================================================================
END MODULE domvvl
