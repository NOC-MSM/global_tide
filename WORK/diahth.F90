MODULE diahth
   !!======================================================================
   !!                       ***  MODULE  diahth  ***
   !! Ocean diagnostics: thermocline and 20 degree depth
   !!======================================================================
   !! History :  OPA  !  1994-09  (J.-P. Boulanger)  Original code
   !!                 !  1996-11  (E. Guilyardi)  OPA8 
   !!                 !  1997-08  (G. Madec)  optimization
   !!                 !  1999-07  (E. Guilyardi)  hd28 + heat content 
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2009-07  (S. Masson) hc300 bugfix + cleaning + add new diag
   !!----------------------------------------------------------------------
#if   defined key_diahth
   !!----------------------------------------------------------------------
   !!   'key_diahth' :                              thermocline depth diag.
   !!----------------------------------------------------------------------
   !!   dia_hth      : Compute varius diagnostics associated with the mixed layer
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library
   USE iom             ! I/O library
   USE timing          ! preformance summary

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_hth       ! routine called by step.F90
   PUBLIC   dia_hth_alloc ! routine called by nemogcm.F90

   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .TRUE.    !: thermocline-20d depths flag
   
   ! note: following variables should move to local variables once iom_put is always used 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hth    !: depth of the max vertical temperature gradient [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd20   !: depth of 20 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hd28   !: depth of 28 C isotherm                         [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   htc3   !: heat content of first 300 m                    [W]

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diahth.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION dia_hth_alloc()
      !!---------------------------------------------------------------------
      INTEGER :: dia_hth_alloc
      !!---------------------------------------------------------------------
      !
      ALLOCATE( hth(jpi,jpj), hd20(jpi,jpj), hd28(jpi,jpj), htc3(jpi,jpj), STAT=dia_hth_alloc )
      !
      IF( lk_mpp           )   CALL mpp_sum ( dia_hth_alloc )
      IF(dia_hth_alloc /= 0)   CALL ctl_warn('dia_hth_alloc: failed to allocate arrays.')
      !
   END FUNCTION dia_hth_alloc


   SUBROUTINE dia_hth( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hth  ***
      !!
      !! ** Purpose : Computes
      !!      the mixing layer depth (turbocline): avt = 5.e-4
      !!      the depth of strongest vertical temperature gradient
      !!      the mixed layer depth with density     criteria: rho = rho(10m or surf) + 0.03(or 0.01)
      !!      the mixed layer depth with temperature criteria: abs( tn - tn(10m) ) = 0.2       
      !!      the top of the thermochine: tn = tn(10m) - ztem2 
      !!      the pycnocline depth with density criteria equivalent to a temperature variation 
      !!                rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC) 
      !!      the barrier layer thickness
      !!      the maximal verical inversion of temperature and its depth max( 0, max of tn - tn(10m) )
      !!      the depth of the 20 degree isotherm (linear interpolation)
      !!      the depth of the 28 degree isotherm (linear interpolation)
      !!      the heat content of first 300 m
      !!
      !! ** Method : 
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER                          ::   ji, jj, jk            ! dummy loop arguments
      INTEGER                          ::   iid, ilevel           ! temporary integers
      INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ik20, ik28  ! levels
      REAL(wp)                         ::   zavt5 = 5.e-4_wp      ! Kz criterion for the turbocline depth
      REAL(wp)                         ::   zrho3 = 0.03_wp       ! density     criterion for mixed layer depth
      REAL(wp)                         ::   zrho1 = 0.01_wp       ! density     criterion for mixed layer depth
      REAL(wp)                         ::   ztem2 = 0.2_wp        ! temperature criterion for mixed layer depth
      REAL(wp)                         ::   zthick_0, zcoef       ! temporary scalars
      REAL(wp)                         ::   zztmp, zzdep          ! temporary scalars inside do loop
      REAL(wp)                         ::   zu, zv, zw, zut, zvt  ! temporary workspace
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zabs2      ! MLD: abs( tn - tn(10m) ) = ztem2 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ztm2       ! Top of thermocline: tn = tn(10m) - ztem2     
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zrho10_3   ! MLD: rho = rho10m + zrho3      
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zpycn      ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ztinv      ! max of temperature inversion
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zdepinv    ! depth of temperature inversion
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zrho0_3    ! MLD rho = rho(surf) = 0.03
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zrho0_1    ! MLD rho = rho(surf) = 0.01
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zmaxdzT    ! max of dT/dz
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zthick     ! vertical integration thickness 
      REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zdelr      ! delta rho equivalent to deltaT = 0.2
      !!----------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('dia_hth')

      IF( kt == nit000 ) THEN
         !                                      ! allocate dia_hth array
         IF( dia_hth_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dia_hth : unable to allocate standard arrays' )

         IF(.NOT. ALLOCATED(ik20) ) THEN
            ALLOCATE(ik20(jpi,jpj), ik28(jpi,jpj), &
               &      zabs2(jpi,jpj),   &
               &      ztm2(jpi,jpj),    &
               &      zrho10_3(jpi,jpj),&
               &      zpycn(jpi,jpj),   &
               &      ztinv(jpi,jpj),   &
               &      zdepinv(jpi,jpj), &
               &      zrho0_3(jpi,jpj), &
               &      zrho0_1(jpi,jpj), &
               &      zmaxdzT(jpi,jpj), &
               &      zthick(jpi,jpj),  &
               &      zdelr(jpi,jpj), STAT=ji)
            IF( lk_mpp  )   CALL mpp_sum(ji)
            IF( ji /= 0 )   CALL ctl_stop( 'STOP', 'dia_hth : unable to allocate standard ocean arrays' )
         END IF

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_hth : diagnostics of the thermocline depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
      ENDIF

      ! initialization
      ztinv  (:,:) = 0._wp  
      zdepinv(:,:) = 0._wp  
      zmaxdzT(:,:) = 0._wp  
      DO jj = 1, jpj
         DO ji = 1, jpi
            zztmp = gdepw_n(ji,jj,mbkt(ji,jj)+1) 
            hth     (ji,jj) = zztmp
            zabs2   (ji,jj) = zztmp
            ztm2    (ji,jj) = zztmp
            zrho10_3(ji,jj) = zztmp
            zpycn   (ji,jj) = zztmp
        END DO
      END DO
      IF( nla10 > 1 ) THEN 
         DO jj = 1, jpj
            DO ji = 1, jpi
               zztmp = gdepw_n(ji,jj,mbkt(ji,jj)+1) 
               zrho0_3(ji,jj) = zztmp
               zrho0_1(ji,jj) = zztmp
            END DO
         END DO
      ENDIF
      
      ! Preliminary computation
      ! computation of zdelr = (dr/dT)(T,S,10m)*(-0.2 degC)
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tmask(ji,jj,nla10) == 1. ) THEN
               zu  =  1779.50 + 11.250 * tsn(ji,jj,nla10,jp_tem) - 3.80   * tsn(ji,jj,nla10,jp_sal)                             &
                  &                                              - 0.0745 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_tem)   &
                  &                                              - 0.0100 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_sal)
               zv  =  5891.00 + 38.000 * tsn(ji,jj,nla10,jp_tem) + 3.00   * tsn(ji,jj,nla10,jp_sal)                             &
                  &                                              - 0.3750 * tsn(ji,jj,nla10,jp_tem) * tsn(ji,jj,nla10,jp_tem)
               zut =    11.25 -  0.149 * tsn(ji,jj,nla10,jp_tem) - 0.01   * tsn(ji,jj,nla10,jp_sal)
               zvt =    38.00 -  0.750 * tsn(ji,jj,nla10,jp_tem)
               zw  = (zu + 0.698*zv) * (zu + 0.698*zv)
               zdelr(ji,jj) = ztem2 * (1000.*(zut*zv - zvt*zu)/zw)
            ELSE
               zdelr(ji,jj) = 0._wp
            ENDIF
         END DO
      END DO

      ! ------------------------------------------------------------- !
      ! thermocline depth: strongest vertical gradient of temperature !
      ! turbocline depth (mixing layer depth): avt = zavt5            !
      ! MLD: rho = rho(1) + zrho3                                     !
      ! MLD: rho = rho(1) + zrho1                                     !
      ! ------------------------------------------------------------- !
      DO jk = jpkm1, 2, -1   ! loop from bottom to 2
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zzdep = gdepw_n(ji,jj,jk)
               zztmp = ( tsn(ji,jj,jk-1,jp_tem) - tsn(ji,jj,jk,jp_tem) ) / zzdep * tmask(ji,jj,jk)   ! vertical gradient of temperature (dT/dz)
               zzdep = zzdep * tmask(ji,jj,1)

               IF( zztmp > zmaxdzT(ji,jj) ) THEN                        
                  zmaxdzT(ji,jj) = zztmp   ;   hth    (ji,jj) = zzdep                ! max and depth of dT/dz
               ENDIF
               
               IF( nla10 > 1 ) THEN 
                  zztmp = rhop(ji,jj,jk) - rhop(ji,jj,1)                             ! delta rho(1)
                  IF( zztmp > zrho3 )          zrho0_3(ji,jj) = zzdep                ! > 0.03
                  IF( zztmp > zrho1 )          zrho0_1(ji,jj) = zzdep                ! > 0.01
               ENDIF

            END DO
         END DO
      END DO
      
      CALL iom_put( "mlddzt", hth )            ! depth of the thermocline
      IF( nla10 > 1 ) THEN 
         CALL iom_put( "mldr0_3", zrho0_3 )   ! MLD delta rho(surf) = 0.03
         CALL iom_put( "mldr0_1", zrho0_1 )   ! MLD delta rho(surf) = 0.01
      ENDIF

      ! ------------------------------------------------------------- !
      ! MLD: abs( tn - tn(10m) ) = ztem2                              !
      ! Top of thermocline: tn = tn(10m) - ztem2                      !
      ! MLD: rho = rho10m + zrho3                                     !
      ! pycnocline: rho = rho10m + (dr/dT)(T,S,10m)*(-0.2 degC)       !
      ! temperature inversion: max( 0, max of tn - tn(10m) )          !
      ! depth of temperature inversion                                !
      ! ------------------------------------------------------------- !
      DO jk = jpkm1, nlb10, -1   ! loop from bottom to nlb10
         DO jj = 1, jpj
            DO ji = 1, jpi
               !
               zzdep = gdepw_n(ji,jj,jk) * tmask(ji,jj,1)
               !
               zztmp = tsn(ji,jj,nla10,jp_tem) - tsn(ji,jj,jk,jp_tem)  ! - delta T(10m)
               IF( ABS(zztmp) > ztem2 )      zabs2   (ji,jj) = zzdep   ! abs > 0.2
               IF(     zztmp  > ztem2 )      ztm2    (ji,jj) = zzdep   ! > 0.2
               zztmp = -zztmp                                          ! delta T(10m)
               IF( zztmp >  ztinv(ji,jj) ) THEN                        ! temperature inversion
                  ztinv(ji,jj) = zztmp   ;   zdepinv (ji,jj) = zzdep   ! max value and depth
               ENDIF

               zztmp = rhop(ji,jj,jk) - rhop(ji,jj,nla10)              ! delta rho(10m)
               IF( zztmp > zrho3        )    zrho10_3(ji,jj) = zzdep   ! > 0.03
               IF( zztmp > zdelr(ji,jj) )    zpycn   (ji,jj) = zzdep   ! > equi. delta T(10m) - 0.2
               !
            END DO
         END DO
      END DO

      CALL iom_put( "mld_dt02", zabs2        )   ! MLD abs(delta t) - 0.2
      CALL iom_put( "topthdep", ztm2         )   ! T(10) - 0.2
      CALL iom_put( "mldr10_3", zrho10_3     )   ! MLD delta rho(10m) = 0.03
      CALL iom_put( "pycndep" , zpycn        )   ! MLD delta rho equi. delta T(10m) = 0.2
      CALL iom_put( "tinv"    , ztinv        )   ! max. temp. inv. (t10 ref) 
      CALL iom_put( "depti"   , zdepinv      )   ! depth of max. temp. inv. (t10 ref) 


      ! ----------------------------------- !
      ! search deepest level above 20C/28C  !
      ! ----------------------------------- !
      ik20(:,:) = 1
      ik28(:,:) = 1
      DO jk = 1, jpkm1   ! beware temperature is not always decreasing with depth => loop from top to bottom
         DO jj = 1, jpj
            DO ji = 1, jpi
               zztmp = tsn(ji,jj,jk,jp_tem)
               IF( zztmp >= 20. )   ik20(ji,jj) = jk
               IF( zztmp >= 28. )   ik28(ji,jj) = jk
            END DO
         END DO
      END DO

      ! --------------------------- !
      !  Depth of 20C/28C isotherm  !
      ! --------------------------- !
      DO jj = 1, jpj
         DO ji = 1, jpi
            !
            zzdep = gdepw_n(ji,jj,mbkt(ji,jj)+1)       ! depth of the oean bottom
            !
            iid = ik20(ji,jj)
            IF( iid /= 1 ) THEN 
               zztmp =      gdept_n(ji,jj,iid  )   &                     ! linear interpolation
                  &  + (    gdept_n(ji,jj,iid+1) - gdept_n(ji,jj,iid)                       )   &
                  &  * ( 20.*tmask(ji,jj,iid+1) - tsn(ji,jj,iid,jp_tem)                       )   &
                  &  / ( tsn(ji,jj,iid+1,jp_tem) - tsn(ji,jj,iid,jp_tem) + (1.-tmask(ji,jj,1)) )
               hd20(ji,jj) = MIN( zztmp , zzdep) * tmask(ji,jj,1)       ! bound by the ocean depth
            ELSE 
               hd20(ji,jj) = 0._wp
            ENDIF
            !
            iid = ik28(ji,jj)
            IF( iid /= 1 ) THEN 
               zztmp =      gdept_n(ji,jj,iid  )   &                     ! linear interpolation
                  &  + (    gdept_n(ji,jj,iid+1) - gdept_n(ji,jj,iid)                       )   &
                  &  * ( 28.*tmask(ji,jj,iid+1) -    tsn(ji,jj,iid,jp_tem)                       )   &
                  &  / (  tsn(ji,jj,iid+1,jp_tem) -    tsn(ji,jj,iid,jp_tem) + (1.-tmask(ji,jj,1)) )
               hd28(ji,jj) = MIN( zztmp , zzdep ) * tmask(ji,jj,1)      ! bound by the ocean depth
            ELSE 
               hd28(ji,jj) = 0._wp
            ENDIF

         END DO
      END DO
      CALL iom_put( "20d", hd20 )   ! depth of the 20 isotherm
      CALL iom_put( "28d", hd28 )   ! depth of the 28 isotherm

      ! ----------------------------- !
      !  Heat content of first 300 m  !
      ! ----------------------------- !

      ! find ilevel with (ilevel+1) the deepest W-level above 300m (we assume we can use e3t_1d to do this search...)
      ilevel   = 0
      zthick_0 = 0._wp
      DO jk = 1, jpkm1                      
         zthick_0 = zthick_0 + e3t_1d(jk)
         IF( zthick_0 < 300. )   ilevel = jk
      END DO
      ! surface boundary condition
      IF( ln_linssh ) THEN   ;   zthick(:,:) = sshn(:,:)   ;   htc3(:,:) = tsn(:,:,1,jp_tem) * sshn(:,:) * tmask(:,:,1)  
      ELSE                   ;   zthick(:,:) = 0._wp       ;   htc3(:,:) = 0._wp                                   
      ENDIF
      ! integration down to ilevel
      DO jk = 1, ilevel
         zthick(:,:) = zthick(:,:) + e3t_n(:,:,jk)
         htc3  (:,:) = htc3  (:,:) + e3t_n(:,:,jk) * tsn(:,:,jk,jp_tem) * tmask(:,:,jk)
      END DO
      ! deepest layer
      zthick(:,:) = 300. - zthick(:,:)   !   remaining thickness to reach 300m
      DO jj = 1, jpj
         DO ji = 1, jpi
            htc3(ji,jj) = htc3(ji,jj) + tsn(ji,jj,ilevel+1,jp_tem)                  &
               &                      * MIN( e3t_n(ji,jj,ilevel+1), zthick(ji,jj) ) * tmask(ji,jj,ilevel+1)
         END DO
      END DO
      ! from temperature to heat contain
      zcoef = rau0 * rcp
      htc3(:,:) = zcoef * htc3(:,:)
      CALL iom_put( "hc300", htc3 )      ! first 300m heat content
      !
      IF( ln_timing )   CALL timing_stop('dia_hth')
      !
   END SUBROUTINE dia_hth

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
CONTAINS
   SUBROUTINE dia_hth( kt )         ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'dia_hth: You should not have seen this print! error?', kt
   END SUBROUTINE dia_hth
#endif

   !!======================================================================
END MODULE diahth
