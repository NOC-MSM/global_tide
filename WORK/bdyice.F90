MODULE bdyice
   !!======================================================================
   !!                       ***  MODULE  bdyice  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (SI3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2012-01 (C. Rousset)  add new sea ice model 
   !!             4.0  !  2018    (C. Rousset)  SI3 compatibility 
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                          SI3 sea ice model
   !!----------------------------------------------------------------------
   !!   bdy_ice        : Application of open boundaries to ice
   !!   bdy_ice_frs    : Application of Flow Relaxation Scheme
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE ice             ! sea-ice: variables
   USE icevar          ! sea-ice: operations
   USE icecor          ! sea-ice: corrections
   USE icectl          ! sea-ice: control prints
   USE phycst          ! physical constant
   USE eosbn2          ! equation of state
   USE par_oce         ! ocean parameters
   USE dom_oce         ! ocean space and time domain variables 
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE bdy_oce         ! ocean open boundary conditions
   !
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! write to numout file
   USE lib_mpp         ! distributed memory computing
   USE lib_fortran     ! to use key_nosignedzero
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_ice     ! routine called in sbcmod
   PUBLIC   bdy_ice_dyn ! routine called in icedyn_rhg_evp

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdyice.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL licence (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_ice( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ice  ***
      !!
      !! ** Purpose : Apply open boundary conditions for sea ice
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      !
      INTEGER ::   jbdy   ! BDY set index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_ice_thd')
      !
      CALL ice_var_glo2eqv
      !
      DO jbdy = 1, nb_bdy
         !
         SELECT CASE( cn_ice(jbdy) )
         CASE('none')   ;   CYCLE
         CASE('frs' )   ;   CALL bdy_ice_frs( idx_bdy(jbdy), dta_bdy(jbdy), kt, jbdy )
         CASE DEFAULT
            CALL ctl_stop( 'bdy_ice : unrecognised option for open boundaries for ice fields' )
         END SELECT
         !
      END DO
      !
      CALL ice_cor( kt , 0 )      ! -- In case categories are out of bounds, do a remapping
      !                           !    i.e. inputs have not the same ice thickness distribution (set by rn_himean)
      !                           !         than the regional simulation
      CALL ice_var_agg(1)
      !
      IF( ln_icectl )   CALL ice_prt( kt, iiceprt, jiceprt, 1, ' - ice thermo bdy - ' )
      IF( ln_timing )   CALL timing_stop('bdy_ice_thd')
      !
   END SUBROUTINE bdy_ice


   SUBROUTINE bdy_ice_frs( idx, dta, kt, jbdy )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_frs  ***
      !!                    
      !! ** Purpose : Apply the Flow Relaxation Scheme for sea-ice fields
      !! 
      !! Reference : Engedahl H., 1995: Use of the flow relaxation scheme in a three-
      !!             dimensional baroclinic ocean model with realistic topography. Tellus, 365-382.
      !!------------------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx     ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta     ! OBC external data
      INTEGER,         INTENT(in) ::   kt      ! main time-step counter
      INTEGER,         INTENT(in) ::   jbdy    ! BDY set index
      !
      INTEGER  ::   jpbound            ! 0 = incoming ice
      !                                ! 1 = outgoing ice
      INTEGER  ::   i_bdy, jgrd        ! dummy loop indices
      INTEGER  ::   ji, jj, jk, jl, ib, jb
      REAL(wp) ::   zwgt, zwgt1        ! local scalar
      REAL(wp) ::   ztmelts, zdh
      !!------------------------------------------------------------------------------
      !
      jgrd = 1      ! Everything is at T-points here
      !
      DO jl = 1, jpl
         DO i_bdy = 1, idx%nblenrim(jgrd)
            ji    = idx%nbi(i_bdy,jgrd)
            jj    = idx%nbj(i_bdy,jgrd)
            zwgt  = idx%nbw(i_bdy,jgrd)
            zwgt1 = 1.e0 - idx%nbw(i_bdy,jgrd)
            a_i(ji,jj,jl) = ( a_i(ji,jj,jl) * zwgt1 + dta%a_i(i_bdy,jl) * zwgt ) * tmask(ji,jj,1)  ! Leads fraction 
            h_i(ji,jj,jl) = ( h_i(ji,jj,jl) * zwgt1 + dta%h_i(i_bdy,jl) * zwgt ) * tmask(ji,jj,1)  ! Ice depth 
            h_s(ji,jj,jl) = ( h_s(ji,jj,jl) * zwgt1 + dta%h_s(i_bdy,jl) * zwgt ) * tmask(ji,jj,1)  ! Snow depth

            ! -----------------
            ! Pathological case
            ! -----------------
            ! In case a) snow load would be in excess or b) ice is coming into a warmer environment that would lead to 
            ! very large transformation from snow to ice (see icethd_dh.F90)

            ! Then, a) transfer the snow excess into the ice (different from icethd_dh)
            zdh = MAX( 0._wp, ( rhos * h_s(ji,jj,jl) + ( rhoi - rau0 ) * h_i(ji,jj,jl) ) * r1_rau0 )
            ! Or, b) transfer all the snow into ice (if incoming ice is likely to melt as it comes into a warmer environment)
            !zdh = MAX( 0._wp, h_s(ji,jj,jl) * rhos / rhoi )

            ! recompute h_i, h_s
            h_i(ji,jj,jl) = MIN( hi_max(jl), h_i(ji,jj,jl) + zdh )
            h_s(ji,jj,jl) = MAX( 0._wp, h_s(ji,jj,jl) - zdh * rhoi / rhos ) 

         ENDDO
      ENDDO
      CALL lbc_bdy_lnk( a_i(:,:,:), 'T', 1., jbdy )
      CALL lbc_bdy_lnk( h_i(:,:,:), 'T', 1., jbdy )
      CALL lbc_bdy_lnk( h_s(:,:,:), 'T', 1., jbdy )

      DO jl = 1, jpl
         DO i_bdy = 1, idx%nblenrim(jgrd)
            ji = idx%nbi(i_bdy,jgrd)
            jj = idx%nbj(i_bdy,jgrd)

            ! condition on ice thickness depends on the ice velocity
            ! if velocity is outward (strictly), then ice thickness, volume... must be equal to adjacent values
            jpbound = 0   ;   ib = ji   ;   jb = jj
            !
            IF( u_ice(ji+1,jj  ) < 0. .AND. umask(ji-1,jj  ,1) == 0. )   jpbound = 1 ; ib = ji+1 ; jb = jj
            IF( u_ice(ji-1,jj  ) > 0. .AND. umask(ji+1,jj  ,1) == 0. )   jpbound = 1 ; ib = ji-1 ; jb = jj
            IF( v_ice(ji  ,jj+1) < 0. .AND. vmask(ji  ,jj-1,1) == 0. )   jpbound = 1 ; ib = ji   ; jb = jj+1
            IF( v_ice(ji  ,jj-1) > 0. .AND. vmask(ji  ,jj+1,1) == 0. )   jpbound = 1 ; ib = ji   ; jb = jj-1
            !
            IF( nn_ice_dta(jbdy) == 0 )   jpbound = 0 ; ib = ji ; jb = jj   ! case ice boundaries = initial conditions
            !                                                               !      do not make state variables dependent on velocity
            !
            IF( a_i(ib,jb,jl) > 0._wp ) THEN   ! there is ice at the boundary
               !
               a_i(ji,jj,jl) = a_i(ib,jb,jl) ! concentration
               h_i(ji,jj,jl) = h_i(ib,jb,jl) ! thickness ice
               h_s(ji,jj,jl) = h_s(ib,jb,jl) ! thickness snw
               !
               SELECT CASE( jpbound )
                  !
               CASE( 0 )   ! velocity is inward
                  !
                  oa_i(ji,jj,  jl) = rn_ice_age(jbdy) * a_i(ji,jj,jl) ! age
                  a_ip(ji,jj,  jl) = 0._wp                            ! pond concentration
                  v_ip(ji,jj,  jl) = 0._wp                            ! pond volume
                  t_su(ji,jj,  jl) = rn_ice_tem(jbdy)                 ! temperature surface
                  t_s (ji,jj,:,jl) = rn_ice_tem(jbdy)                 ! temperature snw
                  t_i (ji,jj,:,jl) = rn_ice_tem(jbdy)                 ! temperature ice
                  s_i (ji,jj,  jl) = rn_ice_sal(jbdy)                 ! salinity
                  sz_i(ji,jj,:,jl) = rn_ice_sal(jbdy)                 ! salinity profile
                  !
               CASE( 1 )   ! velocity is outward
                  !
                  oa_i(ji,jj,  jl) = oa_i(ib,jb,  jl) ! age
                  a_ip(ji,jj,  jl) = a_ip(ib,jb,  jl) ! pond concentration
                  v_ip(ji,jj,  jl) = v_ip(ib,jb,  jl) ! pond volume
                  t_su(ji,jj,  jl) = t_su(ib,jb,  jl) ! temperature surface
                  t_s (ji,jj,:,jl) = t_s (ib,jb,:,jl) ! temperature snw
                  t_i (ji,jj,:,jl) = t_i (ib,jb,:,jl) ! temperature ice
                  s_i (ji,jj,  jl) = s_i (ib,jb,  jl) ! salinity
                  sz_i(ji,jj,:,jl) = sz_i(ib,jb,:,jl) ! salinity profile
                  !
               END SELECT
               !
               IF( nn_icesal == 1 ) THEN     ! if constant salinity
                  s_i (ji,jj  ,jl) = rn_icesal
                  sz_i(ji,jj,:,jl) = rn_icesal
               ENDIF
               !
               ! global fields
               v_i (ji,jj,jl) = h_i(ji,jj,jl) * a_i(ji,jj,jl)                       ! volume ice
               v_s (ji,jj,jl) = h_s(ji,jj,jl) * a_i(ji,jj,jl)                       ! volume snw
               sv_i(ji,jj,jl) = MIN( s_i(ji,jj,jl) , sss_m(ji,jj) ) * v_i(ji,jj,jl) ! salt content
               DO jk = 1, nlay_s
                  e_s(ji,jj,jk,jl) = rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl) ) + rLfus )   ! enthalpy in J/m3
                  e_s(ji,jj,jk,jl) = e_s(ji,jj,jk,jl) * v_s(ji,jj,jl) * r1_nlay_s           ! enthalpy in J/m2
               END DO               
               DO jk = 1, nlay_i
                  ztmelts          = - rTmlt  * sz_i(ji,jj,jk,jl)             ! Melting temperature in C
                  t_i(ji,jj,jk,jl) = MIN( t_i(ji,jj,jk,jl), ztmelts + rt0 )   ! Force t_i to be lower than melting point => likely conservation issue
                  !
                  e_i(ji,jj,jk,jl) = rhoi * ( rcpi  * ( ztmelts - ( t_i(ji,jj,jk,jl) - rt0 ) )           &   ! enthalpy in J/m3
                     &                      + rLfus * ( 1._wp - ztmelts / ( t_i(ji,jj,jk,jl) - rt0 ) )   &
                     &                      - rcp   *   ztmelts )                  
                  e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i                            ! enthalpy in J/m2
               END DO
               !
            ELSE   ! no ice at the boundary
               !
               a_i (ji,jj,  jl) = 0._wp
               h_i (ji,jj,  jl) = 0._wp
               h_s (ji,jj,  jl) = 0._wp
               oa_i(ji,jj,  jl) = 0._wp
               a_ip(ji,jj,  jl) = 0._wp
               v_ip(ji,jj,  jl) = 0._wp
               t_su(ji,jj,  jl) = rt0
               t_s (ji,jj,:,jl) = rt0
               t_i (ji,jj,:,jl) = rt0 
               
               IF( nn_icesal == 1 ) THEN     ! if constant salinity
                  s_i (ji,jj  ,jl) = rn_icesal
                  sz_i(ji,jj,:,jl) = rn_icesal
               ELSE                          ! if variable salinity
                  s_i (ji,jj,jl)   = rn_simin
                  sz_i(ji,jj,:,jl) = rn_simin
               ENDIF
               !
               ! global fields
               v_i (ji,jj,  jl) = 0._wp
               v_s (ji,jj,  jl) = 0._wp
               sv_i(ji,jj,  jl) = 0._wp
               e_s (ji,jj,:,jl) = 0._wp
               e_i (ji,jj,:,jl) = 0._wp

            ENDIF
                        
         END DO
         !
      END DO ! jl

      CALL lbc_bdy_lnk( a_i (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( h_i (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( h_s (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( oa_i(:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( a_ip(:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( v_ip(:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( s_i (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( t_su(:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( v_i (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( v_s (:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( sv_i(:,:,:)  , 'T', 1., jbdy )
      CALL lbc_bdy_lnk( t_s (:,:,:,:), 'T', 1., jbdy )
      CALL lbc_bdy_lnk( e_s (:,:,:,:), 'T', 1., jbdy )
      CALL lbc_bdy_lnk( t_i (:,:,:,:), 'T', 1., jbdy )
      CALL lbc_bdy_lnk( e_i (:,:,:,:), 'T', 1., jbdy )
      !      
   END SUBROUTINE bdy_ice_frs


   SUBROUTINE bdy_ice_dyn( cd_type )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_ice_dyn  ***
      !!                    
      !! ** Purpose : Apply dynamics boundary conditions for sea-ice.
      !!
      !! ** Method :  if this adjacent grid point is not ice free, then u_ice and v_ice take its value
      !!              if                          is     ice free, then u_ice and v_ice are unchanged by BDY
      !!                                                           they keep values calculated in rheology
      !!
      !!------------------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in)  ::   cd_type   ! nature of velocity grid-points
      !
      INTEGER  ::   i_bdy, jgrd      ! dummy loop indices
      INTEGER  ::   ji, jj           ! local scalar
      INTEGER  ::   jbdy             ! BDY set index
      REAL(wp) ::   zmsk1, zmsk2, zflag
      !!------------------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('bdy_ice_dyn')
      !
      DO jbdy=1, nb_bdy
         !
         SELECT CASE( cn_ice(jbdy) )
         !
         CASE('none')
            CYCLE
            !
         CASE('frs')
            !
            IF( nn_ice_dta(jbdy) == 0 ) CYCLE            ! case ice boundaries = initial conditions 
            !                                            !      do not change ice velocity (it is only computed by rheology)
            SELECT CASE ( cd_type )
            !     
            CASE ( 'U' )  
               jgrd = 2      ! u velocity
               DO i_bdy = 1, idx_bdy(jbdy)%nblenrim(jgrd)
                  ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
                  jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
                  zflag = idx_bdy(jbdy)%flagu(i_bdy,jgrd)
                  !
                  IF ( ABS( zflag ) == 1. ) THEN  ! eastern and western boundaries
                     ! one of the two zmsk is always 0 (because of zflag)
                     zmsk1 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji+1,jj) ) ) ! 0 if no ice
                     zmsk2 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji-1,jj) ) ) ! 0 if no ice
                     !  
                     ! u_ice = u_ice of the adjacent grid point except if this grid point is ice-free (then do not change u_ice)
                     u_ice (ji,jj) = u_ice(ji+1,jj) * 0.5_wp * ABS( zflag + 1._wp ) * zmsk1 + &
                        &            u_ice(ji-1,jj) * 0.5_wp * ABS( zflag - 1._wp ) * zmsk2 + &
                        &            u_ice(ji  ,jj) * ( 1._wp - MIN( 1._wp, zmsk1 + zmsk2 ) )
                  ELSE                             ! everywhere else
                     u_ice(ji,jj) = 0._wp
                  ENDIF
                  !
               END DO
               CALL lbc_bdy_lnk( u_ice(:,:), 'U', -1., jbdy )
               !
            CASE ( 'V' )
               jgrd = 3      ! v velocity
               DO i_bdy = 1, idx_bdy(jbdy)%nblenrim(jgrd)
                  ji    = idx_bdy(jbdy)%nbi(i_bdy,jgrd)
                  jj    = idx_bdy(jbdy)%nbj(i_bdy,jgrd)
                  zflag = idx_bdy(jbdy)%flagv(i_bdy,jgrd)
                  !
                  IF ( ABS( zflag ) == 1. ) THEN  ! northern and southern boundaries
                     ! one of the two zmsk is always 0 (because of zflag)
                     zmsk1 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji,jj+1) ) ) ! 0 if no ice
                     zmsk2 = 1._wp - MAX( 0.0_wp, SIGN ( 1.0_wp , - vt_i(ji,jj-1) ) ) ! 0 if no ice
                     !  
                     ! v_ice = v_ice of the adjacent grid point except if this grid point is ice-free (then do not change v_ice)
                     v_ice (ji,jj) = v_ice(ji,jj+1) * 0.5_wp * ABS( zflag + 1._wp ) * zmsk1 + &
                        &            v_ice(ji,jj-1) * 0.5_wp * ABS( zflag - 1._wp ) * zmsk2 + &
                        &            v_ice(ji,jj  ) * ( 1._wp - MIN( 1._wp, zmsk1 + zmsk2 ) )
                  ELSE                             ! everywhere else
                     v_ice(ji,jj) = 0._wp
                  ENDIF
                  !
               END DO
               CALL lbc_bdy_lnk( v_ice(:,:), 'V', -1., jbdy )
               !
            END SELECT
            !
         CASE DEFAULT
            CALL ctl_stop( 'bdy_ice_dyn : unrecognised option for open boundaries for ice fields' )
         END SELECT
         !
      END DO
      !
      IF( ln_timing )   CALL timing_stop('bdy_ice_dyn')
      !
    END SUBROUTINE bdy_ice_dyn

#else
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_ice( kt )      ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'bdy_ice: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_ice
#endif

   !!=================================================================================
END MODULE bdyice
