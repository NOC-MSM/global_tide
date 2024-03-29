MODULE sshwzv   
   !!==============================================================================
   !!                       ***  MODULE  sshwzv  ***
   !! Ocean dynamics : sea surface height and vertical velocity
   !!==============================================================================
   !! History :  3.1  !  2009-02  (G. Madec, M. Leclair)  Original code
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  modified LF-RA 
   !!             -   !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea) Assimilation interface
   !!             -   !  2010-09  (D.Storkey and E.O'Dea) bug fixes for BDY module
   !!            3.3  !  2011-10  (M. Leclair) split former ssh_wzv routine and remove all vvl related work
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ssh_nxt       : after ssh
   !!   ssh_swp       : filter ans swap the ssh arrays
   !!   wzv           : compute now vertical velocity
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE sbc_oce        ! surface boundary condition: ocean
   USE domvvl         ! Variable volume
   USE divhor         ! horizontal divergence
   USE phycst         ! physical constants
   USE bdy_oce , ONLY : ln_bdy, bdytmask   ! Open BounDarY
   USE bdydyn2d       ! bdy_ssh routine
#if defined key_agrif
   USE agrif_oce_interp
#endif
   !
   USE in_out_manager ! I/O manager
   USE restart        ! only for lrst_oce
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   USE wet_dry        ! Wetting/Drying flux limiting

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ssh_nxt    ! called by step.F90
   PUBLIC   wzv        ! called by step.F90
   PUBLIC   ssh_swp    ! called by step.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sshwzv.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ssh_nxt( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ssh_nxt  ***
      !!                   
      !! ** Purpose :   compute the after ssh (ssha)
      !!
      !! ** Method  : - Using the incompressibility hypothesis, the ssh increment
      !!      is computed by integrating the horizontal divergence and multiply by
      !!      by the time step.
      !!
      !! ** action  :   ssha, after sea surface height
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step
      ! 
      INTEGER  ::   jk            ! dummy loop indice
      REAL(wp) ::   z2dt, zcoef   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zhdiv   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ssh_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_nxt : after sea surface height'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      z2dt = 2._wp * rdt                          ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
      zcoef = 0.5_wp * r1_rau0

      !                                           !------------------------------!
      !                                           !   After Sea Surface Height   !
      !                                           !------------------------------!
      IF(ln_wd_il) THEN
         CALL wad_lmt(sshb, zcoef * (emp_b(:,:) + emp(:,:)), z2dt)
      ENDIF

      CALL div_hor( kt )                               ! Horizontal divergence
      !
      zhdiv(:,:) = 0._wp
      DO jk = 1, jpkm1                                 ! Horizontal divergence of barotropic transports
        zhdiv(:,:) = zhdiv(:,:) + e3t_n(:,:,jk) * hdivn(:,:,jk)
      END DO
      !                                                ! Sea surface elevation time stepping
      ! In time-split case we need a first guess of the ssh after (using the baroclinic timestep) in order to
      ! compute the vertical velocity which can be used to compute the non-linear terms of the momentum equations.
      ! 
      ssha(:,:) = (  sshb(:,:) - z2dt * ( zcoef * ( emp_b(:,:) + emp(:,:) ) + zhdiv(:,:) )  ) * ssmask(:,:)
      !
#if defined key_agrif
      CALL agrif_ssh( kt )
#endif
      !
      IF ( .NOT.ln_dynspg_ts ) THEN
         IF( ln_bdy ) THEN
            CALL lbc_lnk( ssha, 'T', 1. )    ! Not sure that's necessary
            CALL bdy_ssh( ssha )             ! Duplicate sea level across open boundaries
         ENDIF
      ENDIF
      !                                           !------------------------------!
      !                                           !           outputs            !
      !                                           !------------------------------!
      !
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=ssha, clinfo1=' ssha  - : ', mask1=tmask )
      !
      IF( ln_timing )   CALL timing_stop('ssh_nxt')
      !
   END SUBROUTINE ssh_nxt

   
   SUBROUTINE wzv( kt, kcall )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE wzv  ***
      !!                   
      !! ** Purpose :   compute the now vertical velocity
      !!
      !! ** Method  : - Using the incompressibility hypothesis, the vertical 
      !!      velocity is computed by integrating the horizontal divergence  
      !!      from the bottom to the surface minus the scale factor evolution.
      !!        The boundary conditions are w=0 at the bottom (no flux) and.
      !!
      !! ** action  :   wn      : now vertical velocity
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step
      INTEGER, INTENT( in ), OPTIONAL     :: kcall   ! optional argument
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      LOGICAL  ::   ll_use_totflx
      REAL(wp) ::   z1_2dt       ! local scalars
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zhdiv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('wzv')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'wzv : now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~ '
         !
         wn(:,:,jpk) = 0._wp                  ! bottom boundary condition: w=0 (set once for all)
      ENDIF
      !
      ! With lagrangian coordinates, consider additionnal diffusive/upstream fluxes only for  
      ! the second call to vertical velocity computation
      ll_use_totflx = .FALSE.
      IF (( PRESENT(kcall) ).AND.(ln_vvl_ztilde .OR. ln_vvl_layer)) THEN
         IF ( kcall==2 ) ll_use_totflx=.TRUE.
      ENDIF
      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
      z1_2dt = 1. / ( 2. * rdt )                         ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 )   z1_2dt = 1. / rdt
      !
      IF (( ln_vvl_ztilde .OR. ln_vvl_layer ).AND.ll_use_totflx) THEN      ! z_tilde and layer cases
         ALLOCATE( zhdiv(jpi,jpj,jpk) ) 
         !
         DO jk = 1, jpkm1
            ! horizontal divergence of thickness diffusion transport ( velocity multiplied by e3t)
            ! - ML - note: computation already done in dom_vvl_sf_nxt. Could be optimized (not critical and clearer this way)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zhdiv(ji,jj,jk) = r1_e1e2t(ji,jj) * ( un_td(ji,jj,jk) - un_td(ji-1,jj,jk) + vn_td(ji,jj,jk) - vn_td(ji,jj-1,jk) )
               END DO
            END DO
         END DO
         CALL lbc_lnk(zhdiv, 'T', 1.)  ! - ML - Perhaps not necessary: not used for horizontal "connexions"
         !                             ! Is it problematic to have a wrong vertical velocity in boundary cells?
         !                             ! Same question holds for hdivn. Perhaps just for security
         DO jk = jpkm1, 1, -1                       ! integrate from the bottom the hor. divergence
            ! computation of w
            wn(:,:,jk) = wn(:,:,jk+1) - (  e3t_n(:,:,jk) * hdivn(:,:,jk) + zhdiv(:,:,jk)    &
               &                         + z1_2dt * ( e3t_a(:,:,jk) - e3t_b(:,:,jk) )     ) * tmask(:,:,jk)
         END DO
         !          IF( ln_vvl_layer ) wn(:,:,:) = 0.e0
         DEALLOCATE( zhdiv ) 
      ELSE   ! z_star and linear free surface cases
         DO jk = jpkm1, 1, -1                       ! integrate from the bottom the hor. divergence
            ! computation of w
            wn(:,:,jk) = wn(:,:,jk+1) - (  e3t_n(:,:,jk) * hdivn(:,:,jk)                 &
               &                         + z1_2dt * ( e3t_a(:,:,jk) - e3t_b(:,:,jk) )  ) * tmask(:,:,jk)
         END DO
      ENDIF

      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wn(:,:,jk) = wn(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !
#if defined key_agrif 
      IF( .NOT. AGRIF_Root() ) THEN 
         IF ((nbondi ==  1).OR.(nbondi == 2)) wn(nlci-1 , :     ,:) = 0.e0      ! east 
         IF ((nbondi == -1).OR.(nbondi == 2)) wn(2      , :     ,:) = 0.e0      ! west 
         IF ((nbondj ==  1).OR.(nbondj == 2)) wn(:      ,nlcj-1 ,:) = 0.e0      ! north 
         IF ((nbondj == -1).OR.(nbondj == 2)) wn(:      ,2      ,:) = 0.e0      ! south 
      ENDIF 
#endif 
      !
      IF( ln_timing )   CALL timing_stop('wzv')
      !
   END SUBROUTINE wzv


   SUBROUTINE ssh_swp( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ssh_nxt  ***
      !!
      !! ** Purpose :   achieve the sea surface  height time stepping by 
      !!              applying Asselin time filter and swapping the arrays
      !!              ssha  already computed in ssh_nxt  
      !!
      !! ** Method  : - apply Asselin time fiter to now ssh (excluding the forcing
      !!              from the filter, see Leclair and Madec 2010) and swap :
      !!                sshn = ssha + atfp * ( sshb -2 sshn + ssha )
      !!                            - atfp * rdt * ( emp_b - emp ) / rau0
      !!                sshn = ssha
      !!
      !! ** action  : - sshb, sshn   : before & now sea surface height
      !!                               ready for the next time step
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      REAL(wp) ::   zcoef   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ssh_swp')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_swp : Asselin time filter and swap of sea surface height'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !              !==  Euler time-stepping: no filter, just swap  ==!
      IF ( neuler == 0 .AND. kt == nit000 ) THEN
         sshn(:,:) = ssha(:,:)                              ! now    <-- after  (before already = now)
         !
      ELSE           !==  Leap-Frog time-stepping: Asselin filter + swap  ==!
         !                                                  ! before <-- now filtered
         sshb(:,:) = sshn(:,:) + atfp * ( sshb(:,:) - 2 * sshn(:,:) + ssha(:,:) )
         IF( .NOT.ln_linssh ) THEN                          ! before <-- with forcing removed
            zcoef = atfp * rdt * r1_rau0
            sshb(:,:) = sshb(:,:) - zcoef * (     emp_b(:,:) - emp   (:,:)   &
               &                             -    rnf_b(:,:) + rnf   (:,:)   &
               &                             + fwfisf_b(:,:) - fwfisf(:,:)   ) * ssmask(:,:)
         ENDIF
         sshn(:,:) = ssha(:,:)                              ! now <-- after
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=sshb, clinfo1=' sshb  - : ', mask1=tmask )
      !
      IF( ln_timing )   CALL timing_stop('ssh_swp')
      !
   END SUBROUTINE ssh_swp

   !!======================================================================
END MODULE sshwzv
