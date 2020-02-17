MODULE bdyvol
   !!======================================================================
   !!                       ***  MODULE  bdyvol  ***
   !! Ocean dynamic :  Volume constraint when unstructured boundary 
   !!                  and filtered free surface are used
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2006-01  (J. Chanut) Bug correction
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE bdy_oce        ! ocean open boundary conditions
   USE sbc_oce        ! ocean surface boundary conditions
   USE dom_oce        ! ocean space and time domain 
   USE phycst         ! physical constants
   USE sbcisf         ! ice shelf
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! for mppsum
   USE lib_fortran    ! Fortran routines library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_vol    ! called by ???

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdyvol.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_vol( kt )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE bdyvol  ***
      !!
      !! ** Purpose :   This routine controls the volume of the system. 
      !!      A correction velocity is calculated to correct the total transport 
      !!      through the unstructured OBC. 
      !!      The total depth used is constant (H0) to be consistent with the 
      !!      linear free surface coded in OPA 8.2    <<<=== !!gm  ???? true ????
      !!
      !! ** Method  :   The correction velocity (zubtpecor here) is defined calculating
      !!      the total transport through all open boundaries (trans_bdy) minus
      !!      the cumulate E-P flux (z_cflxemp) divided by the total lateral 
      !!      surface (bdysurftot) of the unstructured boundary. 
      !!         zubtpecor = [trans_bdy - z_cflxemp ]*(1./bdysurftot)
      !!      with z_cflxemp => sum of (Evaporation minus Precipitation)
      !!                       over all the domain in m3/s at each time step.
      !!      z_cflxemp < 0 when precipitation dominate
      !!      z_cflxemp > 0 when evaporation dominate
      !!
      !!      There are 2 options (user's desiderata): 
      !!         1/ The volume changes according to E-P, this is the default
      !!            option. In this case the cumulate E-P flux are setting to
      !!            zero (z_cflxemp=0) to calculate the correction velocity. So
      !!            it will only balance the flux through open boundaries.
      !!            (set nn_volctl to 0 in tne namelist for this option)
      !!         2/ The volume is constant even with E-P flux. In this case
      !!            the correction velocity must balance both the flux 
      !!            through open boundaries and the ones through the free
      !!            surface. 
      !!            (set nn_volctl to 1 in tne namelist for this option)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk, jb, jgrd
      INTEGER  ::   ib_bdy, ii, ij
      REAL(wp) ::   zubtpecor, z_cflxemp, ztranst
      TYPE(OBC_INDEX), POINTER :: idx
      !!-----------------------------------------------------------------------------
      !
      IF( ln_vol ) THEN
      !
      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'bdy_vol : Correction of velocities along unstructured OBC'
         IF(lwp) WRITE(numout,*)'~~~~~~~'
      END IF 

      ! Calculate the cumulate surface Flux z_cflxemp (m3/s) over all the domain
      ! -----------------------------------------------------------------------
!!gm replace these lines :
      z_cflxemp = SUM ( ( emp(:,:) - rnf(:,:) + fwfisf(:,:) ) * bdytmask(:,:) * e1e2t(:,:) ) / rau0
      IF( lk_mpp )   CALL mpp_sum( z_cflxemp )     ! sum over the global domain
!!gm   by :
!!gm      z_cflxemp = glob_sum(  ( emp(:,:)-rnf(:,:)+fwfisf(:,:) ) * bdytmask(:,:) * e1e2t(:,:)  ) / rau0
!!gm

      ! Transport through the unstructured open boundary
      ! ------------------------------------------------
      zubtpecor = 0._wp
      DO ib_bdy = 1, nb_bdy
         idx => idx_bdy(ib_bdy)
         !
         jgrd = 2                               ! cumulate u component contribution first 
         DO jb = 1, idx%nblenrim(jgrd)
            DO jk = 1, jpkm1
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               zubtpecor = zubtpecor + idx%flagu(jb,jgrd) * ua(ii,ij, jk) * e2u(ii,ij) * e3u_n(ii,ij,jk)
            END DO
         END DO
         jgrd = 3                               ! then add v component contribution
         DO jb = 1, idx%nblenrim(jgrd)
            DO jk = 1, jpkm1
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               zubtpecor = zubtpecor + idx%flagv(jb,jgrd) * va(ii,ij, jk) * e1v(ii,ij) * e3v_n(ii,ij,jk) 
            END DO
         END DO
         !
      END DO
      IF( lk_mpp )   CALL mpp_sum( zubtpecor )   ! sum over the global domain

      ! The normal velocity correction
      ! ------------------------------
      IF( nn_volctl==1 ) THEN   ;   zubtpecor = ( zubtpecor - z_cflxemp ) / bdysurftot 
      ELSE                      ;   zubtpecor =   zubtpecor               / bdysurftot
      END IF

      ! Correction of the total velocity on the unstructured boundary to respect the mass flux conservation
      ! -------------------------------------------------------------
      ztranst = 0._wp
      DO ib_bdy = 1, nb_bdy
         idx => idx_bdy(ib_bdy)
         !
         jgrd = 2                               ! correct u component
         DO jb = 1, idx%nblenrim(jgrd)
            DO jk = 1, jpkm1
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               ua(ii,ij,jk) = ua(ii,ij,jk) - idx%flagu(jb,jgrd) * zubtpecor * umask(ii,ij,jk)
               ztranst = ztranst + idx%flagu(jb,jgrd) * ua(ii,ij,jk) * e2u(ii,ij) * e3u_n(ii,ij,jk)
            END DO
         END DO
         jgrd = 3                              ! correct v component
         DO jb = 1, idx%nblenrim(jgrd)
            DO jk = 1, jpkm1
               ii = idx%nbi(jb,jgrd)
               ij = idx%nbj(jb,jgrd)
               va(ii,ij,jk) = va(ii,ij,jk) -idx%flagv(jb,jgrd) * zubtpecor * vmask(ii,ij,jk)
               ztranst = ztranst + idx%flagv(jb,jgrd) * va(ii,ij,jk) * e1v(ii,ij) * e3v_n(ii,ij,jk)
            END DO
         END DO
         !
      END DO
      IF( lk_mpp )   CALL mpp_sum( ztranst )   ! sum over the global domain
 
      ! Check the cumulated transport through unstructured OBC once barotropic velocities corrected
      ! ------------------------------------------------------
      IF( lwp .AND. MOD( kt, nwrite ) == 0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'bdy_vol : time step :', kt
         IF(lwp) WRITE(numout,*)'~~~~~~~ '
         IF(lwp) WRITE(numout,*)'          cumulate flux EMP             =', z_cflxemp  , ' (m3/s)'
         IF(lwp) WRITE(numout,*)'          total lateral surface of OBC  =', bdysurftot, '(m2)'
         IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor =', zubtpecor , '(m/s)'
         IF(lwp) WRITE(numout,*)'          cumulated transport ztranst   =', ztranst   , '(m3/s)'
      END IF 
      !
      END IF ! ln_vol
      !
   END SUBROUTINE bdy_vol

   !!======================================================================
END MODULE bdyvol
