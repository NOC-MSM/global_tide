MODULE bdydyn2d
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to barotropic solution
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!----------------------------------------------------------------------
   !!   bdy_dyn2d          : Apply open boundary conditions to barotropic variables.
   !!   bdy_dyn2d_frs      : Apply Flow Relaxation Scheme 
   !!   bdy_dyn2d_fla      : Apply Flather condition
   !!   bdy_dyn2d_orlanski : Orlanski Radiation
   !!   bdy_ssh            : Duplicate sea level across open boundaries
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE bdy_oce         ! ocean open boundary conditions
   USE bdylib          ! BDY library routines
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wet_dry         ! Use wet dry to get reference ssh level
   USE in_out_manager  !

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dyn2d   ! routine called in dynspg_ts and bdy_dyn
   PUBLIC   bdy_ssh       ! routine called in dynspg_ts or sshwzv

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdydyn2d.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_dyn2d( kt, pua2d, pva2d, pub2d, pvb2d, phur, phvr, pssh  )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn2d  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for barotropic variables
      !!
      !!----------------------------------------------------------------------
      INTEGER,                      INTENT(in) ::   kt   ! Main time step counter
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d 
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pub2d, pvb2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: phur, phvr
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pssh
      !!
      INTEGER                                  ::   ib_bdy ! Loop counter

      DO ib_bdy=1, nb_bdy

         SELECT CASE( cn_dyn2d(ib_bdy) )
         CASE('none')
            CYCLE
         CASE('frs')
            CALL bdy_dyn2d_frs( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d )
         CASE('flather')
            CALL bdy_dyn2d_fla( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, pua2d, pva2d, pssh, phur, phvr )
         CASE('orlanski')
            CALL bdy_dyn2d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, &
                                     & pua2d, pva2d, pub2d, pvb2d, ll_npo=.false.)
         CASE('orlanski_npo')
            CALL bdy_dyn2d_orlanski( idx_bdy(ib_bdy), dta_bdy(ib_bdy), ib_bdy, &
                                     & pua2d, pva2d, pub2d, pvb2d, ll_npo=.true. )
         CASE DEFAULT
            CALL ctl_stop( 'bdy_dyn2d : unrecognised option for open boundaries for barotropic variables' )
         END SELECT
      ENDDO

   END SUBROUTINE bdy_dyn2d

   SUBROUTINE bdy_dyn2d_frs( idx, dta, ib_bdy, pua2d, pva2d )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_dyn2d_frs  ***
      !!
      !! ** Purpose : - Apply the Flow Relaxation Scheme for barotropic velocities
      !!                at open boundaries.
      !!
      !! References :- Engedahl H., 1995: Use of the flow relaxation scheme in 
      !!               a three-dimensional baroclinic ocean model with realistic
      !!               topography. Tellus, 365-382.
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX), INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),  INTENT(in) ::   dta  ! OBC external data
      INTEGER,         INTENT(in) ::   ib_bdy  ! BDY set index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d 
      !!
      INTEGER  ::   jb, jk         ! dummy loop indices
      INTEGER  ::   ii, ij, igrd   ! local integers
      REAL(wp) ::   zwgt           ! boundary weight
      !!----------------------------------------------------------------------
      !
      igrd = 2                      ! Relaxation of zonal velocity
      DO jb = 1, idx%nblen(igrd)
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
         zwgt = idx%nbw(jb,igrd)
         pua2d(ii,ij) = ( pua2d(ii,ij) + zwgt * ( dta%u2d(jb) - pua2d(ii,ij) ) ) * umask(ii,ij,1)
      END DO
      !
      igrd = 3                      ! Relaxation of meridional velocity
      DO jb = 1, idx%nblen(igrd)
         ii   = idx%nbi(jb,igrd)
         ij   = idx%nbj(jb,igrd)
         zwgt = idx%nbw(jb,igrd)
         pva2d(ii,ij) = ( pva2d(ii,ij) + zwgt * ( dta%v2d(jb) - pva2d(ii,ij) ) ) * vmask(ii,ij,1)
      END DO 
      CALL lbc_bdy_lnk( pua2d, 'U', -1., ib_bdy ) 
      CALL lbc_bdy_lnk( pva2d, 'V', -1., ib_bdy)   ! Boundary points should be updated
      !
   END SUBROUTINE bdy_dyn2d_frs


   SUBROUTINE bdy_dyn2d_fla( idx, dta, ib_bdy, pua2d, pva2d, pssh, phur, phvr )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn2d_fla  ***
      !!             
      !!              - Apply Flather boundary conditions on normal barotropic velocities 
      !!
      !! ** WARNINGS about FLATHER implementation:
      !!1. According to Palma and Matano, 1998 "after ssh" is used. 
      !!   In ROMS and POM implementations, it is "now ssh". In the current 
      !!   implementation (tested only in the EEL-R5 conf.), both cases were unstable. 
      !!   So I use "before ssh" in the following.
      !!
      !!2. We assume that the normal ssh gradient at the bdy is zero. As a matter of 
      !!   fact, the model ssh just inside the dynamical boundary is used (the outside  
      !!   ssh in the code is not updated).
      !!
      !! References:  Flather, R. A., 1976: A tidal model of the northwest European
      !!              continental shelf. Mem. Soc. R. Sci. Liege, Ser. 6,10, 141-164.     
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! BDY set index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pssh, phur, phvr 

      INTEGER  ::   jb, igrd                         ! dummy loop indices
      INTEGER  ::   ii, ij, iim1, iip1, ijm1, ijp1   ! 2D addresses
      REAL(wp), POINTER :: flagu, flagv              ! short cuts
      REAL(wp) ::   zcorr                            ! Flather correction
      REAL(wp) ::   zforc                            ! temporary scalar
      REAL(wp) ::   zflag, z1_2                      !    "        "
      !!----------------------------------------------------------------------

      z1_2 = 0.5_wp

      ! ---------------------------------!
      ! Flather boundary conditions     :!
      ! ---------------------------------! 
     
!!! REPLACE spgu with nemo_wrk work space

      ! Fill temporary array with ssh data (here spgu):
      igrd = 1
      spgu(:,:) = 0.0
      DO jb = 1, idx%nblenrim(igrd)
         ii = idx%nbi(jb,igrd)
         ij = idx%nbj(jb,igrd)
         IF( ll_wd ) THEN
            spgu(ii, ij) = dta%ssh(jb)  - ssh_ref 
         ELSE
            spgu(ii, ij) = dta%ssh(jb)
         ENDIF
      END DO

      CALL lbc_bdy_lnk( spgu(:,:), 'T', 1., ib_bdy )
      !
      igrd = 2      ! Flather bc on u-velocity; 
      !             ! remember that flagu=-1 if normal velocity direction is outward
      !             ! I think we should rather use after ssh ?
      DO jb = 1, idx%nblenrim(igrd)
         ii  = idx%nbi(jb,igrd)
         ij  = idx%nbj(jb,igrd) 
         flagu => idx%flagu(jb,igrd)
         iim1 = ii + MAX( 0, INT( flagu ) )   ! T pts i-indice inside the boundary
         iip1 = ii - MIN( 0, INT( flagu ) )   ! T pts i-indice outside the boundary 
         !
         zcorr = - flagu * SQRT( grav * phur(ii, ij) ) * ( pssh(iim1, ij) - spgu(iip1,ij) )

         ! jchanut tschanges: Set zflag to 0 below to revert to Flather scheme
         ! Use characteristics method instead
         zflag = ABS(flagu)
         zforc = dta%u2d(jb) * (1._wp - z1_2*zflag) + z1_2 * zflag * pua2d(iim1,ij)
         pua2d(ii,ij) = zforc + (1._wp - z1_2*zflag) * zcorr * umask(ii,ij,1) 
      END DO
      !
      igrd = 3      ! Flather bc on v-velocity
      !             ! remember that flagv=-1 if normal velocity direction is outward
      DO jb = 1, idx%nblenrim(igrd)
         ii  = idx%nbi(jb,igrd)
         ij  = idx%nbj(jb,igrd) 
         flagv => idx%flagv(jb,igrd)
         ijm1 = ij + MAX( 0, INT( flagv ) )   ! T pts j-indice inside the boundary
         ijp1 = ij - MIN( 0, INT( flagv ) )   ! T pts j-indice outside the boundary 
         !
         zcorr = - flagv * SQRT( grav * phvr(ii, ij) ) * ( pssh(ii, ijm1) - spgu(ii,ijp1) )

         ! jchanut tschanges: Set zflag to 0 below to revert to std Flather scheme
         ! Use characteristics method instead
         zflag = ABS(flagv)
         zforc  = dta%v2d(jb) * (1._wp - z1_2*zflag) + z1_2 * zflag * pva2d(ii,ijm1)
         pva2d(ii,ij) = zforc + (1._wp - z1_2*zflag) * zcorr * vmask(ii,ij,1)
      END DO
      CALL lbc_bdy_lnk( pua2d, 'U', -1., ib_bdy )   ! Boundary points should be updated
      CALL lbc_bdy_lnk( pva2d, 'V', -1., ib_bdy )   !
      !
   END SUBROUTINE bdy_dyn2d_fla


   SUBROUTINE bdy_dyn2d_orlanski( idx, dta, ib_bdy, pua2d, pva2d, pub2d, pvb2d, ll_npo )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_dyn2d_orlanski  ***
      !!             
      !!              - Apply Orlanski radiation condition adaptively:
      !!                  - radiation plus weak nudging at outflow points
      !!                  - no radiation and strong nudging at inflow points
      !! 
      !!
      !! References:  Marchesiello, McWilliams and Shchepetkin, Ocean Modelling vol. 3 (2001)    
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),              INTENT(in) ::   idx  ! OBC indices
      TYPE(OBC_DATA),               INTENT(in) ::   dta  ! OBC external data
      INTEGER,                      INTENT(in) ::   ib_bdy  ! number of current open boundary set
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pua2d, pva2d
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pub2d, pvb2d 
      LOGICAL,                      INTENT(in) ::   ll_npo  ! flag for NPO version

      INTEGER  ::   ib, igrd                               ! dummy loop indices
      INTEGER  ::   ii, ij, iibm1, ijbm1                   ! indices
      !!----------------------------------------------------------------------
      !
      igrd = 2      ! Orlanski bc on u-velocity; 
      !            
      CALL bdy_orlanski_2d( idx, igrd, pub2d, pua2d, dta%u2d, ll_npo )

      igrd = 3      ! Orlanski bc on v-velocity
      !  
      CALL bdy_orlanski_2d( idx, igrd, pvb2d, pva2d, dta%v2d, ll_npo )
      !
      CALL lbc_bdy_lnk( pua2d, 'U', -1., ib_bdy )   ! Boundary points should be updated
      CALL lbc_bdy_lnk( pva2d, 'V', -1., ib_bdy )   !
      !
   END SUBROUTINE bdy_dyn2d_orlanski


   SUBROUTINE bdy_ssh( zssh )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_ssh  ***
      !!
      !! ** Purpose : Duplicate sea level across open boundaries
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   zssh ! Sea level
      !!
      INTEGER  ::   ib_bdy, ib, igrd                        ! local integers
      INTEGER  ::   ii, ij, zcoef, zcoef1, zcoef2, ip, jp   !   "       "

      igrd = 1                       ! Everything is at T-points here

      DO ib_bdy = 1, nb_bdy
         DO ib = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
            ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
            ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
            ! Set gradient direction:
            zcoef1 = bdytmask(ii-1,ij  ) +  bdytmask(ii+1,ij  )
            zcoef2 = bdytmask(ii  ,ij-1) +  bdytmask(ii  ,ij+1)
            IF ( zcoef1+zcoef2 == 0 ) THEN
               ! corner
!               zcoef = tmask(ii-1,ij,1) + tmask(ii+1,ij,1) +  tmask(ii,ij-1,1) +  tmask(ii,ij+1,1)
!               zssh(ii,ij) = zssh(ii-1,ij  ) * tmask(ii-1,ij  ,1) + &
!                 &           zssh(ii+1,ij  ) * tmask(ii+1,ij  ,1) + &
!                 &           zssh(ii  ,ij-1) * tmask(ii  ,ij-1,1) + &
!                 &           zssh(ii  ,ij+1) * tmask(ii  ,ij+1,1)
               zcoef = bdytmask(ii-1,ij) + bdytmask(ii+1,ij) +  bdytmask(ii,ij-1) +  bdytmask(ii,ij+1)
               zssh(ii,ij) = zssh(ii-1,ij  ) * bdytmask(ii-1,ij  ) + &
                 &           zssh(ii+1,ij  ) * bdytmask(ii+1,ij  ) + &
                 &           zssh(ii  ,ij-1) * bdytmask(ii  ,ij-1) + &
                 &           zssh(ii  ,ij+1) * bdytmask(ii  ,ij+1)
               zssh(ii,ij) = ( zssh(ii,ij) / MAX( 1, zcoef) ) * tmask(ii,ij,1)
            ELSE
               ip = bdytmask(ii+1,ij  ) - bdytmask(ii-1,ij  )
               jp = bdytmask(ii  ,ij+1) - bdytmask(ii  ,ij-1)
               zssh(ii,ij) = zssh(ii+ip,ij+jp) * tmask(ii+ip,ij+jp,1)
            ENDIF
         END DO

         ! Boundary points should be updated
         CALL lbc_bdy_lnk( zssh(:,:), 'T', 1., ib_bdy )
      END DO

   END SUBROUTINE bdy_ssh

   !!======================================================================
END MODULE bdydyn2d

