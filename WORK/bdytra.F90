MODULE bdytra
   !!======================================================================
   !!                       ***  MODULE  bdytra  ***
   !! Ocean tracers:   Apply boundary conditions for tracers
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.5  !  2012     (S. Mocavero, I. Epicoco) Optimization of BDY communications
   !!            4.0  !  2016     (T. Lovato) Generalize OBC structure
   !!----------------------------------------------------------------------
   !!   bdy_tra       : Apply open boundary conditions & damping to T and S
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables 
   USE bdy_oce        ! ocean open boundary conditions
   USE bdylib         ! for orlanski library routines
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   ! Local structure to rearrange tracers data
   TYPE, PUBLIC ::   ztrabdy
      REAL(wp), POINTER, DIMENSION(:,:) ::  tra
   END TYPE

   PUBLIC   bdy_tra      ! called in tranxt.F90 
   PUBLIC   bdy_tra_dmp  ! called in step.F90 

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdytra.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE bdy_tra( kt )
      !!----------------------------------------------------------------------
      !!                  ***  SUBROUTINE bdy_tra  ***
      !!
      !! ** Purpose : - Apply open boundary conditions for temperature and salinity
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! Main time step counter
      !
      INTEGER                        :: ib_bdy, jn, igrd   ! Loop indeces
      TYPE(ztrabdy), DIMENSION(jpts) :: zdta               ! Temporary data structure
      !!----------------------------------------------------------------------
      igrd = 1 

      DO ib_bdy=1, nb_bdy
         !
         zdta(1)%tra => dta_bdy(ib_bdy)%tem
         zdta(2)%tra => dta_bdy(ib_bdy)%sal
         !
         DO jn = 1, jpts
            !
            SELECT CASE( TRIM(cn_tra(ib_bdy)) )
            CASE('none'        )   ;   CYCLE
            CASE('frs'         )   ;   CALL bdy_frs ( idx_bdy(ib_bdy),                tsa(:,:,:,jn), zdta(jn)%tra )
            CASE('specified'   )   ;   CALL bdy_spe ( idx_bdy(ib_bdy),                tsa(:,:,:,jn), zdta(jn)%tra )
            CASE('neumann'     )   ;   CALL bdy_nmn ( idx_bdy(ib_bdy), igrd         , tsa(:,:,:,jn)               )
            CASE('orlanski'    )   ;   CALL bdy_orl ( idx_bdy(ib_bdy), tsb(:,:,:,jn), tsa(:,:,:,jn), zdta(jn)%tra, ll_npo=.false. )
            CASE('orlanski_npo')   ;   CALL bdy_orl ( idx_bdy(ib_bdy), tsb(:,:,:,jn), tsa(:,:,:,jn), zdta(jn)%tra, ll_npo=.true. )
            CASE('runoff'      )   ;   CALL bdy_rnf ( idx_bdy(ib_bdy),                tsa(:,:,:,jn),               jn )
            CASE DEFAULT           ;   CALL ctl_stop( 'bdy_tra : unrecognised option for open boundaries for T and S' )
            END SELECT
            ! Boundary points should be updated
            CALL lbc_bdy_lnk( tsa(:,:,:,jn), 'T', 1., ib_bdy )
            ! 
         END DO
      END DO
      !
   END SUBROUTINE bdy_tra


   SUBROUTINE bdy_rnf( idx, pta, jpa )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_rnf  ***
      !!                    
      !! ** Purpose : Specialized routine to apply TRA runoff values at OBs:
      !!                  - duplicate the neighbour value for the temperature
      !!                  - specified to 0.1 PSU for the salinity
      !! 
      !!----------------------------------------------------------------------
      TYPE(OBC_INDEX),                     INTENT(in) ::   idx  ! OBC indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pta  ! tracer trend
      INTEGER,                             INTENT(in) ::   jpa  ! TRA index
      !
      REAL(wp) ::   zwgt           ! boundary weight
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij, ip, jp ! 2D addresses
      !!----------------------------------------------------------------------
      !
      igrd = 1                       ! Everything is at T-points here
      DO ib = 1, idx%nblenrim(igrd)
         ii = idx%nbi(ib,igrd)
         ij = idx%nbj(ib,igrd)
         DO ik = 1, jpkm1
            ip = bdytmask(ii+1,ij  ) - bdytmask(ii-1,ij  )
            jp = bdytmask(ii  ,ij+1) - bdytmask(ii  ,ij-1)
            if (jpa == jp_tem) pta(ii,ij,ik) = pta(ii+ip,ij+jp,ik) * tmask(ii,ij,ik)
            if (jpa == jp_sal) pta(ii,ij,ik) =                 0.1 * tmask(ii,ij,ik)
         END DO
      END DO
      !
   END SUBROUTINE bdy_rnf


   SUBROUTINE bdy_tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                 ***  SUBROUTINE bdy_tra_dmp  ***
      !!                    
      !! ** Purpose : Apply damping for tracers at open boundaries.
      !! 
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   !
      !
      REAL(wp) ::   zwgt           ! boundary weight
      REAL(wp) ::   zta, zsa, ztime
      INTEGER  ::   ib, ik, igrd   ! dummy loop indices
      INTEGER  ::   ii, ij         ! 2D addresses
      INTEGER  ::   ib_bdy         ! Loop index
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_tra_dmp')
      !
      DO ib_bdy = 1, nb_bdy
         IF( ln_tra_dmp(ib_bdy) ) THEN
            igrd = 1                       ! Everything is at T-points here
            DO ib = 1, idx_bdy(ib_bdy)%nblen(igrd)
               ii = idx_bdy(ib_bdy)%nbi(ib,igrd)
               ij = idx_bdy(ib_bdy)%nbj(ib,igrd)
               zwgt = idx_bdy(ib_bdy)%nbd(ib,igrd)
               DO ik = 1, jpkm1
                  zta = zwgt * ( dta_bdy(ib_bdy)%tem(ib,ik) - tsb(ii,ij,ik,jp_tem) ) * tmask(ii,ij,ik)
                  zsa = zwgt * ( dta_bdy(ib_bdy)%sal(ib,ik) - tsb(ii,ij,ik,jp_sal) ) * tmask(ii,ij,ik)
                  tsa(ii,ij,ik,jp_tem) = tsa(ii,ij,ik,jp_tem) + zta
                  tsa(ii,ij,ik,jp_sal) = tsa(ii,ij,ik,jp_sal) + zsa
               END DO
            END DO
         ENDIF
      END DO
      !
      IF( ln_timing )   CALL timing_stop('bdy_tra_dmp')
      !
   END SUBROUTINE bdy_tra_dmp
 
   !!======================================================================
END MODULE bdytra
