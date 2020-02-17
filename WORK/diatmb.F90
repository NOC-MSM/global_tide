MODULE diatmb 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.6  !  08-2014  (E O'Dea)  Original code
   !!            3.7  !  05-2016  (G. Madec)  use mbkt, mikt to account for ocean cavities
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   !
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library

   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_diatmb     !: Top Middle and Bottom output 
   PUBLIC   dia_tmb_init            ! routine called by nemogcm.F90
   PUBLIC   dia_tmb                 ! routine called by diawri.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: diatmb.F90 10126 2018-09-13 14:59:46Z jchanut $
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_tmb_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb_init  ***
      !!     
      !! ** Purpose :   Initialization of tmb namelist 
      !!        
      !! ** Method  :   Read namelist
      !!---------------------------------------------------------------------------
      INTEGER ::   ios                 ! Local integer output status for namelist read
      !
      NAMELIST/nam_diatmb/ ln_diatmb
      !!----------------------------------------------------------------------
      !
      REWIND( numnam_ref )              ! Read Namelist nam_diatmb in reference namelist : TMB diagnostics
      READ  ( numnam_ref, nam_diatmb, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diatmb in reference namelist', lwp )
 
      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_diatmb, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nam_diatmb in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diatmb )

      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dia_tmb_init : Output Top, Middle, Bottom Diagnostics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_diatmb : set tmb outputs '
         WRITE(numout,*) '      Switch for TMB diagnostics (T) or not (F)  ln_diatmb  = ', ln_diatmb
      ENDIF
      !
   END SUBROUTINE dia_tmb_init


   SUBROUTINE dia_calctmb( pfield, ptmb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb  ***
      !!                   
      !! ** Purpose :    Find the Top, Mid and Bottom fields of water Column
      !!
      !! ** Method  :    use mbkt, mikt to find surface, mid and bottom of 
      !!              model levels due to potential existence of ocean cavities
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi, jpj, jpk), INTENT(in   ) ::   pfield   ! Input 3D field and mask
      REAL(wp), DIMENSION(jpi, jpj,  3 ), INTENT(  out) ::   ptmb     ! top, middle, bottom extracted from pfield
      !
      INTEGER ::   ji, jj   ! Dummy loop indices
      INTEGER ::   itop, imid, ibot   ! local integers
      REAL(wp)::   zmdi = 1.e+20_wp   ! land value
      !!---------------------------------------------------------------------
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            itop = mikt(ji,jj)                        ! top    ocean 
            ibot = mbkt(ji,jj)                        ! bottom ocean 
            imid =  itop + ( ibot - itop + 1 ) / 2    ! middle ocean          
            !                    
            ptmb(ji,jj,1) = pfield(ji,jj,itop)*tmask(ji,jj,itop) + zmdi*( 1._wp-tmask(ji,jj,itop) )
            ptmb(ji,jj,2) = pfield(ji,jj,imid)*tmask(ji,jj,imid) + zmdi*( 1._wp-tmask(ji,jj,imid) )
            ptmb(ji,jj,3) = pfield(ji,jj,ibot)*tmask(ji,jj,ibot) + zmdi*( 1._wp-tmask(ji,jj,ibot) )
         END DO
      END DO
      !
   END SUBROUTINE dia_calctmb


   SUBROUTINE dia_tmb
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_tmb  ***
      !! ** Purpose :   Write diagnostics for Top, Mid and Bottom of water Column
      !!
      !! ** Method  :  use mikt,mbkt to find surface, mid and bottom of model levels
      !!      calls calctmb to retrieve TMB values before sending to iom_put
      !!
      !!--------------------------------------------------------------------
      REAL(wp) ::   zmdi =1.e+20     ! land value
      REAL(wp), DIMENSION(jpi,jpj,3) :: zwtmb    ! workspace 
      !!--------------------------------------------------------------------
      !
      CALL dia_calctmb( tsn(:,:,:,jp_tem), zwtmb )
      !ssh already output but here we output it masked
      CALL iom_put( "sshnmasked", sshn(:,:)*tmask(:,:,1) + zmdi*(1.0 - tmask(:,:,1)) )
      CALL iom_put( "top_temp"  , zwtmb(:,:,1) )    ! tmb Temperature
      CALL iom_put( "mid_temp"  , zwtmb(:,:,2) )    ! tmb Temperature
      CALL iom_put( "bot_temp"  , zwtmb(:,:,3) )    ! tmb Temperature
      !
      CALL dia_calctmb( tsn(:,:,:,jp_sal), zwtmb )
      CALL iom_put( "top_sal"   , zwtmb(:,:,1) )    ! tmb Salinity 
      CALL iom_put( "mid_sal"   , zwtmb(:,:,2) )    ! tmb Salinity
      CALL iom_put( "bot_sal"   , zwtmb(:,:,3) )    ! tmb Salinity
      !
      CALL dia_calctmb( un(:,:,:), zwtmb )
      CALL iom_put( "top_u"     , zwtmb(:,:,1) )    ! tmb  U Velocity
      CALL iom_put( "mid_u"     , zwtmb(:,:,2) )    ! tmb  U Velocity
      CALL iom_put( "bot_u"     , zwtmb(:,:,3) )    ! tmb  U Velocity
      !
      CALL dia_calctmb( vn(:,:,:), zwtmb )
      CALL iom_put( "top_v"     , zwtmb(:,:,1) )    ! tmb  V Velocity
      CALL iom_put( "mid_v"     , zwtmb(:,:,2) )    ! tmb  V Velocity
      CALL iom_put( "bot_v"     , zwtmb(:,:,3) )    ! tmb  V Velocity
      !
   END SUBROUTINE dia_tmb

   !!======================================================================
END MODULE diatmb
