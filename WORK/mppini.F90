MODULE mppini
   !!======================================================================
   !!                       ***  MODULE mppini   ***
   !! Ocean initialization : distributed memory computing initialization
   !!======================================================================
   !! History :  6.0  !  1994-11  (M. Guyon)  Original code
   !!  OPA       7.0  !  1995-04  (J. Escobar, M. Imbard)
   !!            8.0  !  1998-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
   !!  NEMO      1.0  !  2004-01  (G. Madec, J.M Molines)  F90 : free form , north fold jpni > 1
   !!            3.4  ! 2011-10  (A. C. Coward, NOCS & J. Donners, PRACE) add mpp_init_nfdcom
   !!            3.   ! 2013-06  (I. Epicoco, S. Mocavero, CMCC) mpp_init_nfdcom: setup avoiding MPI communication 
   !!            4.0  !  2016-06  (G. Madec)  use domain configuration file instead of bathymetry file
   !!            4.0  !  2017-06  (J.M. Molines, T. Lovato) merge of mppini and mppini_2
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  mpp_init          : Lay out the global domain over processors with/without land processor elimination
   !!  mpp_init_mask     : Read global bathymetric information to facilitate land suppression
   !!  mpp_init_ioipsl   : IOIPSL initialization in mpp 
   !!  mpp_init_partition: Calculate MPP domain decomposition
   !!  factorise         : Calculate the factors of the no. of MPI processes
   !!  mpp_init_nfdcom   : Setup for north fold exchanges with explicit point-to-point messaging
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE bdy_oce        ! open BounDarY  
   !
   USE lbcnfd  , ONLY : isendto, nsndto, nfsloop, nfeloop   ! Setup of north fold exchanges 
   USE lib_mpp        ! distribued memory computing library
   USE iom            ! nemo I/O library 
   USE ioipsl         ! I/O IPSL library
   USE in_out_manager ! I/O Manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC mpp_init       ! called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: mppini.F90 10126 2018-09-13 14:59:46Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if ! defined key_mpp_mpi
   !!----------------------------------------------------------------------
   !!   Default option :                            shared memory computing
   !!----------------------------------------------------------------------

   SUBROUTINE mpp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init  ***
      !!
      !! ** Purpose :   Lay out the global domain over processors.
      !!
      !! ** Method  :   Shared memory computing, set the local processor
      !!              variables to the value of the global domain
      !!----------------------------------------------------------------------
      !
      jpimax = jpiglo
      jpjmax = jpjglo
      jpi    = jpiglo
      jpj    = jpjglo
      jpk    = jpkglo
      jpim1  = jpi-1                                            ! inner domain indices
      jpjm1  = jpj-1                                            !   "           "
      jpkm1  = MAX( 1, jpk-1 )                                  !   "           "
      jpij   = jpi*jpj
      jpni   = 1
      jpnj   = 1
      jpnij  = jpni*jpnj
      nimpp  = 1           ! 
      njmpp  = 1
      nlci   = jpi
      nlcj   = jpj
      nldi   = 1
      nldj   = 1
      nlei   = jpi
      nlej   = jpj
      nbondi = 2
      nbondj = 2
      nidom  = FLIO_DOM_NONE
      npolj = jperio
      l_Iperio = jpni == 1 .AND. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7)
      l_Jperio = jpnj == 1 .AND. (jperio == 2 .OR. jperio == 7)
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'mpp_init : NO massively parallel processing'
         WRITE(numout,*) '~~~~~~~~ '
         WRITE(numout,*) '   l_Iperio = ', l_Iperio, '    l_Jperio = ', l_Jperio 
         WRITE(numout,*) '     npolj  = ',   npolj , '      njmpp  = ', njmpp
      ENDIF
      !
      IF(  jpni /= 1 .OR. jpnj /= 1 .OR. jpnij /= 1 )                                     &
         CALL ctl_stop( 'mpp_init: equality  jpni = jpnj = jpnij = 1 is not satisfied',   &
            &           'the domain is lay out for distributed memory computing!' )
         !
   END SUBROUTINE mpp_init

#else
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'                     MPI massively parallel processing
   !!----------------------------------------------------------------------


   SUBROUTINE mpp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init  ***
      !!                    
      !! ** Purpose :   Lay out the global domain over processors.
      !!      If land processors are to be eliminated, this program requires the
      !!      presence of the domain configuration file. Land processors elimination
      !!      is performed if jpni x jpnj /= jpnij. In this case, using the MPP_PREP
      !!      preprocessing tool, help for defining the best cutting out.
      !!
      !! ** Method  :   Global domain is distributed in smaller local domains.
      !!      Periodic condition is a function of the local domain position
      !!      (global boundary or neighbouring domain) and of the global
      !!      periodic
      !!      Type :         jperio global periodic condition
      !!
      !! ** Action : - set domain parameters
      !!                    nimpp     : longitudinal index 
      !!                    njmpp     : latitudinal  index
      !!                    narea     : number for local area
      !!                    nlci      : first dimension
      !!                    nlcj      : second dimension
      !!                    nbondi    : mark for "east-west local boundary"
      !!                    nbondj    : mark for "north-south local boundary"
      !!                    nproc     : number for local processor
      !!                    noea      : number for local neighboring processor
      !!                    nowe      : number for local neighboring processor
      !!                    noso      : number for local neighboring processor
      !!                    nono      : number for local neighboring processor
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jn, jproc, jarea   ! dummy loop indices
      INTEGER ::   inum                       ! local logical unit
      INTEGER ::   idir, ifreq, icont, isurf  ! local integers
      INTEGER ::   ii, il1, ili, imil         !   -       -
      INTEGER ::   ij, il2, ilj, ijm1         !   -       -
      INTEGER ::   iino, ijno, iiso, ijso     !   -       -
      INTEGER ::   iiea, ijea, iiwe, ijwe     !   -       -
      INTEGER ::   iresti, irestj, iarea0     !   -       -
      INTEGER ::   ierr                       ! local logical unit
      REAL(wp)::   zidom, zjdom               ! local scalars
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   iin, ii_nono, ii_noea          ! 1D workspace
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   ijn, ii_noso, ii_nowe          !  -     -
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::   iimppt, ilci, ibondi, ipproc   ! 2D workspace
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::   ijmppt, ilcj, ibondj, ipolj    !  -     -
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::   ilei, ildi, iono, ioea         !  -     -
      INTEGER, ALLOCATABLE, DIMENSION(:,:) ::   ilej, ildj, ioso, iowe         !  -     -
      INTEGER, DIMENSION(jpiglo,jpjglo) ::   imask   ! 2D global domain workspace
      !!----------------------------------------------------------------------

      ! If dimensions of processor grid weren't specified in the namelist file
      ! then we calculate them here now that we have our communicator size
      IF( jpni < 1 .OR. jpnj < 1 )   CALL mpp_init_partition( mppsize )
      !
#if defined key_agrif
      IF( jpnij /= jpni*jpnj ) CALL ctl_stop( 'STOP', 'Cannot remove land proc with AGRIF' )
#endif
      !
      ALLOCATE(  nfiimpp(jpni,jpnj), nfipproc(jpni,jpnj), nfilcit(jpni,jpnj) ,    &
         &       nimppt(jpnij) , ibonit(jpnij) , nlcit(jpnij) , nlcjt(jpnij) ,    &
         &       njmppt(jpnij) , ibonjt(jpnij) , nldit(jpnij) , nldjt(jpnij) ,    &
         &                                       nleit(jpnij) , nlejt(jpnij) ,    &
         &       iin(jpnij), ii_nono(jpnij), ii_noea(jpnij),   &
         &       ijn(jpnij), ii_noso(jpnij), ii_nowe(jpnij),   &
         &       iimppt(jpni,jpnj), ilci(jpni,jpnj), ibondi(jpni,jpnj), ipproc(jpni,jpnj),   &
         &       ijmppt(jpni,jpnj), ilcj(jpni,jpnj), ibondj(jpni,jpnj), ipolj(jpni,jpnj),   &
         &       ilei(jpni,jpnj), ildi(jpni,jpnj), iono(jpni,jpnj), ioea(jpni,jpnj),   &
         &       ilej(jpni,jpnj), ildj(jpni,jpnj), ioso(jpni,jpnj), iowe(jpni,jpnj),   &
         &       STAT=ierr )
      CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'mpp_init: unable to allocate standard ocean arrays' )
      
      !
#if defined key_agrif
      IF( .NOT. Agrif_Root() ) THEN       ! AGRIF children: specific setting (cf. agrif_user.F90)
         IF( jpiglo /= nbcellsx + 2 + 2*nbghostcells )   &
            CALL ctl_stop( 'STOP', 'mpp_init: Agrif children requires jpiglo == nbcellsx + 2 + 2*nbghostcells' )
         IF( jpjglo /= nbcellsy + 2 + 2*nbghostcells )   &
            CALL ctl_stop( 'STOP', 'mpp_init: Agrif children requires jpjglo == nbcellsy + 2 + 2*nbghostcells' )
         IF( ln_use_jattr )   CALL ctl_stop( 'STOP', 'mpp_init:Agrif children requires ln_use_jattr = .false. ' )
      ENDIF
#endif

#if defined key_nemocice_decomp
      jpimax = ( nx_global+2-2*nn_hls + (jpni-1) ) / jpni + 2*nn_hls    ! first  dim.
      jpjmax = ( ny_global+2-2*nn_hls + (jpnj-1) ) / jpnj + 2*nn_hls    ! second dim. 
#else
      jpimax = ( jpiglo - 2*nn_hls + (jpni-1) ) / jpni + 2*nn_hls    ! first  dim.
      jpjmax = ( jpjglo - 2*nn_hls + (jpnj-1) ) / jpnj + 2*nn_hls    ! second dim.
#endif

      !
      IF ( jpni * jpnj == jpnij ) THEN    ! regular domain lay out over processors
         imask(:,:) = 1               
      ELSEIF ( jpni*jpnj > jpnij ) THEN   ! remove land-only processor (i.e. where imask(:,:)=0)
         CALL mpp_init_mask( imask )   
      ELSE                                ! error
         CALL ctl_stop( 'mpp_init: jpnij > jpni x jpnj. Check namelist setting!' )
      ENDIF
      !
      !  1. Dimension arrays for subdomains
      ! -----------------------------------
      !  Computation of local domain sizes ilci() ilcj()
      !  These dimensions depend on global sizes jpni,jpnj and jpiglo,jpjglo
      !  The subdomains are squares lesser than or equal to the global
      !  dimensions divided by the number of processors minus the overlap array.
      !
      nreci = 2 * nn_hls
      nrecj = 2 * nn_hls
      iresti = 1 + MOD( jpiglo - nreci -1 , jpni )
      irestj = 1 + MOD( jpjglo - nrecj -1 , jpnj )
      !
      !  Need to use jpimax and jpjmax here since jpi and jpj not yet defined
#if defined key_nemocice_decomp
      ! Change padding to be consistent with CICE
      ilci(1:jpni-1      ,:) = jpimax
      ilci(jpni          ,:) = jpiglo - (jpni - 1) * (jpimax - nreci)
      !
      ilcj(:,      1:jpnj-1) = jpjmax
      ilcj(:,          jpnj) = jpjglo - (jpnj - 1) * (jpjmax - nrecj)
#else
      ilci(1:iresti      ,:) = jpimax
      ilci(iresti+1:jpni ,:) = jpimax-1

      ilcj(:,      1:irestj) = jpjmax
      ilcj(:, irestj+1:jpnj) = jpjmax-1
#endif
      !
      zidom = nreci + sum(ilci(:,1) - nreci )
      zjdom = nrecj + sum(ilcj(1,:) - nrecj )
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'mpp_init : MPI Message Passing MPI - domain lay out over processors'
         WRITE(numout,*) '~~~~~~~~ '
         WRITE(numout,*) '   defines mpp subdomains'
         WRITE(numout,*) '      iresti = ', iresti, ' jpni = ', jpni  
         WRITE(numout,*) '      irestj = ', irestj, ' jpnj = ', jpnj
         WRITE(numout,*)
         WRITE(numout,*) '      sum ilci(i,1) = ', zidom, ' jpiglo = ', jpiglo
         WRITE(numout,*) '      sum ilcj(1,j) = ', zjdom, ' jpjglo = ', jpjglo
      ENDIF

      !  2. Index arrays for subdomains
      ! -------------------------------
      iimppt(:,:) =  1
      ijmppt(:,:) =  1
      ipproc(:,:) = -1
      !
      IF( jpni > 1 ) THEN
         DO jj = 1, jpnj
            DO ji = 2, jpni
               iimppt(ji,jj) = iimppt(ji-1,jj) + ilci(ji-1,jj) - nreci
            END DO
         END DO
      ENDIF
      nfiimpp(:,:) = iimppt(:,:)
      !
      IF( jpnj > 1 )THEN
         DO jj = 2, jpnj
            DO ji = 1, jpni
               ijmppt(ji,jj) = ijmppt(ji,jj-1) + ilcj(ji,jj-1) - nrecj
            END DO
         END DO
      ENDIF

      ! 3. Subdomain description in the Regular Case
      ! --------------------------------------------
      ! specific cases where there is no communication -> must do the periodicity by itself
      ! Warning: because of potential land-area suppression, do not use nbond[ij] == 2  
      l_Iperio = jpni == 1 .AND. (jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7)
      l_Jperio = jpnj == 1 .AND. (jperio == 2 .OR. jperio == 7)
      
      icont = -1
      DO jarea = 1, jpni*jpnj
         iarea0 = jarea - 1
         ii = 1 + MOD(iarea0,jpni)
         ij = 1 +     iarea0/jpni
         ili = ilci(ii,ij)
         ilj = ilcj(ii,ij)
         ibondi(ii,ij) = 0                         ! default: has e-w neighbours
         IF( ii   ==    1 )   ibondi(ii,ij) = -1   ! first column, has only e neighbour
         IF( ii   == jpni )   ibondi(ii,ij) =  1   ! last column,  has only w neighbour
         IF( jpni ==    1 )   ibondi(ii,ij) =  2   ! has no e-w neighbour
         ibondj(ii,ij) = 0                         ! default: has n-s neighbours
         IF( ij   ==    1 )   ibondj(ii,ij) = -1   ! first row, has only n neighbour
         IF( ij   == jpnj )   ibondj(ii,ij) =  1   ! last row,  has only s neighbour
         IF( jpnj ==    1 )   ibondj(ii,ij) =  2   ! has no n-s neighbour

         ! Subdomain neighbors (get their zone number): default definition
         ioso(ii,ij) = iarea0 - jpni
         iowe(ii,ij) = iarea0 - 1
         ioea(ii,ij) = iarea0 + 1
         iono(ii,ij) = iarea0 + jpni
         ildi(ii,ij) =  1  + nn_hls
         ilei(ii,ij) = ili - nn_hls
         ildj(ii,ij) =  1  + nn_hls
         ilej(ii,ij) = ilj - nn_hls

         ! East-West periodicity: change ibondi, ioea, iowe
         IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7 ) THEN
            IF( jpni  /= 1 )   ibondi(ii,ij) = 0                        ! redefine: all have e-w neighbours
            IF( ii ==    1 )   iowe(ii,ij) = iarea0 +        (jpni-1)   ! redefine: first column, address of w neighbour
            IF( ii == jpni )   ioea(ii,ij) = iarea0 -        (jpni-1)   ! redefine: last column,  address of e neighbour
         ENDIF

         ! Simple North-South periodicity: change ibondj, ioso, iono
         IF( jperio == 2 .OR. jperio == 7 ) THEN
            IF( jpnj  /= 1 )   ibondj(ii,ij) = 0                        ! redefine: all have n-s neighbours
            IF( ij ==    1 )   ioso(ii,ij) = iarea0 + jpni * (jpnj-1)   ! redefine: first row, address of s neighbour
            IF( ij == jpnj )   iono(ii,ij) = iarea0 - jpni * (jpnj-1)   ! redefine: last row,  address of n neighbour
         ENDIF

         ! North fold: define ipolj, change iono. Warning: we do not change ibondj...
         ipolj(ii,ij) = 0
         IF( jperio == 3 .OR. jperio == 4 ) THEN
            ijm1 = jpni*(jpnj-1)
            imil = ijm1+(jpni+1)/2
            IF( jarea > ijm1 ) ipolj(ii,ij) = 3
            IF( MOD(jpni,2) == 1 .AND. jarea == imil ) ipolj(ii,ij) = 4
            IF( ipolj(ii,ij) == 3 ) iono(ii,ij) = jpni*jpnj-jarea+ijm1   ! MPI rank of northern neighbour
         ENDIF
         IF( jperio == 5 .OR. jperio == 6 ) THEN
            ijm1 = jpni*(jpnj-1)
            imil = ijm1+(jpni+1)/2
            IF( jarea > ijm1) ipolj(ii,ij) = 5
            IF( MOD(jpni,2) == 1 .AND. jarea == imil ) ipolj(ii,ij) = 6
            IF( ipolj(ii,ij) == 5) iono(ii,ij) = jpni*jpnj-jarea+ijm1    ! MPI rank of northern neighbour
         ENDIF
         !
         ! Check wet points over the entire domain to preserve the MPI communication stencil
         isurf = 0
         DO jj = 1, ilj
            DO  ji = 1, ili
               IF( imask(ji+iimppt(ii,ij)-1, jj+ijmppt(ii,ij)-1) == 1)   isurf = isurf+1
            END DO
         END DO
         !
         IF( isurf /= 0 ) THEN
            icont = icont + 1
            ipproc(ii,ij) = icont
            iin(icont+1) = ii
            ijn(icont+1) = ij
         ENDIF
      END DO
      !
      nfipproc(:,:) = ipproc(:,:)

      ! Check potential error
      IF( icont+1 /= jpnij ) THEN
         WRITE(ctmp1,*) ' jpni =',jpni,' jpnj =',jpnj
         WRITE(ctmp2,*) ' jpnij =',jpnij, '< jpni x jpnj' 
         WRITE(ctmp3,*) ' ***********, mpp_init2 finds jpnij=',icont+1
         CALL ctl_stop( 'STOP', 'mpp_init: Eliminate land processors algorithm', '', ctmp1, ctmp2, '', ctmp3 )
      ENDIF

      ! 4. Subdomain print
      ! ------------------
      IF(lwp) THEN
         ifreq = 4
         il1 = 1
         DO jn = 1, (jpni-1)/ifreq+1
            il2 = MIN(jpni,il1+ifreq-1)
            WRITE(numout,*)
            WRITE(numout,9400) ('***',ji=il1,il2-1)
            DO jj = jpnj, 1, -1
               WRITE(numout,9403) ('   ',ji=il1,il2-1)
               WRITE(numout,9402) jj, (ilci(ji,jj),ilcj(ji,jj),ji=il1,il2)
               WRITE(numout,9404) (ipproc(ji,jj),ji=il1,il2)
               WRITE(numout,9403) ('   ',ji=il1,il2-1)
               WRITE(numout,9400) ('***',ji=il1,il2-1)
            END DO
            WRITE(numout,9401) (ji,ji=il1,il2)
            il1 = il1+ifreq
         END DO
 9400    FORMAT('           ***'   ,20('*************',a3)    )
 9403    FORMAT('           *     ',20('         *   ',a3)    )
 9401    FORMAT('              '   ,20('   ',i3,'          ') )
 9402    FORMAT('       ',i3,' *  ',20(i3,'  x',i3,'   *   ') )
 9404    FORMAT('           *  '   ,20('      ',i3,'   *   ') )
      ENDIF

      ! 5. neighbour treatment: change ibondi, ibondj if next to a land zone
      ! ----------------------
      DO jarea = 1, jpni*jpnj
         ii = 1 + MOD( jarea-1  , jpni )
         ij = 1 +     (jarea-1) / jpni
         ! land-only area with an active n neigbour
         IF ( ipproc(ii,ij) == -1 .AND. 0 <= iono(ii,ij) .AND. iono(ii,ij) <= jpni*jpnj-1 ) THEN
            iino = 1 + MOD( iono(ii,ij) , jpni )                    ! ii index of this n neigbour
            ijno = 1 +      iono(ii,ij) / jpni                      ! ij index of this n neigbour
            ! In case of north fold exchange: I am the n neigbour of my n neigbour!! (#1057)
            ! --> for northern neighbours of northern row processors (in case of north-fold)
            !     need to reverse the LOGICAL direction of communication 
            idir = 1                                           ! we are indeed the s neigbour of this n neigbour
            IF( ij == jpnj .AND. ijno == jpnj )   idir = -1    ! both are on the last row, we are in fact the n neigbour
            IF( ibondj(iino,ijno) == idir     )   ibondj(iino,ijno) =   2     ! this n neigbour had only a s/n neigbour -> no more
            IF( ibondj(iino,ijno) == 0        )   ibondj(iino,ijno) = -idir   ! this n neigbour had both, n-s neighbours -> keep 1
         ENDIF
         ! land-only area with an active s neigbour
         IF( ipproc(ii,ij) == -1 .AND. 0 <= ioso(ii,ij) .AND. ioso(ii,ij) <= jpni*jpnj-1 ) THEN
            iiso = 1 + MOD( ioso(ii,ij) , jpni )                    ! ii index of this s neigbour
            ijso = 1 +      ioso(ii,ij) / jpni                      ! ij index of this s neigbour
            IF( ibondj(iiso,ijso) == -1 )   ibondj(iiso,ijso) = 2   ! this s neigbour had only a n neigbour    -> no more neigbour
            IF( ibondj(iiso,ijso) ==  0 )   ibondj(iiso,ijso) = 1   ! this s neigbour had both, n-s neighbours -> keep s neigbour
         ENDIF
         ! land-only area with an active e neigbour
         IF( ipproc(ii,ij) == -1 .AND. 0 <= ioea(ii,ij) .AND. ioea(ii,ij) <= jpni*jpnj-1 ) THEN
            iiea = 1 + MOD( ioea(ii,ij) , jpni )                    ! ii index of this e neigbour
            ijea = 1 +      ioea(ii,ij) / jpni                      ! ij index of this e neigbour
            IF( ibondi(iiea,ijea) == 1 )   ibondi(iiea,ijea) =  2   ! this e neigbour had only a w neigbour    -> no more neigbour
            IF( ibondi(iiea,ijea) == 0 )   ibondi(iiea,ijea) = -1   ! this e neigbour had both, e-w neighbours -> keep e neigbour
         ENDIF
         ! land-only area with an active w neigbour
         IF( ipproc(ii,ij) == -1 .AND. 0 <= iowe(ii,ij) .AND. iowe(ii,ij) <= jpni*jpnj-1) THEN
            iiwe = 1 + MOD( iowe(ii,ij) , jpni )                    ! ii index of this w neigbour
            ijwe = 1 +      iowe(ii,ij) / jpni                      ! ij index of this w neigbour
            IF( ibondi(iiwe,ijwe) == -1 )   ibondi(iiwe,ijwe) = 2   ! this w neigbour had only a e neigbour    -> no more neigbour
            IF( ibondi(iiwe,ijwe) ==  0 )   ibondi(iiwe,ijwe) = 1   ! this w neigbour had both, e-w neighbours -> keep w neigbour
         ENDIF
      END DO

      ! Update il[de][ij] according to modified ibond[ij]
      ! ----------------------
      DO jproc = 1, jpnij
         ii = iin(jproc)
         ij = ijn(jproc)
         IF( ibondi(ii,ij) == -1 .OR. ibondi(ii,ij) == 2 ) ildi(ii,ij) =  1
         IF( ibondi(ii,ij) ==  1 .OR. ibondi(ii,ij) == 2 ) ilei(ii,ij) = ilci(ii,ij)
         IF( ibondj(ii,ij) == -1 .OR. ibondj(ii,ij) == 2 ) ildj(ii,ij) =  1
         IF( ibondj(ii,ij) ==  1 .OR. ibondj(ii,ij) == 2 ) ilej(ii,ij) = ilcj(ii,ij)
      END DO
         
      ! just to save nono etc for all proc
      ! warning ii*ij (zone) /= nproc (processors)!
      ! ioso = zone number, ii_noso = proc number
      ii_noso(:) = -1
      ii_nono(:) = -1
      ii_noea(:) = -1
      ii_nowe(:) = -1 
      DO jproc = 1, jpnij
         ii = iin(jproc)
         ij = ijn(jproc)
         IF( 0 <= ioso(ii,ij) .AND. ioso(ii,ij) <= (jpni*jpnj-1) ) THEN
            iiso = 1 + MOD( ioso(ii,ij) , jpni )
            ijso = 1 +      ioso(ii,ij) / jpni
            ii_noso(jproc) = ipproc(iiso,ijso)
         ENDIF
         IF( 0 <= iowe(ii,ij) .AND. iowe(ii,ij) <= (jpni*jpnj-1) ) THEN
          iiwe = 1 + MOD( iowe(ii,ij) , jpni )
          ijwe = 1 +      iowe(ii,ij) / jpni
          ii_nowe(jproc) = ipproc(iiwe,ijwe)
         ENDIF
         IF( 0 <= ioea(ii,ij) .AND. ioea(ii,ij) <= (jpni*jpnj-1) ) THEN
            iiea = 1 + MOD( ioea(ii,ij) , jpni )
            ijea = 1 +      ioea(ii,ij) / jpni
            ii_noea(jproc)= ipproc(iiea,ijea)
         ENDIF
         IF( 0 <= iono(ii,ij) .AND. iono(ii,ij) <= (jpni*jpnj-1) ) THEN
            iino = 1 + MOD( iono(ii,ij) , jpni )
            ijno = 1 +      iono(ii,ij) / jpni
            ii_nono(jproc)= ipproc(iino,ijno)
         ENDIF
      END DO
    
      ! 6. Change processor name
      ! ------------------------
      ii = iin(narea)
      ij = ijn(narea)
      !
      ! set default neighbours
      noso = ii_noso(narea)
      nowe = ii_nowe(narea)
      noea = ii_noea(narea)
      nono = ii_nono(narea)
      nlci = ilci(ii,ij)  
      nldi = ildi(ii,ij)
      nlei = ilei(ii,ij)
      nlcj = ilcj(ii,ij)  
      nldj = ildj(ii,ij)
      nlej = ilej(ii,ij)
      nbondi = ibondi(ii,ij)
      nbondj = ibondj(ii,ij)
      nimpp = iimppt(ii,ij)  
      njmpp = ijmppt(ii,ij)
      jpi = nlci
      jpj = nlcj
      jpk = jpkglo                                             ! third dim
#if defined key_agrif
      ! simple trick to use same vertical grid as parent but different number of levels: 
      ! Save maximum number of levels in jpkglo, then define all vertical grids with this number.
      ! Suppress once vertical online interpolation is ok
!!$      IF(.NOT.Agrif_Root())   jpkglo = Agrif_Parent( jpkglo )
#endif
      jpim1 = jpi-1                                            ! inner domain indices
      jpjm1 = jpj-1                                            !   "           "
      jpkm1 = MAX( 1, jpk-1 )                                  !   "           "
      jpij  = jpi*jpj                                          !  jpi x j
      DO jproc = 1, jpnij
         ii = iin(jproc)
         ij = ijn(jproc)
         nlcit(jproc) = ilci(ii,ij)
         nldit(jproc) = ildi(ii,ij)
         nleit(jproc) = ilei(ii,ij)
         nlcjt(jproc) = ilcj(ii,ij)
         nldjt(jproc) = ildj(ii,ij)
         nlejt(jproc) = ilej(ii,ij)
         ibonit(jproc) = ibondi(ii,ij)
         ibonjt(jproc) = ibondj(ii,ij)
         nimppt(jproc) = iimppt(ii,ij)  
         njmppt(jproc) = ijmppt(ii,ij) 
      END DO
      nfilcit(:,:) = ilci(:,:)

      ! Save processor layout in ascii file
      IF (lwp) THEN
         CALL ctl_opn( inum, 'layout.dat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, .FALSE., narea )
         WRITE(inum,'(a)') '   jpnij   jpimax  jpjmax    jpk  jpiglo  jpjglo'//&
   &           ' ( local:    narea     jpi     jpj )'
         WRITE(inum,'(6i8,a,3i8,a)') jpnij,jpimax,jpjmax,jpk,jpiglo,jpjglo,&
   &           ' ( local: ',narea,jpi,jpj,' )'
         WRITE(inum,'(a)') 'nproc nlci nlcj nldi nldj nlei nlej nimp njmp nono noso nowe noea nbondi nbondj '

         DO jproc = 1, jpnij
            WRITE(inum,'(13i5,2i7)')   jproc-1, nlcit  (jproc), nlcjt  (jproc),   &
               &                                nldit  (jproc), nldjt  (jproc),   &
               &                                nleit  (jproc), nlejt  (jproc),   &
               &                                nimppt (jproc), njmppt (jproc),   & 
               &                                ii_nono(jproc), ii_noso(jproc),   &
               &                                ii_nowe(jproc), ii_noea(jproc),   &
               &                                ibonit (jproc), ibonjt (jproc) 
         END DO
      END IF

      !                          ! north fold parameter
      ! Defined npolj, either 0, 3 , 4 , 5 , 6
      ! In this case the important thing is that npolj /= 0
      ! Because if we go through these line it is because jpni >1 and thus
      ! we must use lbcnorthmpp, which tests only npolj =0 or npolj /= 0
      npolj = 0
      ij = ijn(narea)
      IF( jperio == 3 .OR. jperio == 4 ) THEN
         IF( ij == jpnj )   npolj = 3
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN
         IF( ij == jpnj )   npolj = 5
      ENDIF
      !
      nproc = narea-1
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '   resulting internal parameters : '
         WRITE(numout,*) '      nproc  = ', nproc
         WRITE(numout,*) '      nowe   = ', nowe  , '   noea  =  ', noea
         WRITE(numout,*) '      nono   = ', nono  , '   noso  =  ', noso
         WRITE(numout,*) '      nbondi = ', nbondi
         WRITE(numout,*) '      nbondj = ', nbondj
         WRITE(numout,*) '      npolj  = ', npolj
         WRITE(numout,*) '    l_Iperio = ', l_Iperio
         WRITE(numout,*) '    l_Jperio = ', l_Jperio
         WRITE(numout,*) '      nlci   = ', nlci
         WRITE(numout,*) '      nlcj   = ', nlcj
         WRITE(numout,*) '      nimpp  = ', nimpp
         WRITE(numout,*) '      njmpp  = ', njmpp
         WRITE(numout,*) '      nreci  = ', nreci  
         WRITE(numout,*) '      nrecj  = ', nrecj  
         WRITE(numout,*) '      nn_hls = ', nn_hls 
      ENDIF

      !                          ! Prepare mpp north fold
      IF( jperio >= 3 .AND. jperio <= 6 .AND. jpni > 1 ) THEN
         CALL mpp_ini_north
         IF (lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) '   ==>>>   North fold boundary prepared for jpni >1'
            ! additional prints in layout.dat
            WRITE(inum,*)
            WRITE(inum,*)
            WRITE(inum,*) 'number of subdomains located along the north fold : ', ndim_rank_north
            WRITE(inum,*) 'Rank of the subdomains located along the north fold : ', ndim_rank_north
            DO jproc = 1, ndim_rank_north, 5
               WRITE(inum,*) nrank_north( jproc:MINVAL( (/jproc+4,ndim_rank_north/) ) )
            END DO
         ENDIF
      ENDIF
      !
      CALL mpp_init_ioipsl       ! Prepare NetCDF output file (if necessary)
      !
      IF( ln_nnogather ) THEN
         CALL mpp_init_nfdcom     ! northfold neighbour lists
         IF (lwp) THEN
            WRITE(inum,*)
            WRITE(inum,*)
            WRITE(inum,*) 'north fold exchanges with explicit point-to-point messaging :'
            WRITE(inum,*) 'nfsloop : ', nfsloop
            WRITE(inum,*) 'nfeloop : ', nfeloop
            WRITE(inum,*) 'nsndto : ', nsndto
            WRITE(inum,*) 'isendto : ', isendto
         ENDIF
      ENDIF
      !
      IF (lwp) CLOSE(inum)   
      !
      DEALLOCATE(iin, ijn, ii_nono, ii_noea, ii_noso, ii_nowe,    &
         &       iimppt, ijmppt, ibondi, ibondj, ipproc, ipolj,   &
         &       ilci, ilcj, ilei, ilej, ildi, ildj,              &
         &       iono, ioea, ioso, iowe)
      !
    END SUBROUTINE mpp_init


    SUBROUTINE mpp_init_mask( kmask )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init_mask  ***
      !!
      !! ** Purpose : Read relevant bathymetric information in a global array
      !!              in order to provide a land/sea mask used for the elimination
      !!              of land domains, in an mpp computation.
      !!
      !! ** Method  : Read the namelist ln_zco and ln_isfcav in namelist namzgr
      !!              in order to choose the correct bathymetric information
      !!              (file and variables)  
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(jpiglo,jpjglo), INTENT(out) ::   kmask   ! global domain 
  
      INTEGER :: inum   !: logical unit for configuration file
      INTEGER :: ios    !: iostat error flag
      INTEGER ::  ijstartrow                   ! temporary integers
      REAL(wp), DIMENSION(jpiglo,jpjglo) ::   zbot, zbdy          ! global workspace
      REAL(wp) ::   zidom , zjdom          ! local scalars
      NAMELIST/nambdy/ ln_bdy, nb_bdy, ln_coords_file, cn_coords_file,           &
           &             ln_mask_file, cn_mask_file, cn_dyn2d, nn_dyn2d_dta,     &
           &             cn_dyn3d, nn_dyn3d_dta, cn_tra, nn_tra_dta,             &  
           &             ln_tra_dmp, ln_dyn3d_dmp, rn_time_dmp, rn_time_dmp_out, &
           &             cn_ice, nn_ice_dta,                                     &
           &             rn_ice_tem, rn_ice_sal, rn_ice_age,                     &
           &             ln_vol, nn_volctl, nn_rimwidth, nb_jpk_bdy
      !!----------------------------------------------------------------------
      ! 0. initialisation
      ! -----------------
      CALL iom_open( cn_domcfg, inum )
      !
      ! ocean bottom level
      CALL iom_get( inum, jpdom_unknown, 'bottom_level' , zbot , lrowattr=ln_use_jattr )  ! nb of ocean T-points
      !
      CALL iom_close( inum )
      !
      ! 2D ocean mask (=1 if at least one level of the water column is ocean, =0 otherwise)
      WHERE( zbot(:,:) > 0 )   ;   kmask(:,:) = 1
      ELSEWHERE                ;   kmask(:,:) = 0
      END WHERE
  
      ! Adjust kmask with bdy_msk if it exists
  
      REWIND( numnam_ref )              ! Namelist nambdy in reference namelist : BDY
      READ  ( numnam_ref, nambdy, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy in reference namelist (mppini)', lwp )
      !
      REWIND( numnam_cfg )              ! Namelist nambdy in configuration namelist : BDY
      READ  ( numnam_cfg, nambdy, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy in configuration namelist (mppini)', lwp )

      IF( ln_bdy .AND. ln_mask_file ) THEN
         CALL iom_open( cn_mask_file, inum )
         CALL iom_get ( inum, jpdom_unknown, 'bdy_msk', zbdy )
         CALL iom_close( inum )
         WHERE ( zbdy(:,:) <= 0. ) kmask = 0
      ENDIF
      !
   END SUBROUTINE mpp_init_mask


   SUBROUTINE mpp_init_ioipsl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE mpp_init_ioipsl  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!
      !! History :
      !!   9.0  !  04-03  (G. Madec )  MPP-IOIPSL 
      !!   " "  !  08-12  (A. Coward)  addition in case of jpni*jpnj < jpnij
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) ::   iglo, iloc, iabsf, iabsl, ihals, ihale, idid
      !!----------------------------------------------------------------------

      ! The domain is split only horizontally along i- or/and j- direction
      ! So we need at the most only 1D arrays with 2 elements.
      ! Set idompar values equivalent to the jpdom_local_noextra definition
      ! used in IOM. This works even if jpnij .ne. jpni*jpnj.
      iglo(1) = jpiglo
      iglo(2) = jpjglo
      iloc(1) = nlci
      iloc(2) = nlcj
      iabsf(1) = nimppt(narea)
      iabsf(2) = njmppt(narea)
      iabsl(:) = iabsf(:) + iloc(:) - 1
      ihals(1) = nldi - 1
      ihals(2) = nldj - 1
      ihale(1) = nlci - nlei
      ihale(2) = nlcj - nlej
      idid(1) = 1
      idid(2) = 2

      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) 'mpp_init_ioipsl :   iloc  = ', iloc (1), iloc (2)
          WRITE(numout,*) '~~~~~~~~~~~~~~~     iabsf = ', iabsf(1), iabsf(2)
          WRITE(numout,*) '                    ihals = ', ihals(1), ihals(2)
          WRITE(numout,*) '                    ihale = ', ihale(1), ihale(2)
      ENDIF
      !
      CALL flio_dom_set ( jpnij, nproc, idid, iglo, iloc, iabsf, iabsl, ihals, ihale, 'BOX', nidom)
      !
   END SUBROUTINE mpp_init_ioipsl  


   SUBROUTINE mpp_init_partition( num_pes )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE mpp_init_partition  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   num_pes   ! The number of MPI processes we have
      !
      INTEGER, PARAMETER :: nfactmax = 20
      INTEGER :: nfact ! The no. of factors returned
      INTEGER :: ierr  ! Error flag
      INTEGER :: ji
      INTEGER :: idiff, mindiff, imin ! For choosing pair of factors that are closest in value
      INTEGER, DIMENSION(nfactmax) :: ifact ! Array of factors
      !!----------------------------------------------------------------------
      !
      ierr = 0
      !
      CALL factorise( ifact, nfactmax, nfact, num_pes, ierr )
      !
      IF( nfact <= 1 ) THEN
         WRITE (numout, *) 'WARNING: factorisation of number of PEs failed'
         WRITE (numout, *) '       : using grid of ',num_pes,' x 1'
         jpnj = 1
         jpni = num_pes
      ELSE
         ! Search through factors for the pair that are closest in value
         mindiff = 1000000
         imin    = 1
         DO ji = 1, nfact-1, 2
            idiff = ABS( ifact(ji) - ifact(ji+1) )
            IF( idiff < mindiff ) THEN
               mindiff = idiff
               imin = ji
            ENDIF
         END DO
         jpnj = ifact(imin)
         jpni = ifact(imin + 1)
      ENDIF
      !
      jpnij = jpni*jpnj
      !
   END SUBROUTINE mpp_init_partition


   SUBROUTINE factorise( kfax, kmaxfax, knfax, kn, kerr )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE factorise  ***
      !!
      !! ** Purpose :   return the prime factors of n.
      !!                knfax factors are returned in array kfax which is of
      !!                maximum dimension kmaxfax.
      !! ** Method  :
      !!----------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kn, kmaxfax
      INTEGER                    , INTENT(  out) ::   kerr, knfax
      INTEGER, DIMENSION(kmaxfax), INTENT(  out) ::   kfax
      !
      INTEGER :: ifac, jl, inu
      INTEGER, PARAMETER :: ntest = 14
      INTEGER, DIMENSION(ntest) ::   ilfax
      !!----------------------------------------------------------------------
      !
      ! lfax contains the set of allowed factors.
      ilfax(:) = (/(2**jl,jl=ntest,1,-1)/)
      !
      ! Clear the error flag and initialise output vars
      kerr  = 0
      kfax  = 1
      knfax = 0
      !
      IF( kn /= 1 ) THEN      ! Find the factors of n
         !
         ! nu holds the unfactorised part of the number.
         ! knfax holds the number of factors found.
         ! l points to the allowed factor list.
         ! ifac holds the current factor.
         !
         inu   = kn
         knfax = 0
         !
         DO jl = ntest, 1, -1
            !
            ifac = ilfax(jl)
            IF( ifac > inu )   CYCLE
            !
            ! Test whether the factor will divide.
            !
            IF( MOD(inu,ifac) == 0 ) THEN
               !
               knfax = knfax + 1            ! Add the factor to the list
               IF( knfax > kmaxfax ) THEN
                  kerr = 6
                  write (*,*) 'FACTOR: insufficient space in factor array ', knfax
                  return
               ENDIF
               kfax(knfax) = ifac
               ! Store the other factor that goes with this one
               knfax = knfax + 1
               kfax(knfax) = inu / ifac
               !WRITE (*,*) 'ARPDBG, factors ',knfax-1,' & ',knfax,' are ', kfax(knfax-1),' and ',kfax(knfax)
            ENDIF
            !
         END DO
         !
      ENDIF
      !
   END SUBROUTINE factorise


   SUBROUTINE mpp_init_nfdcom
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE  mpp_init_nfdcom  ***
      !! ** Purpose :   Setup for north fold exchanges with explicit 
      !!                point-to-point messaging
      !!
      !! ** Method :   Initialization of the northern neighbours lists.
      !!----------------------------------------------------------------------
      !!    1.0  ! 2011-10  (A. C. Coward, NOCS & J. Donners, PRACE)
      !!    2.0  ! 2013-06 Setup avoiding MPI communication (I. Epicoco, S. Mocavero, CMCC) 
      !!----------------------------------------------------------------------
      INTEGER  ::   sxM, dxM, sxT, dxT, jn
      INTEGER  ::   njmppmax
      !!----------------------------------------------------------------------
      !
      njmppmax = MAXVAL( njmppt )
      !
      !initializes the north-fold communication variables
      isendto(:) = 0
      nsndto     = 0
      !
      IF ( njmpp == njmppmax ) THEN      ! if I am a process in the north
         !
         !sxM is the first point (in the global domain) needed to compute the north-fold for the current process
         sxM = jpiglo - nimppt(narea) - nlcit(narea) + 1
         !dxM is the last point (in the global domain) needed to compute the north-fold for the current process
         dxM = jpiglo - nimppt(narea) + 2
         !
         ! loop over the other north-fold processes to find the processes
         ! managing the points belonging to the sxT-dxT range
         !
         DO jn = 1, jpni
            !
            sxT = nfiimpp(jn, jpnj)                            ! sxT = 1st  point (in the global domain) of the jn process
            dxT = nfiimpp(jn, jpnj) + nfilcit(jn, jpnj) - 1    ! dxT = last point (in the global domain) of the jn process
            !
            IF    ( sxT < sxM  .AND.  sxM < dxT ) THEN
               nsndto          = nsndto + 1
               isendto(nsndto) = jn
            ELSEIF( sxM <= sxT  .AND.  dxM >= dxT ) THEN
               nsndto          = nsndto + 1
               isendto(nsndto) = jn
            ELSEIF( dxM <  dxT  .AND.  sxT <  dxM ) THEN
               nsndto          = nsndto + 1
               isendto(nsndto) = jn
            ENDIF
            !
         END DO
         nfsloop = 1
         nfeloop = nlci
         DO jn = 2,jpni-1
            IF( nfipproc(jn,jpnj) == (narea - 1) ) THEN
               IF( nfipproc(jn-1,jpnj) == -1 )   nfsloop = nldi
               IF( nfipproc(jn+1,jpnj) == -1 )   nfeloop = nlei
            ENDIF
         END DO
         !
      ENDIF
      l_north_nogather = .TRUE.
      !
   END SUBROUTINE mpp_init_nfdcom

   
#endif

   !!======================================================================
END MODULE mppini
