MODULE trcext_blingv0

#if defined key_bling

   !!======================================================================
   !!                     ***  MODULE trcext_bling  ***
   !! TOP :  Iron external inputs
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014
   !!======================================================================

   USE oce_trc
   USE par_trc
   USE trc
   USE fldread 
   USE iom

   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_ext_bling
   PUBLIC trc_ext_init_bling

   INTEGER :: numphpext, numfedext, numoxyext

   TYPE(FLD), ALLOCATABLE, DIMENSION(:)     :: sf_dust_bling

# include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_ext_bling (kt) 

      INTEGER, INTENT( IN ) ::   kt   ! ocean time step

      INTEGER  :: ji, jj, jk, ikb

      REAL(wp) :: sss, sst, o2fact, kfact
      REAL(wp) :: tt, tk, ts, ts2, ts3, ts4, ts5
      REAL(wp) :: sat_o2, sch_o2, sch_no_term, alpha_o2, csurf_o2

      REAL(wp) :: zrfact, f_oxy, fe_2_p_sed
      REAL(wp) :: sumbfpo4_glob, sumtffed_glob, sumbffed_glob, sumtfoxy_glob, sumbfoxy_glob

      REAL(wp), POINTER, DIMENSION(:,:) :: bfpo4, bffed, bfoxy, o2flx

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_ext_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      CALL wrk_alloc( jpi, jpj, bfpo4, bffed, bfoxy, o2flx )

      !---------------------------------------------------------------------
      ! Calculate air-sea exhange for O2
      !---------------------------------------------------------------------

      kfact  = 0.01_wp / 3600._wp ! convert from cm/h to m/s (piston velocity)

      DO jj=1,jpj
        DO ji=1,jpi

          zrfact = rfact / fse3t(ji,jj,1)
          o2fact = 1027.d0 /22391.6d0     ! convert from ml/l to mol/m3 w/ density 1027kg/m3

          sst=MIN( 35.,tsn(ji,jj,1,jp_tem) )
          sss=tsn(ji,jj,1,jp_sal)

          tt=298.15d0-sst
          tk=273.15d0+sst
          ts=LOG(tt/tk)
          ts2=ts *ts
          ts3=ts2*ts
          ts4=ts3*ts
          ts5=ts4*ts

          sat_o2=o2fact*EXP(        a_0+a_1*ts+a_2*ts2+a_3*ts3+a_4*ts4+a_5*ts5 &
                          + sss*( b_0+b_1*ts+b_2*ts2+b_3*ts3 + c_0*sss)       )          

          sch_o2=a_1_o2+sst*(a_2_o2+sst*(a_3_o2+sst*a_4_o2))

          ! air-sea gas transfer velocity (piston velocity). Units are cm/h
          sch_no_term=0.39 * (1.-fr_i(ji,jj)) * wndm(ji,jj)**2 * SQRT(660.d0/sch_o2)

          ! units are (m/s)*(mol/m3)
          alpha_o2=sat_o2*sch_no_term*kfact
          csurf_o2=trn(ji,jj,1,jpOxy_bling)*sch_no_term*kfact

          ! Flux of oxygen (mol/m2/s)
          o2flx(ji,jj)   =alpha_o2-csurf_o2

          trn(ji,jj,1,jpOxy_bling)=trn(ji,jj,1,jpOxy_bling)+o2flx(ji,jj)*zrfact
          !tra(ji,jj,1,jpOxy_bling)=tra(ji,jj,1,jpOxy_bling)+o2flx(ji,jj)/fse3t(ji,jj,1)

          !dum1(ji,jj,11)=dum1(ji,jj,11)+o2flx*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)
          dum1(ji,jj,11)=o2flx(ji,jj)*zrfact

          sumtfoxy=sumtfoxy+o2flx(ji,jj)*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)

        ENDDO
      ENDDO

      !---------------------------------------------------------------------
      ! Calculate external fluxes for iron. 
      !---------------------------------------------------------------------

      !Get dust field
      IF( ln_dust_bling ) THEN
         CALL fld_read( kt, 1, sf_dust_bling )
         dust_bling(:,:) = sf_dust_bling(1)%fnow(:,:,1)
      ENDIF

      fe_2_p_sed=106.0e-4

      DO jj = 1, jpj
         DO ji = 1, jpi
            zrfact  = rfact / fse3t(ji,jj,1)
            trn(ji,jj,1,jpFed_bling) =  trn(ji,jj,1,jpFed_bling) + dust_bling(ji,jj)*zrfact
            !tra(ji,jj,1,jpFed_bling) = tra(ji,jj,1,jpFed_bling) + dust_bling(ji,jj) / fse3t(ji,jj,1)

            dum1(ji,jj,6)=dust_bling(ji,jj)*zrfact

            sumtffed=sumtffed+dust_bling(ji,jj)*zrfact*cvol(ji,jj,1)*tmask(ji,jj,1)

         ENDDO
      ENDDO


            ! Exchange with sediments
            ! -------------------------------------
      DO jj = 1, jpj
        DO ji = 1, jpi


            ! mbkt is a matrix containing the vertical index of the
            ! bottom layer at each horizontal point
            ikb     = mbkt(ji,jj)

            f_oxy   = trn(ji,jj,ikb,jpOxy_bling)
            zrfact  = rfact / fse3t(ji,jj,ikb)

            ! Phosphate
            bfpo4(ji,jj)=fpop_b(ji,jj)

            ! Oxygen
            bfoxy(ji,jj) = -oxy2p*fpop_b(ji,jj)

            ! Iron
            IF (f_oxy>oxy_min) THEN
               bffed(ji,jj)= fe_2_p_sed*fpop_b(ji,jj)
            ELSE
               bffed(ji,jj)=(fe_2_p_sed*fpop_b(ji,jj)+fpofe_b(ji,jj))
            ENDIF

            ! Add the bottom flux trend
            trn(ji,jj,ikb,jpPO4_bling) = trn(ji,jj,ikb,jpPO4_bling) + bfpo4(ji,jj)*zrfact
            trn(ji,jj,ikb,jpFed_bling) = trn(ji,jj,ikb,jpFed_bling) + bffed(ji,jj)*zrfact
            trn(ji,jj,ikb,jpOxy_bling) = trn(ji,jj,ikb,jpOxy_bling) + bfoxy(ji,jj)*zrfact

            dum1(ji,jj, 3)=bfpo4  (ji,jj)*zrfact
            dum1(ji,jj, 7)=bffed  (ji,jj)*zrfact                                                 
            dum1(ji,jj, 8)=fpofe_b(ji,jj)*zrfact
            dum1(ji,jj,12)=bfoxy  (ji,jj)*zrfact

            bfpo4(ji,jj)=bfpo4(ji,jj)*tmask(ji,jj,ikb)
            bffed(ji,jj)=bffed(ji,jj)*tmask(ji,jj,ikb)
            bfoxy(ji,jj)=bfoxy(ji,jj)*tmask(ji,jj,ikb)

         END DO
      END DO

      IF (ln_bling_ext) THEN

        !IF (narea==1) THEN

          IF( kt == nittrc000 ) THEN
            CALL ctl_opn( numphpext, 'php.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numfedext, 'fed.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
            CALL ctl_opn( numoxyext, 'oxy.extflx' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
          ENDIF

          WRITE(numphpext,9500) kt,  glob_sum( bfpo4(:,:)*e1e2t(:,:) )
          WRITE(numfedext,9501) kt,  glob_sum( bffed(:,:)*e1e2t(:,:) ), glob_sum( dust_bling(:,:)*e1e2t(:,:)*tmask(:,:,1) )
          WRITE(numoxyext,9501) kt,  glob_sum( bfoxy(:,:)*e1e2t(:,:) ), glob_sum( o2flx     (:,:)*e1e2t(:,:)*tmask(:,:,1) )

        !ENDIF

9500  FORMAT(i10,e18.10)
9501  FORMAT(i10,2(e18.10))

      ENDIF

      IF (ln_diatrc) THEN

        IF (lk_iomput) THEN
          CALL iom_put(  "jp_sed", bfpo4  (:,:) )
          CALL iom_put( "jfe_sed", bffed  (:,:) )
          CALL iom_put( "jfe_bur", fpofe_b(:,:) )
          CALL iom_put( "jox_sed", bfoxy  (:,:) )
        ENDIF

        IF (  .NOT. lk_iomput ) THEN
          trc2d( :,:,jp_bling0_2d   ) = bfpo4  (:,:)*1.d3
          trc2d( :,:,jp_bling0_2d+1 ) = bffed  (:,:)*1.d3
          trc2d( :,:,jp_bling0_2d+2 ) = fpofe_b(:,:)*1.d3
          trc2d( :,:,jp_bling0_2d+3 ) = bfoxy  (:,:)*1.d3
        ENDIF

      ENDIF

      CALL wrk_dealloc( jpi, jpj, bfpo4, bffed, bfoxy, o2flx )

   END SUBROUTINE trc_ext_bling

   SUBROUTINE trc_ext_init_bling

      INTEGER  :: jm, ierr
      INTEGER  :: numdust, ntimes_dust
      REAL(wp), DIMENSION(12) :: zsteps                 ! times records

      REAL(wp) , ALLOCATABLE, DIMENSION(:,:,:) :: zdust
      !!----------------------------------------------------------------------

      ALLOCATE( sf_dust_bling(1), STAT=ierr )    
      IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'trc_ext_init_bling: unable to allocate sf_apr structure' )

      CALL fld_fill( sf_dust_bling, (/ sn_dust_bling /), cn_dir_bling, 'trc_ext_init_bling', 'Iron from sediment ', 'namblingiron' )

      ALLOCATE( sf_dust_bling(1)%fnow(jpi,jpj,1) )
      IF( sn_dust_bling%ln_tint ) ALLOCATE( sf_dust_bling(1)%fdta(jpi,jpj,1,2) )

      CALL iom_open (  TRIM( sn_dust_bling%clname ) , numdust )

      CALL iom_gettime( numdust, zsteps, kntime=ntimes_dust)  ! get number of record in file

      ALLOCATE( zdust(jpi,jpj,ntimes_dust) )
      DO jm = 1, ntimes_dust
         CALL iom_get( numdust, jpdom_data, TRIM( sn_dust_bling%clvar ), zdust(:,:,jm), jm )
      ENDDO
       
      CALL iom_close( numdust )
      DEALLOCATE( zdust)

   END SUBROUTINE trc_ext_init_bling
#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                        No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ext_bling 
      WRITE(*,*) 'trc_ext_bling: You should not have seen this print',  kt
   END SUBROUTINE trc_ext_bling

#endif

END MODULE trcext_blingv0
