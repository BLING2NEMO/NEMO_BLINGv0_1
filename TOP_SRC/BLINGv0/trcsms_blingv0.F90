MODULE trcsms_blingv0

#if defined key_bling
   !!======================================================================
   !!                     ***  MODULE trcsms_bling  ***
   !! TOP :  BLINGv0 model main routine. Computes the sources and sinks
   !!======================================================================
   !! History :  MClaret@McGill@04-07/2014. Ported from BLINGv0 in GFDL
   !!-----------------------------------------------------------------------
   USE par_trc         ! TOP parameters
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE trdmod_oce
   USE trdmod_trc
   USE iom

   ! BLINGv0 specific modules
   USE trcopt_blingv0
   USE trcext_blingv0
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_bling       ! called by trcsms.F90 module

   INTEGER :: numsms, numsms2, numphp, numfed, numoxy, numine ! logical unit for phosphate budget

#  include "top_substitute.h90"

CONTAINS

   SUBROUTINE trc_sms_bling( kt )
      !!------------------------------------------------------------------------
      !!                     ***  trc_sms_bling  ***
      !!
      !! ** History : default leap-frog scheme(MY_TRC) changed to FORWARD scheme
      !!------------------------------------------------------------------------
      
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      INTEGER  :: ji, jj, jk, jn

      REAL(wp) :: f_po4, f_dop, f_fer, f_oxy
      REAL(wp) :: ztra
      ! Irradiance k
      REAL(wp) :: po4_up, thetamax_fe, alpha_chl 
      ! Production
      REAL(wp) :: pc_tot, biomass_p_ts, theta, chl_dia, mulamb0expkT
      ! Phosphorous
      REAL(wp) :: frac_pop, jdop, fe2p_up, fpopkm1
      REAL(wp) :: zzz, wsink, oxy_up

      ! Iron
      REAL(wp) :: jpofe, fpofekm1
      REAL(wp) :: dum5, dum2, dum3

      REAL(wp), POINTER, DIMENSION(:,:,:) :: expkT, irr_inst, irr_mix, irrk, pc_m, mu
      REAL(wp), POINTER, DIMENSION(:,:,:) :: def_fe, feprime, kfe_eq_lig, fpofe
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jpop, fpop, zremin
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jp_uptake, jp_remin, jp_recycle
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jfe_uptake, jfe_remin, jfe_recycle
      REAL(wp), POINTER, DIMENSION(:,:,:) :: jfe_ads_inorg, jfe_ads_org
      REAL(wp), POINTER, DIMENSION(:,:,:) :: j_po4, j_dop, j_fed, j_oxy

      REAL(wp), POINTER, DIMENSION(:)     :: dum4
      REAL(wp), POINTER, DIMENSION(:,:,:) :: xnegtr
      REAL(wp), POINTER, DIMENSION(:,:,:) :: wrk1, wrk2, wrk3, wrk4

      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('trc_sms_bling')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_bling:  BLINGv0 model'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! allocate matrix variables
      CALL wrk_alloc( jpi, jpj, jpk, expkT, irr_inst, irr_mix, irrk, pc_m, mu )
      CALL wrk_alloc( jpi, jpj, jpk, def_fe, feprime, kfe_eq_lig, fpofe )
      CALL wrk_alloc( jpi, jpj, jpk, jpop, fpop, zremin )
      CALL wrk_alloc( jpi, jpj, jpk, jp_uptake,  jp_remin,  jp_recycle )
      CALL wrk_alloc( jpi, jpj, jpk, jfe_uptake, jfe_remin, jfe_recycle )
      CALL wrk_alloc( jpi, jpj, jpk, jfe_ads_inorg, jfe_ads_org )
      CALL wrk_alloc( jpi, jpj, jpk, j_po4, j_dop, j_fed, j_oxy )

      CALL wrk_alloc( jpk, dum4 )
      CALL wrk_alloc( jpi, jpj, jpk, xnegtr )
      CALL wrk_alloc( jpi, jpj, jpk, wrk1, wrk2, wrk3, wrk4 )

      !IF( (narea==1) .and. ln_bling_mass ) THEN      !   Write values for phosphate budget
      IF( ln_bling_mass ) THEN      !   Write values for phosphate budget
        CALL trc_sms_bling_mass_conserv (kt)
      ENDIF

      !!----------------------------------------------------------------------
      ! BLING model

      CALL trc_opt_bling_rb  (kt, irr_inst, irr_mix)  ! optical model (red and blue wavelengths)
      !wrk1=irr_inst
      !CALL trc_opt_bling_rgb (kt, irr_inst, irr_mix)  ! optical model (red, blue, and green wavelengths)
      !wrk1=irr_inst
      !wrk2=irr_mix

      DO ji=1, jpi
          DO jj=1, jpj

              dum4(:) = 0.d0

              DO jk=1, jpk

               ! ----------------------------------------------------------
               ! negative trophic variables DO not contribute to the fluxes
               ! ----------------------------------------------------------

               f_po4 = MAX( 0.e0, trn(ji,jj,jk,jpPO4_bling) )
               f_dop = MAX( 0.e0, trn(ji,jj,jk,jpDOP_bling) )
               f_fer = MAX( 0.e0, trn(ji,jj,jk,jpFed_bling) )
               f_oxy = trn(ji,jj,jk,jpOxy_bling)

               ! ----------------------------------------------------------
               ! TEMPERATURE DEPENDENCE
               ! NB: The temperature effect of Eppley (1972) is used instead 
               !     of that in Geider et al (1997) for both simplicity and 
               !     to incorporate combined effects on uptake, incorporation
               !     into organic matter and photorespiration.  Values of PCmax
               !     are normalized to 0C rather than 20C in Geider et al.(1997)
               ! ----------------------------------------------------------

               expkT(ji,jj,jk)=EXP(kappa_eppley*tsn(ji,jj,jk,jp_tem))

               ! ----------------------------------------------------------
               ! Phytoplankton are assumed to grow according to the general properties 
               ! described in Geider (1997). This formulation gives a biomass-specific 
               ! growthrate as a function of light, nutrient limitation, and 
               ! temperature. We modify this relationship slightly here, as described 
               ! below, and also use the assumption of steady state growth vs. loss to 
               ! derive a simple relationship between growth rate, biomass and uptake.
               ! ----------------------------------------------------------
               ! First, we calculate the limitation terms for PO4 and Fe, and the 
               ! Fe-limited Chl:C maximum.
               ! The light-saturated maximal photosynthesis rate term (pc_m) is simply 
               ! the product of a prescribed maximal photosynthesis rate (pc_0), the 
               ! Eppley temperature dependence, and a Liebig limitation (the minimum
               ! of Michaelis-Menton PO4-limitation, or iron-limitation).
               ! The iron limitation term has a lower limit of def_fe_min 
               ! and is scaled by (k_fe_2_p + fe_2_p_max) / fe_2_p_max
               ! so that it approaches 1 as fed approaches infinity. Thus, 
               ! it's of comparable magnitude to the PO4 limitation term.
               !
               ! Fe limitation acts in two additional mechanisms:
               ! 1. By reducing the maximum achievable Chl:C ratio 
               ! (theta) below a prescribed, Fe-replete maximum value (thetamax), to 
               ! approach a prescribed minimum Chl:C (thetamin) under extreme
               ! Fe-limitation.
               ! 2. By reducing alpha (the initial slope of the P-I curve) under Fe-
               ! limitation.
               ! ----------------------------------------------------------

               ! Iron uptake

               fe2p_up = fe2p_max * f_fer / (kfe + f_fer)
               def_fe(ji,jj,jk)  = def_fe_min + (1.d0-def_fe_min) &
                                  *fe2p_up/(kfe2p_up+fe2p_up)*(kfe2p_up+fe2p_max)/fe2p_max 

               ! Phosphate uptake
               po4_up = f_po4 /( kpo4 + f_po4 )

               ! Maximum production 
               pc_m(ji,jj,jk) = pc_0 * expkT(ji,jj,jk) * MIN(po4_up,def_fe(ji,jj,jk)) 


               ! Iron limitation on photosyntesis machinery
               thetamax_fe=thetamax_lo + (thetamax_hi - thetamax_lo)*def_fe(ji,jj,jk)
               alpha_chl  =alpha_min   + (alpha_max   - alpha_min  )*def_fe(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Next, the nutrient-limited efficiency of algal photosystems, Irrk, is
               ! calculated. This requires a prescribed quantum yield, alpha.
               ! The iron deficiency term is included here as a multiplier of the 
               ! thetamax_fe to represent the importance of Fe in forming chlorophyll
               ! accessory antennae, which do not affect the Chl:C but still affect the
               ! phytoplankton ability to use light (eg Stzrepek & Harrison Nature 
               ! 2004).
               !-----------------------------------------------------------------------

               ! I_k
               irrk(ji,jj,jk) = pc_m(ji,jj,jk) / (epsln + alpha_chl*thetamax_fe) + irr_mem(ji,jj,jk)*0.5d0

               !-----------------------------------------------------------------------
               ! Now we can calculate the carbon-specific photosynthesis rate, pc_tot.
               !-----------------------------------------------------------------------

               pc_tot = pc_m(ji,jj,jk)*(1.d0-EXP(-irr_mix(ji,jj,jk)/(irrk(ji,jj,jk)+epsln)))

               !-----------------------------------------------------------------------
               ! Next, we account for the maintenance effort that phytoplankton must 
               ! exert in order to combat decay. This is prescribed as a fraction of the
               ! light-saturated photosynthesis rate, resp_frac. The result of this is 
               ! to set a level of energy availability below which net growth (and 
               ! therefore nutrient uptake) is zero, given by resp_frac * pc_m.
               !-----------------------------------------------------------------------

               ! Net total production
               mu(ji,jj,jk) = MAX (0.d0,pc_tot-resp_frac*pc_m(ji,jj,jk))

               !-----------------------------------------------------------------------
               ! We now must convert this net carbon-specific growth rate to nutrient 
               ! uptake rates, the quantities we are interested in. Since we have no 
               ! explicit biomass tracer, we use the result of Dunne et al. (GBC, 2005) 
               ! to calculate an implicit biomass from the uptake rate through the  
               ! application of a simple idealized grazing law. This has the effect of 
               ! reducing uptake in low growth-rate regimes and increasing uptake in 
               ! high growth-rate regimes - essentially a non-linear amplification of 
               ! the growth rate variability. The result is:
               !-----------------------------------------------------------------------

               ! Biomass
               mulamb0expkT = mu(ji,jj,jk)/(lambda0*expkT(ji,jj,jk))
               biomass_p_ts = p_star*mulamb0expkT*(1.d0+(mulamb0expkT)**2)

               IF (kt==nittrc000) biomass_p(ji,jj,jk)=biomass_p_ts

               biomass_p(ji,jj,jk) =   biomass_p(ji,jj,jk) &
                                    + (biomass_p_ts-biomass_p(ji,jj,jk))*MIN(1.d0,gam_biomass*rfact)*tmask(ji,jj,jk)

               !if (ji==80 .and. jj==60 .and. jk==1) write(*,'(I3,5(1X,E11.4))') kt, &
               !pc_0, expkT(ji,jj,jk), po4_up, def_fe(ji,jj,jk), pc_tot,mu(ji,jj,jk),biomass_p(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! We can now use the diagnostic biomass to calculate the chlorophyll
               ! concentration:
               !-----------------------------------------------------------------------

               ! Chl:C ration
               theta   = thetamax_fe / (1.d0 + (thetamax_fe*alpha_chl*irr_mem(ji,jj,jk))/(2.d0*pc_m(ji,jj,jk)+epsln) )

               ! Chl biomass
               chl_dia = biomass_p(ji,jj,jk) * c2p * 12.011e+6 * theta * tmask(ji,jj,jk) 
               chl_bling(ji,jj,jk) = MAX(chl_min, chl_dia)

               !--------------------------------------------------
               ! PHOSPHORUS CYCLE
               !--------------------------------------------------
               ! The uptake of nutrients is assumed to contribute to the growth of
               ! phytoplankton, which subsequently die and are consumed by heterotrophs.
               ! This can involve the transfer of nutrient elements between many
               ! organic pools, both particulate and dissolved, with complex histories.
               ! We take a simple approach here, partitioning the total uptake into two
               ! fractions - sinking and non-sinking - as a function of temperature, 
               ! following Dunne et al. (2005). 
               ! Then, the non-sinking fraction is further subdivided, such that the 
               ! majority is recycled instantaneously to the inorganic nutrient pool,
               ! representing the fast turnover of labile dissolved organic matter via
               ! the microbial loop, and the remainder is converted to semi-labile
               ! dissolved organic matter. Iron and phosphorus are treated identically 
               ! for the first step, but all iron is recycled instantaneously in the
               ! second step (i.e. there is no dissolved organic iron pool).
               !-----------------------------------------------------------------------

               ! Phosphorous uptake flux
               jp_uptake(ji,jj,jk)=biomass_p(ji,jj,jk)*mu(ji,jj,jk)

               ! 
               frac_pop=(phi_sm+phi_lg*(mulamb0expkT)**2) / (1+(mulamb0expkT)**2) * EXP(kappa_remin*tsn(ji,jj,jk,jp_tem))

               !
               jpop(ji,jj,jk)=frac_pop*jp_uptake(ji,jj,jk)

               !-----------------------------------------------------------------------
               ! Then the remainder is divided between instantaneously recycled and
               ! long-lived dissolved organic matter,
               !-----------------------------------------------------------------------
               !
               jdop=phi_dop*(jp_uptake(ji,jj,jk)-jpop(ji,jj,jk))

               ! 
               jp_recycle(ji,jj,jk)=jp_uptake(ji,jj,jk)-jpop(ji,jj,jk)-jdop

               !---------------------------------------------------------------------
               ! IRON
               !---------------------------------------------------------------------
               ! Iron is then taken up as a function of PO4 uptake and iron limitation,
               ! with a maximum Fe:P uptake ratio of fe2p_max:
               !-----------------------------------------------------------------------

               ! Iron uptake flux
               jfe_uptake(ji,jj,jk) =jp_uptake(ji,jj,jk)*fe2p_up
               jpofe                =frac_pop*jfe_uptake(ji,jj,jk)
               jfe_recycle(ji,jj,jk)=jfe_uptake(ji,jj,jk)-jpofe

               !-----------------------------------------------------------------------
               ! Calculate free and inorganically associated iron concentration for
               ! scavenging.
               ! We assume that there is a 
               ! spectrum of iron ligands present in seawater, with varying binding
               ! strengths and whose composition varies with light and iron 
               ! concentrations. For example, photodissocation of ligand complexes 
               ! occurs under bright light, weakening the binding strength 
               ! (e.g. Barbeau et al., Nature 2001), while at very low iron 
               ! concentrations (order kfe_eq_lig_femin), siderophores are thought
               ! to be produced as a response to extreme
               ! iron stress.
               ! In anoxic waters, iron should be reduced, and therefore mostly 
               ! immune to scavenging. Easiest way to do this is to skip the feprime
               ! calculation if oxygen is less than 0.
               !-----------------------------------------------------------------------

               if (f_oxy > oxy_min ) then
                 dum5       = irr_inst(ji,jj,jk)**2/(kfe_eq_lig_irr**2+irr_inst(ji,jj,jk)**2)
                 dum2       = max(  0.e0, min( 1.e0, 1.2e0*(f_fer-kfe_eq_lig_femin)/(epsln+f_fer) )  )
                 kfe_eq_lig(ji,jj,jk) = kfe_eq_lig_max -(kfe_eq_lig_max-kfe_eq_lig_min)*dum5*dum2
                 feprime(ji,jj,jk)= 1.e0 + kfe_eq_lig(ji,jj,jk)*(felig_bkg - f_fer)
                 feprime(ji,jj,jk)= (-feprime(ji,jj,jk) + sqrt(feprime(ji,jj,jk)**2 + 4.e0*kfe_eq_lig(ji,jj,jk)*f_fer) )/(2.e0*kfe_eq_lig(ji,jj,jk))
               else
                 feprime(ji,jj,jk) = 0.e0
               endif

               jfe_ads_inorg(ji,jj,jk) = min( 0.5d0/rfact, kfe_inorg*sqrt(feprime(ji,jj,jk)) )*feprime(ji,jj,jk)

               dum4(jk) = jpofe + jfe_ads_inorg(ji,jj,jk)

               !--------------------------------------------------
               ! COMPUTE TRENDS w/o remineralization processes
               !--------------------------------------------------
               !
               j_po4(ji,jj,jk) =   jp_recycle(ji,jj,jk) + gamma_dop*f_dop -jp_uptake(ji,jj,jk)
               j_dop(ji,jj,jk) = - gamma_dop*f_dop + phi_dop*(jp_uptake(ji,jj,jk)-jpop(ji,jj,jk))
               j_fed(ji,jj,jk) =   jfe_recycle(ji,jj,jk)-jfe_uptake(ji,jj,jk)-jfe_ads_inorg(ji,jj,jk)

wrk1(ji,jj,jk)=gamma_dop*f_dop

               ! checks
               ! wrk1(ji,jj,jk)=mulamb0expkT 
               ! wrk2(ji,jj,jk)=irrk
               ! wrk3(ji,jj,jk)=f_po4+f_dop
               ! wrk4(ji,jj,jk)=0.d0

            ENDDO


            !-----------------------------------------------------------------------
            ! SINKING AND REMINERALIZATION
            !-----------------------------------------------------------------------
            ! Calculate the remineralization lengthscale matrix, zremin, a function 
            ! of z. Sinking rate (wsink) is constant over the upper wsink0_z metres,
            ! then  increases linearly with depth.
            ! The remineralization rate is a function of oxygen concentrations,
            ! following a Holling type 2 dependence, decreasing to a minimum value
            ! of remin_min. This is ad hoc, following work by Bianchi, Sarmiento,
            ! Galbraith and Kwon (unpublished).
            !-----------------------------------------------------------------------
            ! In general, the flux at the bottom of a grid cell should equal
            ! Fb = (Ft + Prod*dz) / (1 + zremin*dz)
            ! where Ft is the flux at the top, and prod*dz is the integrated 
            ! production of new sinking particles within the layer.
            ! Since Ft=0 in the first layer,
            !---------------------------------------------------------------------

            ! k=1: surface layer
            zzz=fse3t(ji,jj,1)

            IF (zzz .lt. wsink0_z) THEN
               wsink=wsink0
            ELSE
               wsink=wsink0+wsink_acc*(zzz-wsink0_z)
            ENDIF

            ! Phosphorous
            f_oxy = trn(ji,jj,1,jpOxy_bling)
            oxy_up =f_oxy**2 / (koxy**2 + f_oxy**2)
            zremin(ji,jj,1) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)

            fpop(ji,jj,1)    = jpop(ji,jj,1)*fse3t(ji,jj,1)/(1.d0+fse3t(ji,jj,1)*zremin(ji,jj,1)) 
            jp_remin(ji,jj,1)=(jpop(ji,jj,1)*fse3t(ji,jj,1)-fpop(ji,jj,1))/(epsln+fse3t(ji,jj,1))

            !-----------------------------------------------------------------------
            ! Now, calculate the Fe adsorption using this fpop:
            ! The absolute first order rate constant is calculated from the 
            ! concentration of organic particles, after Parekh et al. (2005). Never
            !  allowed to be greater than 1/2dt for numerical stability.
            !-----------------------------------------------------------------------

            !Iron
            dum3                 = (fpop(ji,jj,1)*c2p*12.011d0/(epsln+wsink))**(0.58)
            jfe_ads_org(ji,jj,1) = min (0.5d0/rfact, kfe_org*dum3)*feprime(ji,jj,1)
            dum4(1)              = (dum4(1)+jfe_ads_org(ji,jj,1))*fse3t(ji,jj,1)
            fpofe(ji,jj,1)       = dum4(1)/(1.d0+fse3t(ji,jj,1)*zremin(ji,jj,1))
            jfe_remin(ji,jj,1)   = (dum4(1)-fpofe(ji,jj,1))/(epsln+fse3t(ji,jj,1))

            ! Add remineralization terms to trends
            j_po4(ji,jj,1)=j_po4(ji,jj,1)+(1.d0-phi_dop)*jp_remin(ji,jj,1)
            j_dop(ji,jj,1)=j_dop(ji,jj,1)+       phi_dop*jp_remin(ji,jj,1)
            j_fed(ji,jj,1)=j_fed(ji,jj,1)+jfe_remin(ji,jj,1)-jfe_ads_org(ji,jj,1)   

            ! k=2:NK: rest of the water column
            DO jk=2, jpk

               fpopkm1 = fpop(ji,jj,jk-1)
               fpofekm1=fpofe(ji,jj,jk-1)
    
               zzz=zzz+fse3t(ji,jj,jk)

               IF (zzz .lt. wsink0_z) THEN
                  wsink=wsink0
               ELSE
                  wsink=wsink0+wsink_acc*(zzz-wsink0_z)
               ENDIF

               ! Phosphorous
               f_oxy = trn(ji,jj,jk,jpOxy_bling)
               oxy_up =f_oxy**2 / (koxy**2 + f_oxy**2)
               zremin(ji,jj,jk) =gamma_pop*(oxy_up*(1.d0-remin_min)+remin_min)/(epsln+wsink)

               fpop(ji,jj,jk)    =(fpopkm1+jpop(ji,jj,jk)*fse3t(ji,jj,jk))/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk)) 
               jp_remin(ji,jj,jk)=(fpopkm1+jpop(ji,jj,jk)*fse3t(ji,jj,jk)-fpop(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))
          
               ! Iron
               dum3            = (fpop(ji,jj,jk)*c2p*12.011d0/(epsln+wsink))**(0.58)
               jfe_ads_org(ji,jj,jk) = min (0.5d0/rfact, kfe_org*dum3)*feprime(ji,jj,jk)
               dum4(jk)        = (dum4(jk)+jfe_ads_org(ji,jj,jk))*fse3t(ji,jj,jk)
               fpofe(ji,jj,jk) = (fpofekm1 + dum4(jk))/(1.d0+fse3t(ji,jj,jk)*zremin(ji,jj,jk))
               jfe_remin(ji,jj,jk) = (fpofekm1+dum4(jk)-fpofe(ji,jj,jk))/(epsln+fse3t(ji,jj,jk))

               ! Save fPOP and fPOFe at the bottom grid cell to compute bottom fluxes
               IF (jk==mbkt(ji,jj)) THEN
                  fpop_b (ji,jj) = fpop(ji,jj,jk)
                  fpofe_b(ji,jj) = fpofe(ji,jj,jk)
               ENDIF
           
               ! Add remineralization terms to trends
               j_po4(ji,jj,jk)=j_po4(ji,jj,jk)+(1.d0-phi_dop)*jp_remin(ji,jj,jk)
               j_dop(ji,jj,jk)=j_dop(ji,jj,jk)+      phi_dop *jp_remin(ji,jj,jk)
               j_fed(ji,jj,jk)=j_fed(ji,jj,jk)+jfe_remin(ji,jj,jk)-jfe_ads_org(ji,jj,jk)

            ENDDO 

            !-----------------------------------------------------------------------
            ! OXYGEN
            !-----------------------------------------------------------------------
            ! Assuming constant P:O ratio.
            ! Optional prevention of negative oxygen (does not conserve ocean 
            ! redox potential) or alternatively it can be allowed to go negative, 
            ! keeping track of an implicit nitrate deficit 
            ! plus sulfate reduction.
            !-----------------------------------------------------------------------

            DO jk=1, jpk
               f_oxy = trn(ji,jj,jk,jpOxy_bling)
               IF ( (ln_prev_o2lt0) .and. (f_oxy<oxy_min) ) then
                   j_oxy(ji,jj,jk)=0.d0
               ELSE
                   j_oxy(ji,jj,jk)=-oxy2p*j_po4(ji,jj,jk)
               ENDIF
            ENDDO

         ENDDO
      ENDDO

!ji=80;jj=60;jk=1
!write(*,'(I3,5E11.4)') kt, j_dop(ji,jj,jk), -wrk1(ji,jj,jk), phi_dop*(jp_uptake(ji,jj,jk)) &
!                        , -phi_dop*jpop(ji,jj,jk), phi_dop *jp_remin(ji,jj,jk)


      !write(*,'(I3,3(1X,E11.4))') kt,tra(ji,jj,jk,jpPO4_bling),j_po4(ji,jj,jk),rfact
      !write(*,'(I3,3(1X,E11.4))') kt,tra(ji,jj,jk,jpDOP_bling),j_dop(ji,jj,jk),rfact

      tra(:,:,:,jpPO4_bling) = tra(:,:,:,jpPO4_bling) + j_po4(:,:,:)*rfact
      tra(:,:,:,jpDOP_bling) = tra(:,:,:,jpDOP_bling) + j_dop(:,:,:)*rfact
      tra(:,:,:,jpFed_bling) = tra(:,:,:,jpFed_bling) + j_fed(:,:,:)*rfact
      tra(:,:,:,jpOxy_bling) = tra(:,:,:,jpOxy_bling) + j_oxy(:,:,:)*rfact
      !tra(:,:,:,jpine_bling) = 10.e0/(86400.e0*365.e0)*rfact
      tra(:,:,:,jpine_bling) = 0.e0

      !tra(:,:,:,jpPO4_bling) = tra(:,:,:,jpPO4_bling) + j_po4(:,:,:)
      !tra(:,:,:,jpDOP_bling) = tra(:,:,:,jpDOP_bling) + j_dop(:,:,:)
      !tra(:,:,:,jpFed_bling) = tra(:,:,:,jpFed_bling) + j_fed(:,:,:)
      !tra(:,:,:,jpOxy_bling) = tra(:,:,:,jpOxy_bling) + j_oxy(:,:,:)
      !tra(:,:,:,jpine_bling) = 10.e0/5760.d0

      dum1(:,:,20)=dum1(:,:,20)+glob_sum ( j_oxy(:,:,:)*cvol(:,:,:) )*rfact

      !test if concentrations fall below 0
      xnegtr(:,:,:) = 1.e0
      DO jn = jp_bling0, jp_bling1
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  IF( ( trn(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) < 0.e0 ) THEN 
                     ztra             = ABS(  ( trn(ji,jj,jk,jn) - rtrn ) &
                                            / ( tra(ji,jj,jk,jn) + rtrn ) )
                     xnegtr(ji,jj,jk) = MIN( xnegtr(ji,jj,jk),  ztra )
                  ENDIF
              END DO
            END DO
         END DO
      END DO

      ! Prgonostic tracer fields
      DO jn = jp_bling0, jp_bling1
         !write(*,'(I3,2(1X,E11.4))') jp_bling0, trn(ji,jj,jk,jn),tra(ji,jj,jk,jn)
         trn(:,:,:,jn) = trn(:,:,:,jn) + xnegtr(:,:,:) * tra(:,:,:,jn)
      END DO

      !DO jn = jp_bling0, jp_bling1
      !   trn(:,:,:,jn) = trn(:,:,:,jn) + tra(:,:,:,jn)
      !   DO jk = 1, jpk
      !      DO jj = 1, jpj
      !         DO ji = 1, jpi
      !           trn(ji,jj,jk,jn)=MAX( 0.e0, trn(ji,jj,jk,jn) )
      !         END DO
      !      END DO
      !   END DO
      !END DO

      ! Copy new arrays to trb (tracer fields before) and set tra to zero
      ! to compute tracer gradients with tracer fields after ecological forcing
      tra(:,:,:,jp_bling0:jp_bling1) = 0.e0

      ! add external fluxes
      CALL trc_ext_bling (kt)

      ! Copy to trb to use BLING updated tracer values to compute transport
      ! trends
      DO jn=jp_bling0, jp_bling1
         trb(:,:,:,jn)=trn(:,:,:,jn)
      ENDDO

      DO jn=jp_bling0, jp_bling1
         CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )
         CALL lbc_lnk( trb(:,:,:,jn), 'T', 1. )
         CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )
      ENDDO

      !checks
      !jk=20
      !write(*,'(I3,3F14.7)'), kt, trb(60,60,jk,jpPO4_bling)*1.d6, &
      !                            trn(60,60,jk,jpPO4_bling)*1.d6, &
      !                            tra(60,60,jk,jpPO4_bling)*1.d6
      !IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
      !   WRITE(charout, FMT="('bio')")
      !   CALL prt_ctl_trc_info(charout)
      !   CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      !ENDIF

wrk2(:,:,:)=cvol(:,:,:)
wrk3(:,:,1)=hmld (:,:)
wrk3(:,:,2)=hmlp(:,:)

      !IF ( ln_bling_mass .and.  narea == 1 ) THEN
      IF ( ln_bling_mass  ) THEN

        IF( kt == nittrc000 ) THEN 
          CALL ctl_opn( numsms ,  'sms.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
          CALL ctl_opn( numsms2, 'sms2.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        ENDIF

        WRITE(UNIT=numsms,FMT='(i10,6(3x,e18.10))')  kt  &
                           , glob_sum (jp_recycle(:,:,:)*cvol(:,:,:))*rfact &
              , glob_sum (jp_remin(:,:,:)*cvol(:,:,:))*(1.d0-phi_dop)*rfact &
                           , glob_sum (      wrk1(:,:,:)*cvol(:,:,:))*rfact &
                           , glob_sum (-jp_uptake(:,:,:)*cvol(:,:,:))*rfact &
      , glob_sum (phi_dop*(jp_uptake(:,:,:)-jpop(:,:,:))*cvol(:,:,:))*rfact &
              , glob_sum (jp_remin(:,:,:)*     (phi_dop)*cvol(:,:,:))*rfact

        WRITE(UNIT=numsms2,FMT='(i10,2(3x,e18.10))')  kt  &
                          , glob_sum (j_po4(:,:,:)*cvol(:,:,:))*rfact &
                          , glob_sum (j_dop(:,:,:)*cvol(:,:,:))*rfact 

      ENDIF

      IF( lk_iomput ) THEN
            CALL iom_put( "wrk1", wrk1(:,:,:)  )  
            CALL iom_put( "wrk2", wrk2(:,:,:)  )  
            CALL iom_put( "wrk3", wrk3(:,:,:)  )  
            CALL iom_put( "wrk4", wrk4(:,:,:)  )  
      ENDIF

      IF (ln_diatrc) THEN
         IF (lk_iomput) THEN
            CALL iom_put(      "expkT"  ,        expkT(:,:,:)*tmask(:,:,:) )
            CALL iom_put(   "irr_inst"  ,     irr_inst(:,:,:)*tmask(:,:,:) )
            CALL iom_put(    "irr_mix"  ,      irr_mix(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "irrk"  ,         irrk(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "pc_m"  ,         pc_m(:,:,:)*tmask(:,:,:) )
            CALL iom_put(         "mu"  ,           mu(:,:,:)*tmask(:,:,:) )
            CALL iom_put(  "biomass_p"  ,    biomass_p(:,:,:)*tmask(:,:,:) )
            CALL iom_put(       "fpop"  ,         fpop(:,:,:)*tmask(:,:,:) )
            CALL iom_put(     "zremin"  ,       zremin(:,:,:)*tmask(:,:,:) )
            CALL iom_put(     "def_fe"  ,       def_fe(:,:,:)*tmask(:,:,:) )
            CALL iom_put(    "feprime"  ,      feprime(:,:,:)*tmask(:,:,:) )
            CALL iom_put( "kfe_eq_lig"  ,   kfe_eq_lig(:,:,:)*tmask(:,:,:) )
            CALL iom_put(      "fpofe"  ,        fpofe(:,:,:)*tmask(:,:,:) )

            CALL iom_put(       "jpop"  ,         jpop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(  "jp_uptake"  ,    jp_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(   "jp_remin"  ,     jp_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put( "jp_recycle"  ,   jp_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                           
            CALL iom_put( "jfe_uptake"  ,   jfe_uptake(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put(  "jfe_remin"  ,    jfe_remin(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put( "jfe_recycle" ,  jfe_recycle(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put( "jfe_ads_org" ,  jfe_ads_org(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                           
            CALL iom_put("jfe_ads_inorg",jfe_ads_inorg(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )                          
            CALL iom_put(        "jpo4" ,        j_po4(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "jdop" ,        j_dop(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "jfed" ,        j_fed(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
            CALL iom_put(        "joxy" ,        j_oxy(:,:,:)*fse3t(:,:,:)*tmask(:,:,:) )
         ENDIF

         IF ( .NOT. lk_iomput )  THEN
           trc3d( :,:,:,jp_bling0_3d    ) =    chl_bling(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+1  ) =        j_po4(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+2  ) =        j_dop(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+3  ) =    jp_uptake(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+4  ) =   jp_recycle(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+5  ) =     jp_remin(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+6  ) =        j_fed(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+7  ) =   jfe_uptake(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+8  ) =  jfe_recycle(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+9  ) =    jfe_remin(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+10 ) =  jfe_ads_org(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+11 ) =jfe_ads_inorg(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+12 ) =        j_oxy(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)

           trc3d( :,:,:,jp_bling0_3d+13 ) =      jpop(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+14 ) =      fpop(:,:,:)*fse3t(:,:,:)*1.d3*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+15 ) =     fpofe(:,:,:)*fse3t(:,:,:)*1.d6*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+16 ) =     expkT(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+17 ) =  irr_inst(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+18 ) =   irr_mix(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+19 ) =      irrk(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+20 ) =      pc_m(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+21 ) =        mu(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+22 ) = biomass_p(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+23 ) =    zremin(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+24 ) =    def_fe(:,:,:)*tmask(:,:,:)
           trc3d( :,:,:,jp_bling0_3d+25 ) =   feprime(:,:,:)*tmask(:,:,:)*1.d6
           trc3d( :,:,:,jp_bling0_3d+26 ) =kfe_eq_lig(:,:,:)*tmask(:,:,:)
         ENDIF
      ENDIF

      CALL wrk_dealloc( jpi, jpj, jpk, expkT, irr_inst, irr_mix, irrk, pc_m, mu )
      CALL wrk_dealloc( jpi, jpj, jpk, def_fe, feprime, kfe_eq_lig, fpofe )
      CALL wrk_dealloc( jpi, jpj, jpk, jpop, fpop, zremin )
      CALL wrk_dealloc( jpi, jpj, jpk, jp_uptake, jp_remin, jp_recycle )
      CALL wrk_dealloc( jpi, jpj, jpk, jfe_uptake, jfe_remin, jfe_recycle )
      CALL wrk_dealloc( jpi, jpj, jpk, jfe_ads_inorg, jfe_ads_org )
      CALL wrk_dealloc( jpi, jpj, jpk, j_po4, j_dop, j_fed, j_oxy )

      CALL wrk_dealloc( jpk, dum4 )
      CALL wrk_dealloc( jpi, jpj, jpk, xnegtr  )
      CALL wrk_dealloc( jpi, jpj, jpk, wrk1, wrk2, wrk3, wrk4 )

      IF( nn_timing == 1 )  CALL timing_stop('trc_sms_bling')
      !
   END SUBROUTINE trc_sms_bling

   SUBROUTINE trc_sms_bling_mass_conserv (kt)

      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      REAL(wp) :: sum_phosp, sum_fed, sum_oxy, sum_inert
      
      !!-------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN 

        CALL ctl_opn( numphp, 'php.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numfed, 'fed.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numoxy, 'oxy.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )
        CALL ctl_opn( numine, 'ine.budget' , 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, 6, .FALSE., narea )

        WRITE(numphp,9500) kt,  areatot
        WRITE(numfed,9500) kt,  areatot
        WRITE(numoxy,9500) kt,  areatot
        WRITE(numine,9500) kt,  areatot

      ENDIF

      ! total mass of phosphate
      sum_phosp = glob_sum ( ( trn(:,:,:,jpPO4_bling) + trn(:,:,:,jpDOP_bling) )*cvol(:,:,:)  )
      sum_fed   = glob_sum (   trn(:,:,:,jpFed_bling)*cvol(:,:,:)  )
      sum_oxy   = glob_sum (   trn(:,:,:,jpOxy_bling)*cvol(:,:,:)  )
      sum_inert = glob_sum (   trn(:,:,:,jpine_bling)*cvol(:,:,:)  )

      dum1(:,:, 1)=sum_phosp
      dum1(:,:, 2)=sum_phosp/areatot
      dum1(:,:, 4)=sum_fed
      dum1(:,:, 5)=sum_fed/areatot
      dum1(:,:, 9)=sum_oxy
      dum1(:,:,10)=sum_oxy/areatot
      dum1(:,:,13)=sum_inert
      dum1(:,:,14)=sum_inert/areatot

      IF( lwp ) THEN
        WRITE(numphp,9500) kt,  sum_phosp
        WRITE(numfed,9500) kt,  sum_fed
        WRITE(numoxy,9500) kt,  sum_oxy
        WRITE(numine,9500) kt,  sum_inert
      ENDIF

9500  FORMAT(i10,e18.10)    

   END SUBROUTINE trc_sms_bling_mass_conserv

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                                       No BLINGv0 model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sms_bling             ! Empty routine
      WRITE(*,*) 'trc_sms_bling: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms_bling
#endif
   !!======================================================================
END MODULE trcsms_blingv0
