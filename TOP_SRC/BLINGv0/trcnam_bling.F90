MODULE trcnam_bling
   !!======================================================================
   !!                      ***  MODULE trcnam_bling  ***
   !! TOP :   initialisation of some run parameters for LOBSTER bio-model
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec) Original code
   !!----------------------------------------------------------------------
#if defined key_bling
   !!----------------------------------------------------------------------
   !!   'key_bling'   :                                       BLINGv0 model
   !!----------------------------------------------------------------------
   !! trc_nam_bling      : BLINGv0 model initialisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE par_trc         ! TOP parameters
   USE trc             ! TOP variables
   USE iom
   USE vars_bling

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_nam_bling   ! called by trcnam.F90 module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id: trcnam_bling.F90 -1   $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nam_bling
      !!----------------------------------------------------------------------
      !!                     ***  trc_nam_bling  ***  
      !!
      !! ** Purpose :   read BLINGv0 namelist
      !!
      !!----------------------------------------------------------------------
      !
      INTEGER :: numnatg, jl, jn

      TYPE(DIAG), DIMENSION(jp_bling_2d )  :: blidia2d
      TYPE(DIAG), DIMENSION(jp_bling_3d )  :: blidia3d
      TYPE(DIAG), DIMENSION(jp_bling_trd ) :: blidiabio

      NAMELIST/namblingrat/    c2n, c2p, oxy2p
      NAMELIST/namblingopt/    xkr0, xkb0, xkrp, xkbp, xlr, xlb, rpig, rcchl, gam_irr_mem
      NAMELIST/namblingprod/   pc_0, kappa_eppley, kpo4, kfe, fe2p_max, kfe2p_up, def_fe_min, thetamax_lo, thetamax_hi
      NAMELIST/namblingprod/   alpha_min, alpha_max, resp_frac, p_star, lambda0, gam_biomass
      NAMELIST/namblingremin/  wsink0_z, wsink0, wsink_acc, koxy, remin_min, phi_dop, phi_sm
      NAMELIST/namblingremin/  phi_lg, kappa_remin, gamma_dop, gamma_pop
      NAMELIST/namblingairsea/ a_0, a_1, a_2, a_3, a_4, a_5, b_0, b_1, b_2, b_3, c_0
      NAMELIST/namblingairsea/ a_1_o2, a_2_o2, a_3_o2, a_4_o2
      NAMELIST/namblingiron/   oxy_min, kfe_eq_lig_max, kfe_eq_lig_min, kfe_eq_lig_irr, kfe_eq_lig_femin
      NAMELIST/namblingiron/   felig_bkg, kfe_inorg, kfe_org, ln_prev_o2lt0
      NAMELIST/namblingiron/   cn_dir_bling, sn_dust_bling, ln_dust_bling
      NAMELIST/namblingbudget/ ln_bling_mass, ln_bling_ext               ! additional diagnostics
      NAMELIST/namblingdia/    blidia3d, blidia2d                        ! additional diagnostics
      NAMELIST/namblingdbi/    blidiabio                                 ! additional diagnostics
      !
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_nam_bling : read BLINGv0 namelists'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

      CALL ctl_opn( numnatg, 'namelist_blingv0', 'OLD', 'FORMATTED','SEQUENTIAL', -1, numout, .FALSE. )

      ! namblingrat  : Stochiometric ratios
      REWIND (numnatg)
      READ   (numnatg, namblingrat)

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingrat:'
         WRITE(numout,*) '                           c2n     = ', c2n
         WRITE(numout,*) '                           c2p     = ', c2p
         WRITE(numout,*) '                           oxy2p   = ', oxy2p
      ENDIF

      ! namblingprod : Production parameters
      REWIND (numnatg)
      READ   (numnatg, namblingprod)

      lambda0    =    lambda0/86400.d0
      gam_biomass=gam_biomass/86400.d0

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingprod:'
         WRITE(numout,*) '    max growth rate at 0 C                          pc_0        = ', pc_0
         WRITE(numout,*) '    T dependence of growth                          kappa_epply = ', kappa_eppley 
         WRITE(numout,*) '    PO4 uptake half-saturation cte                  kpo4        = ', kpo4
         WRITE(numout,*) '    dissolved Fe uptake half-saturation cte         kfe         = ', kfe
         WRITE(numout,*) '    minimum value for iron deficiency term          def_fe_min  = ', def_fe_min
         WRITE(numout,*) '    maximum Fe:P uptake ratio                       fe2p_up_max = ', fe2p_max 
         WRITE(numout,*) '    half-saturation cellular Fe:P                   kfe2p_up    = ', kfe2p_up
         WRITE(numout,*) '    max Chl:C ratio, extreme iron limitation        thetamax_lo = ', thetamax_lo
         WRITE(numout,*) '    max Chl:C ratio, abundant iron                  thetamax_hi = ', thetamax_hi
         WRITE(numout,*) '    quantum yield under low light, iron limited     alpha_min   = ', alpha_min 
         WRITE(numout,*) '    quantum yield under low light, abundant iron    alpha_max   = ', alpha_max
         WRITE(numout,*) '    fraction of gross production respirated         resp_frac   = ', resp_frac
         WRITE(numout,*) '    pivotal phytoplankton biomass                   p_star      = ', p_star
         WRITE(numout,*) '    carbon-specific phytoplankton mortality rate    lambda0     = ', lambda0
         WRITE(numout,*) '    biomass adjustment time constant                gam_biomass = ', gam_biomass
      ENDIF

      ! namblingopt : Optical parameters
      REWIND (numnatg)
      READ   (numnatg, namblingopt)

      gam_irr_mem=gam_irr_mem/86400.d0

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingopt:'
         WRITE(numout,*) '    green   water absorption coeff                       xkg0  = ', xkb0
         WRITE(numout,*) '    red water absorption coeff                           xkr0  = ', xkr0
         WRITE(numout,*) '    pigment red absorption coeff                         xkrp  = ', xkrp
         WRITE(numout,*) '    pigment green absorption coeff                       xkgp  = ', xkbp
         WRITE(numout,*) '    green chl exposant                                   xlg   = ', xlb
         WRITE(numout,*) '    red   chl exposant                                   xlr   = ', xlr
         WRITE(numout,*) '    chla/chla+phea ratio                                 rpig  = ', rpig
         WRITE(numout,*) '    Photoadaptation time constant                 gam_irr_mem  = ', gam_irr_mem
      ENDIF

      ! namblingremin : Remineralization parameters
      REWIND (numnatg)
      READ   (numnatg, namblingremin)

      wsink0   =wsink0/86400.d0
      wsink_acc=wsink_acc/86400.d0
      gamma_dop=gamma_dop/365.25/86400.d0
      gamma_pop=gamma_pop/86400.d0

      IF (lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingremin:'
         WRITE(numout,*) '    depth at which sinking rate starts increasing    wsink0_z    =', wsink0_z 
         WRITE(numout,*) '    initial sinking rate                             wsink0      =', wsink0
         WRITE(numout,*) '    accerelation rate of sinking with depth          wsink_acc   =', wsink_acc
         WRITE(numout,*) '    half saturation const for aerobic respiration    koxy        =', koxy
         WRITE(numout,*) '    minimum anaerobic respiration rate               remin_min   =', remin_min 
         WRITE(numout,*) '    fraction of non-particulate uptake to DOM        phi_dop     =', phi_dop
         WRITE(numout,*) '    detritus production by small phyto               phi_sm      =', phi_sm
         WRITE(numout,*) '    detritus production by large phyto               phi_lg      =', phi_lg
         WRITE(numout,*) '    T dependence of particulate production           kappa_remin =', kappa_remin
         WRITE(numout,*) '    decay timescale of DOM                           gamma_dop   =', gamma_dop
         WRITE(numout,*) '    remineralization rate of sinking POM             gamma_pop   =', gamma_pop
         !WRITE(numout,*) ' '
      ENDIF

      ! namblingairsea : Air-sea interaction parameters
      REWIND (numnatg)
      READ   (numnatg, namblingairsea)

      IF (lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingairsea:'
         WRITE(numout,*) '    saturation coefficient for O2    a_0 =', a_0
         WRITE(numout,*) '    saturation coefficient for O2    a_1 =', a_1
         WRITE(numout,*) '    saturation coefficient for O2    a_2 =', a_2
         WRITE(numout,*) '    saturation coefficient for O2    a_3 =', a_3
         WRITE(numout,*) '    saturation coefficient for O2    a_4 =', a_4
         WRITE(numout,*) '    saturation coefficient for O2    a_5 =', a_5
         WRITE(numout,*) '    saturation coefficient for O2    b_0 =', b_0
         WRITE(numout,*) '    saturation coefficient for O2    b_1 =', b_1
         WRITE(numout,*) '    saturation coefficient for O2    b_2 =', b_2
         WRITE(numout,*) '    saturation coefficient for O2    b_3 =', b_3
         WRITE(numout,*) '    saturation coefficient for O2    c_0 =', c_0
         WRITE(numout,*) '       Schmidt coefficient for O2 a_1_o2 =', a_1_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_2_o2 =', a_2_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_3_o2 =', a_3_o2
         WRITE(numout,*) '       Schmidt coefficient for O2 a_4_o2 =', a_4_o2
      ENDIF

      ! namblingiron : Iron cycle parameters
      REWIND (numnatg)
      READ   (numnatg, namblingiron)

      kfe_inorg = kfe_inorg/86400.d0
      kfe_org   = kfe_org/86400.d0

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist namblingiron'
         WRITE(numout,*) '    Constant for iron binding with ligands  kfe_eq_lig_max   = ', kfe_eq_lig_max
         WRITE(numout,*) '    Minimum ligand strength                 kfe_eq_lig_min   = ', kfe_eq_lig_min
         WRITE(numout,*) '    Irradiance scaling cte for Kfe          kfe_eq_lig_irr   = ', kfe_eq_lig_irr
         WRITE(numout,*) '    Low-iron threshold for lig stability    kfe_eq_lig_femin = ', kfe_eq_lig_femin
         WRITE(numout,*) '    Global uniform iron ligand concentration     felig_bkg   = ', felig_bkg
         WRITE(numout,*) '    1.5-order iron scavenging                    kfe_inorg   = ', kfe_inorg
         WRITE(numout,*) '    Adsorption rate coefficient for detritus       kfe_org   = ', kfe_org
         WRITE(numout,*) '    Minimum [o2] for oxic remineralization         oxy_min   = ', oxy_min        
         WRITE(numout,*) '    Prevent oxygen from becoming negative    ln_prev_o2lt0   = ', ln_prev_o2lt0
      ENDIF

      ! namblingbudget : Spatially integrated values for tracer budget
      REWIND (numnatg)
      READ   (numnatg, namblingbudget)

      IF(lwp) THEN                         
         WRITE(numout,*)
         WRITE(numout,*) ' Namelist parameter for mass conservation checking BLING'
         WRITE(numout,*) '    Flag to check mass conservation of phosphate = ', ln_bling_mass
         WRITE(numout,*) '    Flag to check external fluxes of tracers     = ', ln_bling_ext
      ENDIF

      !
      IF( .NOT.lk_iomput .AND. ln_diatrc ) THEN
         !
         ! Namelist namblingdia
         ! -------------------
         DO jl = 1, jp_bling_2d
            WRITE(blidia2d(jl)%sname,'("2D_",I1)') jl                      ! short name
            WRITE(blidia2d(jl)%lname,'("2D DIAGNOSTIC NUMBER ",I2)') jl    ! long name
            blidia2d(jl)%units = ' '                                        ! units
         END DO
         !                                 ! 3D output arrays
         DO jl = 1, jp_bling_3d
           IF    (  jl <  10 ) THEN; WRITE(blidia3d(jl)%sname,'("3D_",I1)') jl     ! short name
           ELSE                    ; WRITE(blidia3d(jl)%sname,'("3D_",I2)') jl     ! short name
           ENDIF
           WRITE(blidia3d(jl)%lname,'("3D DIAGNOSTIC NUMBER ",I2)') jl    ! long name
           blidia3d(jl)%units = ' '                                        ! units
         END DO

         REWIND( numnatg )               ! read natrtd
         READ  ( numnatg, namblingdia )

         DO jl = 1, jp_bling_2d
            jn = jp_bling0_2d + jl - 1
            ctrc2d(jn) = blidia2d(jl)%sname
            ctrc2l(jn) = blidia2d(jl)%lname
            ctrc2u(jn) = blidia2d(jl)%units
         END DO

         DO jl = 1, jp_bling_3d
            jn = jp_bling0_3d + jl - 1
            ctrc3d(jn) = blidia3d(jl)%sname
            ctrc3l(jn) = blidia3d(jl)%lname
            ctrc3u(jn) = blidia3d(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : natadd'
            DO jl = 1, jp_bling_3d
               jn = jp_bling0_3d + jl - 1
               WRITE(numout,*) '  3d diag nb : ', jn, '    short name : ', ctrc3d(jn), &
                 &             '  long name  : ', ctrc3l(jn), '   unit : ', ctrc3u(jn)
            END DO
            WRITE(numout,*) ' '

            DO jl = 1, jp_bling_2d
               jn = jp_bling0_2d + jl - 1
               WRITE(numout,*) '  2d diag nb : ', jn, '    short name : ', ctrc2d(jn), &
                 &             '  long name  : ', ctrc2l(jn), '   unit : ', ctrc2u(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF

      IF ( .NOT.lk_iomput .AND. ln_diabio )  THEN
         !
         ! Namelist namblingdbi
         ! -------------------
         DO jl = 1, jp_bling_trd
            IF(     jl <  10 ) THEN   ;   WRITE (blidiabio(jl)%sname,'("BIO_",I1)') jl      ! short name
            ELSEIF (jl < 100 ) THEN   ;   WRITE (blidiabio(jl)%sname,'("BIO_",I2)') jl  
            ELSE                      ;   WRITE (blidiabio(jl)%sname,'("BIO_",I3)') jl
            ENDIF
            WRITE(blidiabio(jl)%lname,'("BIOLOGICAL TREND NUMBER ",I2)') jl                 ! long name
            blidiabio(jl)%units = 'mmolP/m2/s '                                            ! units
         END DO

         REWIND( numnatg )
         READ  ( numnatg, namblingdbi ) 
 
         DO jl = 1, jp_bling_trd
            jn = jp_bling0_trd + jl - 1
            ctrbio(jl) = blidiabio(jl)%sname
            ctrbil(jl) = blidiabio(jl)%lname
            ctrbiu(jl) = blidiabio(jl)%units
         END DO

         IF(lwp) THEN                   ! control print
            WRITE(numout,*)
            WRITE(numout,*) ' Namelist : namblingdbi'
            DO jl = 1, jp_bling_trd
               jn = jp_bling0_trd + jl - 1
print*, 'control...', jl, jn
               WRITE(numout,*) '  biological trend No : ', jn, '    short name : ', ctrbio(jn), &
                 &             '  long name  : ', ctrbil(jn), '   unit : ', ctrbiu(jn)
            END DO
            WRITE(numout,*) ' '
         ENDIF
         !
      ENDIF
                
   END  SUBROUTINE  trc_nam_bling
   
#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                             No BLINGv0
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nam_bling                      ! Empty routine
   END  SUBROUTINE  trc_nam_bling
#endif  

   !!======================================================================
END MODULE trcnam_bling
