      Subroutine  beaminit
      Implicit none
      include 'gconst.inc'          !geant
      include 'uggeom.inc'          !local
      include 'rescom.inc'
      include 'uevent.inc'           !local
      include 'res.inc'
      include 'beamcom.inc'


      REAL*4 beamz,tlif,ubuf(1),resz,etot
      REAL*4 recoilenerg,beammom,refmom  !In MeV
      REAL*4 gamma,beta,totmass,eint,excit,ereccm,toten,momm,betacm,
     +       gamcm,
     +       erec, trec, treco, eloss, e0rec, pcm, neutmass, treco0,
     +       Bref, Eref, cmag, Btun, Etun, Atun
      Integer*4 itrktyp, nubuf, i
      character*20 state
C
C===
C     Emittances from Laxtal area divided by 4 pi
C     delx, dely are 2 sigma values for spot size
C===
      REAL ex0/ 2.7e-4 /, ey0/ 2.7e-4 /, el0/ 5.0e-12 /
      REAL betagamma, betagamma0 / 0.01851 /
      REAL delx/ 0.25 /, dely/ 0.25 /
C
      INTEGER imate, partid, ixst
      REAL dedx, pcut(5) , beamm
C.
      REAL*4 etaref,etatune



C
C If the resonance is being placed in the center of the target...
C ex. (FFKEY: BEAM .eq. 0)
        If(beamenerg.eq.0) then


C First calculate nominal beam and recoil product energies at
C center of target
        CALL gfpart(80,state,itrktyp,beammass,beamz,tlif,ubuf,nubuf)
        CALL gfpart(81,state,itrktyp,resmass,resz,tlif,ubuf,nubuf)
        print*, '***', targmass, resmass, beammass, resenerg, prodm
        if(lkine.ne.13) then
           eres0 = (prodm/targmass)*.001 !! non-relativistic transformation from lab to CM 
        else
           eres0 = ((beammass+targmass)/targmass)*.001 !! non-relativistic transformation from lab to CM 1.213E-3
        endif
        print*, prodm,targmass,prodm/targmass,eres0
        e0beam = (resmass-beammass-targmass)*(resmass+beammass+targmass)
     &   /(2.0*targmass)
        print*, e0beam
        etot = e0beam+beammass
        e0recoil=(resmass**2+prodm**2)*(etot+targmass)/(2*resmass**2)
     &    - prodm   !ErecoilCM * gamma = value for 90deg CM gamma
C       e0recoil=5.6
C       print*, '***RECOIL ENERGY: ', e0recoil
!  This is the place any angular distribution equations would come in.
!  Lorentz boosts need to be made for each angle in 3-d kinematics. 
!
!  Now corrections for energy loss in the target gas
        imate  = mtarg
        partid = 80
C.
C.
C     Divide target into 100 slices, add up energyloss to get beam energy
c        eloss = 0.
c        beamenerg = e0beam
c        DO i = 1, 100
c           beamenerg = beamenerg + eloss
          
           CALL gftmat(imate,partid,'LOSS',1,e0beam,dedx,pcut,ixst)
c           eloss = (entdens/100.)*dedx * 0.001
           
c        ENDDO
c        beamenerg = beamenerg*1000.

        
C.==
C.=     entdens - The gas thickness to target center
C.=     beamenerg, dedx in MeV/cm , e0beam in GeV
C.==
        
c        If(beamenerg.ne.0.)then
c         continue    
c        Else !If(beamenerg.eq.0) then

        beamenerg = e0beam*1000. + dedx * entdens
c        Endif
        beammom=sqrt(beamenerg*(beamenerg+2000.*beammass))  !in Mev/c        
        beamvel = clight*beammom*.001/beammass
        gamma = 1.0/sqrt(1.0 - (beamvel/clight)**2)
        betagamma = beamvel/clight*gamma  !sign wrong on original?
C.
        write(6,*) '+++++++++++++++++BEAM AND TARGET+++++++++++++++'
        write(6,*) 'Beam energy       ',beamenerg,' MeV', ' Momentum ',
     &  beammom, ' MeV/c'
        write(6,*) 'Gas half thickness',entdens,' cm'
        write(6,*) 'dE/dx in target   ',dedx   ,' MeV/cm'
      write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        beamo = beamenerg - dedx*entdens*2.
        print*, 'Beam Mean Energy leaving target' , beamo, ' MeV'
C.
C.-->   Calculate initial beam distribution parameters
C.
        sigx = delx/2.
        sigy = dely/2.
C.
C.-->   assumed Gaussian buncht = 1 sigma
C.
        buncht = 1.E-9
        bunchl = buncht*beamvel
C.
C.-->   scale ex and ey
C.
        ex = ex0*betagamma0/betagamma
        ey = ey0*betagamma0/betagamma
        el = el0*betagamma0/betagamma
C.
        If(sigx .eq. 0.)then
          amax = 0.
          bmax = 0.
          emax = 0.
        Else
          amax = ex/sigx
          bmax = ey/sigy
          emax = el/buncht
        Endif
!  Now energy loss for the recoil
        partid = irecoil
C.
C.
        e0rec = e0recoil
        CALL gftmat(imate,partid,'LOSS',1,e0recoil,dedx,pcut,ixst)
        e0recoil = e0recoil -0.001*dedx*exitdens
        recoilenerg = e0recoil*1000.    !in MeV
        recoilenerg = recoilenerg*(1.+energscale)
        recoilmom = sqrt(recoilenerg*(recoilenerg+2000.*prodm))
        write(6,*) '++++++++++++++++RECOIL+++++++++++++++++++++++' 
        write(6,*)'Recoil Mean Energy from reaction',
     +            e0rec*1000.,' MeV'
        write(6,*)  'Gas half thickness',exitdens,' cm'
        write(6,*)   'dE/dx in target   ',dedx   ,' MeV/cm'
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        write(6,*) 'Recoil Mean Energy leaving target ',e0recoil*1000.,
     + ' MeV',   ' Momentum',recoilmom,' MeV/c' 
        write(6,*) 'Recoil energy tuned to ', recoilenerg,' MeV'
!
!  Determine scaling parameters for this reaction c/w the reference tune
!
        etaref=refenerg/(2*(refatno*amumev)**2)
        etatune=recoilenerg/(2*(prodm*1000.)**2)
        refmom = sqrt(refenerg*(refenerg+2*refatno*amumev))
        bscale = recoilmom/fkine(2)/ (refmom/refq)
        escale = recoilenerg*refq/(refenerg*fkine(2))*
     +     (1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref)
C. ----- Trick to make beam particles get through if needed 
C. ----- (rescale Electric Dipoles)
C.       escale = escale*prodm/beammass
C. ----- 
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        write(6,*) 'Magnetic element scale factor ',bscale
        write(6,*) 'Electric element scale factor ',escale
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'

        Else
           if(alpha)then
              print*, '--------- alpha source ----------'
            CALL gfpart(80,state,itrktyp,beammass,beamz,tlif,ubuf,nubuf)
              beammom=sqrt(beamenerg*(beamenerg+2000.*beammass))
              
              prodm = beammass
              recoilenerg = beamenerg*(1.+energscale)
              recoilmom = sqrt(recoilenerg*(recoilenerg+2000.*prodm))
              etaref=refenerg/(2*(refatno*amumev)**2)
              etatune=recoilenerg/(2*(prodm*1000.)**2)
              refmom = sqrt(refenerg*(refenerg+2*refatno*amumev))
              bscale = recoilmom/fkine(2)/ (refmom/refq)
              escale = recoilenerg*refq/(refenerg*fkine(2))*
     +             (1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref)
            write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
              write(6,*) 'Magnetic element scale factor ',bscale
              write(6,*) 'Electric element scale factor ',escale
            write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
            print*, 'beaminit.f'
         else                   !If(beamenerg.eq.0)
              print*, '********** Using BEAM FFCARD **********'
            CALL gfpart(80,state,itrktyp,beammass,beamz,tlif,ubuf,nubuf)
              CALL gfpart(81,state,itrktyp,resmass,resz,tlif,ubuf,nubuf)
         if(lkine.ne.13) then
            eres0 = prodm/targmass*0.001
         else
            eres0 = ((beammass+targmass)/targmass)*.001 !! non-relativistic transformation from lab to CM 1.213E-3
         endif
         imate = mtarg
         partid = 80
       CALL gftmat(imate,partid,'LOSS',1,beamenerg*.001,dedx,pcut,ixst)
         beammom=sqrt(beamenerg*(beamenerg+2000.*beammass))
         beamvel = clight*beammom*.001/beammass
         gamma = 1.0/sqrt(1.0-(beamvel/clight)**2)
         betagamma = beamvel/clight*gamma
         beamo = beamenerg - dedx*entdens*2.
         beamm = beamenerg - dedx*entdens
C.
         write(6,*) '+++++++++++++++BEAM AND TARGET+++++++++++++++++'
         write(6,*) 'Beam energy    ',beamenerg,' MeV',' Momentum ',
     &   beammom, ' MeV/c'
         write(6,*) 'Gas half thickness',entdens, ' cm'
         write(6,*) 'dE/dx in target  ',dedx,' MeV/cm'
         write(6,*) 'Beam energy at target exit ',beamo,' MeV'
         write(6,*) 'Beam energy at target centre ',beamm,' MeV'
         write(6,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
C.
C.       Initial beam distribution parameters
C.       
         sigx = delx/2.
         sigy = dely/2.
C.
C.-->   assumed Gaussian buncht = 1 sigma
C.
        buncht = 1.E-9
        bunchl = buncht*beamvel
C.
C.-->   scale ex and ey
C.
        ex = ex0*betagamma0/betagamma
        ey = ey0*betagamma0/betagamma
        el = el0*betagamma0/betagamma
C.
        If(sigx .eq. 0.)then
          amax = 0.
          bmax = 0.
          emax = 0.
        Else
          amax = ex/sigx
          bmax = ey/sigy
          emax = el/buncht
        Endif

C.
C.--> Recoil energy at center of target
C     .double pcm = sqrt((pow(S - M3*M3 - M4*M4, 2) - 4*M3*M3*M4*M4) / (4*S));
        
        totmass = sqrt( (beammass+targmass)**2 +
     +     2.*targmass*(beamm/1000.))

        eint = totmass - beammass - targmass
C       different excit for (a,n)
        If(lkine.eq.13) then
           neutmass = 0.9395654
           excit = (beammass+targmass-prodm-neutmass) + eint
        Else
           excit = (beammass+targmass-prodm) + eint
        Endif
        ereccm = sqrt( prodm**2 + excit**2 )
        toten = (beamm/1000.) + beammass
        momm = sqrt( toten**2 - beammass**2 )
        betacm = momm/(toten+targmass)
        gamcm = (toten+targmass)/sqrt((toten+targmass)**2-momm**2)
        erec = gamcm*ereccm
        trec = erec - prodm
C        print*, "totmass, beammass, targmass, eint, excit, ereccm"
C        print*, "momm, betacm, gamcm, erec, trec"
C        print*, totmass, beammass, targmass, eint, excit, ereccm
C        print*, momm, betacm, gamcm, erec, trec
        partid = irecoil
        CALL gftmat(imate,partid,'LOSS',1,trec,dedx,pcut,ixst)
        treco = trec - 0.001*dedx*exitdens
		treco0 = treco
        treco = treco*(1.+energscale)
        recoilmom = 1000*(sqrt((treco + prodm)**2 - prodm**2))
        write(6,*) '++++++++++++++++RECOIL+++++++++++++++++++++++' 
        write(6,*) 'CM, Lab resonance energy', resenerg,
     +            (resenerg*(beammass+targmass)/targmass)
        write(6,*)'Recoil Mean Energy at target centre (90 deg gamma)',
     +            trec*1000.,' MeV'
        write(6,*)  'Gas half thickness',exitdens,' cm'
        write(6,*)   'dE/dx in target   ',dedx   ,' MeV/cm'
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        write(6,*)   'Recoil Mean Energy leaving target [90 deg CM]'
     + , treco0*1000., ' MeV',   ' Momentum',recoilmom,' MeV/c'
        recoilenerg = treco*1000.
		write(6,*) 'Energy mistune (%) ', energscale*100
        write(6,*) 'Recoil energy tuned to ', recoilenerg,' MeV'

!
!  Determine scaling parameters for this reaction c/w the reference tune
!
        write(6,*) ' Initial (b,e)scale:', bscale, escale
        etaref=refenerg/(2*(refatno*amumev)**2)
        etatune=recoilenerg/(2*(prodm*1000.)**2)
        refmom = sqrt(refenerg*(refenerg+2*refatno*amumev))
        bscale = recoilmom/fkine(2)/ (refmom/refq)
        escale = recoilenerg*refq/(refenerg*fkine(2))*
     +     (1+etatune)/(1+2*etatune)*(1+2*etaref)/(1+etaref)
C
C  Calculate B field (hall probe) and E field (setpoint) for
C  reference tune, E=1.8885, A=19, q=4
C  B from DAH 2012 NIM paper
C  E from printout in DRAGON counting room
        cmag = 48.23 !design MeV/T^2
        Bref = 
     +  sqrt((refenerg/refatno + (refenerg/refatno)**2/(2*amumev))
     +       * (refatno/refq)**2 / cmag) * 1E4
        Eref = 2468.*(Bref/1E4)**2/(refatno/refq)
C
C  Calculate B field (hall probe) and E field (setpoint) for
C  scaled recoil tune
        Atun = (prodm*1000)/amumev
        Btun = 
     +  sqrt((recoilenerg/Atun + (recoilenerg/Atun)**2/(2*amumev))
     +       * (Atun/refq)**2 / cmag) * 1E4
        Etun = 2468.*(Btun/1E4)**2/(Atun/fkine(2))
        
        
C.----- Trick to make beam particles get through if needed 
C.  ----- (rescale Electric Dipoles)
C.        escale = escale*prodm/beammass
C. ----- 
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        write(6,*) 'Reference E, A, Q, eta, P, B, E',
     +       refenerg, refatno, refq, etaref, refmom, Bref, Eref
        write(6,*) 'Scaled E, A, Q, eta, P, B, E',
     +   recoilenerg, Atun, fkine(2), etatune, recoilmom,
     +   Btun, Etun
        write(6,*) 'Magnetic element scale factor ',bscale
        write(6,*) 'Bref * bscale ', Bref*bscale
        write(6,*) 'Electric element scale factor ',escale
        write(6,*) 'Eref * escale ', Eref*escale
        write(6,*)'++++++++++++++++++++++++++++++++++++++++++++++++'
        Endif
        endif
        print*, 'beaminit,f'
        return
        end

























