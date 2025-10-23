C
      SUBROUTINE ureact
C
      IMPLICIT none
C
      include 'gckine.inc'          !geant
C
      include 'uevent.inc'          !local
      include 'rescom.inc'          !local
      include 'beamcom.inc'
      include 'res.inc'             !local
      include 'params.inc'          !local
C
      INTEGER ireaction
C
      INTEGER nubuf, i, j
      PARAMETER ( nubuf = 1)
      REAL    ubuf(nubuf)
C
      REAL hbar
      REAL amugev, hmass, hemass, deutmass, c12mass, neutmass
C
      PARAMETER ( hbar   = 6.582122   E-22 )
      PARAMETER ( amugev = 0.93149432 E+00 )
      PARAMETER ( hmass  = 1.007825032 * amugev)
      PARAMETER ( hemass = 4.002603250 * amugev)
      PARAMETER (deutmass = 2.0141018 * amugev )
      PARAMETER (neutmass = 1.0086649 * amugev )
      PARAMETER ( c12mass = 12.0 * amugev )
C
      REAL devmass,aamass,aamass1,zbeam, elevel, aprod,  tlif
C
      INTEGER mode(10)
      REAL    brat(10)
      CHARACTER*2 num(31)
      CHARACTER*120 CARDNAME
      CHARACTER*5 INPUT
      CHARACTER*3 targtyp
      LOGICAL fexist
      REAL brsum
      DATA num/ '1 ','2 ','3 ','4 ','5 ','6 ','7 ','8 ','9 ','10', 
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24','25','26','27','28','29','30','31'/


      
      
      REAL rndm(2)
      INTEGER ilevel


      INTEGER MAXBR, NLEV
      PARAMETER (MAXBR=6, NLEV=18)

      INTEGER St43Ca(0:NLEV-1, MAXBR)
      REAL    Br43Ca(0:NLEV-1, MAXBR)
      REAL    Ex43Ca(0:NLEV-1)

C--   initialize padding
      INTEGER LL, JJ
      DO LL = 0, NLEV-1
         DO JJ = 1, MAXBR
            St43Ca(LL,JJ) = -1
            Br43Ca(LL,JJ) = 0.0
         END DO
      END DO

C-- level energies (keV)
      DATA Ex43Ca /
     &   0.000,   372.762,  593.394,  990.257, 1394.473,
     &   1677.840, 1901.990, 1931.530, 1957.400, 2046.210,
     &   2067.210, 2093.810, 2102.700, 2223.900, 2248.000,
     &   2249.010, 2272.800, 2409.680 / !, 2523.000 /

      DO JJ=1,NLEV-1 
         Ex43Ca(JJ) = Ex43Ca(JJ)/1E3
      ENDDO
      
C==================================================================
C     Gamma decay branches for 43Ca (levels 0–18)
C==================================================================

C--   level 0 : ground state (no gammas)

C--   level 1 : 372.8 keV
      St43Ca(1,1) = 0
      Br43Ca(1,1) = 1.0

C--   level 2 : 593.4 keV
      St43Ca(2,1) = 1           ! → 372.8
      Br43Ca(2,1) = 0.423
      St43Ca(2,2) = 0           ! → g.s.
      Br43Ca(2,2) = 1.0

C--   level 3 : 990.3 keV
      St43Ca(3,1) = 1           ! → 372.8
      Br43Ca(3,1) = 1.0
      St43Ca(3,2) = 2           ! → 593.4
      Br43Ca(3,2) = 0.1493
      St43Ca(3,3) = 0           ! → g.s.
      Br43Ca(3,3) = 0.0036

C--   level 4 : 1394.5 keV
      St43Ca(4,1) = 1           ! → 372.8
      Br43Ca(4,1) = 1.0
      St43Ca(4,2) = 3           ! → 990
      Br43Ca(4,2) = 0.187
      St43Ca(4,3) = 2           ! → 593
      Br43Ca(4,3) = 0.075
      St43Ca(4,4) = 0           ! → g.s.
      Br43Ca(4,4) = 0.09

C--   level 5 : 1677.8 keV
      St43Ca(5,1) = 0
      Br43Ca(5,1) = 1.0

C--   level 6 : 1902.0 keV
      St43Ca(6,1) = 0
      Br43Ca(6,1) = 1.0
      St43Ca(6,2) = 4           ! → 1394
      Br43Ca(6,2) = 0.24
      St43Ca(6,3) = 3           ! → 990
      Br43Ca(6,3) = 0.19

C--   level 7 : 1931.5 keV
      St43Ca(7,1) = 0
      Br43Ca(7,1) = 1.0
      St43Ca(7,2) = 1           ! → 372.8
      Br43Ca(7,2) = 0.59
      St43Ca(7,3) = 2           ! → 593
      Br43Ca(7,3) = 0.117

C--   level 8 : 1957.4 keV
      St43Ca(8,1) = 0           ! → 593
      Br43Ca(8,1) = 1.0
      St43Ca(8,2) = 3           ! → 990
      Br43Ca(8,2) = 0.28

C--   level 9 : 2046.2 keV
      St43Ca(9,1) = 0
      Br43Ca(9,1) = 1.0
      St43Ca(9,2) = 1           ! → 372.8
      Br43Ca(9,2) = 0.32
      St43Ca(9,3) = 2           ! → 593
      Br43Ca(9,3) = 0.13
      St43Ca(9,4) = 3           ! → 990
      Br43Ca(9,4) = 0.11
      St43Ca(9,5) = 4           ! → 1394
      Br43Ca(9,5) = 0.02

C--   level 10 : 2067.2 keV
      St43Ca(10,1) = 0
      Br43Ca(10,1) = 1.0
      St43Ca(10,2) = 1
      Br43Ca(10,2) = 0.28

C--   level 11 : 2093.8 keV
      St43Ca(11,1) = 0
      Br43Ca(11,1) = 1.0

C--   level 12 : 2102.7 keV
      St43Ca(12,1) = 1          ! → 372.8
      Br43Ca(12,1) = 1.0
      St43Ca(12,2) = 2          ! → 593
      Br43Ca(12,2) = 0.50
      St43Ca(12,3) = 0          ! → g.s.
      Br43Ca(12,3) = 0.50

C--   level 13 : 2223.9 keV
      St43Ca(13,1) = 2          ! → 593
      Br43Ca(13,1) = 1.0
      St43Ca(13,2) = 1          ! → 372.8
      Br43Ca(13,2) = 0.745

C--   level 14 : 2248.0 keV
C     (no adopted gammas → remains padded)
      St43Ca(14,1) = 0
      Br43Ca(14,1) = 1.0

C--   level 15 : 2249.0 keV
      St43Ca(15,1) = 0
      Br43Ca(15,1) = 1.0
      St43Ca(15,2) = 1
      Br43Ca(15,2) = 0.126
      St43Ca(15,3) = 5          ! → 1677
      Br43Ca(15,3) = 0.023

C--   level 16 : 2272.8 keV
      St43Ca(16,1) = 3          ! → 990
      Br43Ca(16,1) = 1.0
      St43Ca(16,2) = 4          ! → 1394
      Br43Ca(16,2) = 0.19

C--   level 17 : 2409.7 keV
      St43Ca(17,1) = 6    ! → 1901.99 keV
      Br43Ca(17,1) = 0.24
      St43Ca(17,2) = 5    ! → 1677.84 keV
      Br43Ca(17,2) = 0.15
      St43Ca(17,3) = 4    ! → 1394.47 keV
      Br43Ca(17,3) = 0.98
      St43Ca(17,4) = 0    ! → g.s.
      Br43Ca(17,4) = 1.0

C--   level 18 : 2523.0 keV
C     (no adopted gammas → remains padded)

C     Convert to branching fractions in per-cent
C     (sum to 100%)
      DO LL = 0, NLEV-1
         brsum = 0.
         DO JJ = 1, MAXBR
            brsum = brsum + Br43Ca(LL,JJ)
         ENDDO
         DO JJ = 1, MAXBR
            IF(St43Ca(LL,JJ).NE.-1)THEN
               Br43Ca(LL,JJ) = 100. * Br43Ca(LL,JJ)/brsum
            ENDIF
         ENDDO
      ENDDO

      print*, "***43CA STATES & BRANCHES***"
      DO LL = 0, NLEV-1
         Print*, "State --",LL,Ex43Ca(LL)," MeV"
         DO JJ = 1, MAXBR
            IF(St43Ca(LL,JJ).NE.-1)THEN
               PRINT*,"---->",St43Ca(LL,JJ),",",
     +              Ex43Ca(JJ)," MeV:",Br43Ca(LL,JJ)
            ENDIF
         ENDDO
      ENDDO


c$$$
c$$$      
c$$$      REAL Ex43Ca(20)
c$$$      DATA Ex43Ca / 
c$$$     &     0.000000, 0.372762, 0.593394, 0.990257, 1.394473,
c$$$     &     1.677840, 1.901990, 1.931530, 1.957400, 2.046210,
c$$$     &     2.067210, 2.093810, 2.102700, 2.223900, 2.248000,
c$$$     &     2.249010, 2.272800, 2.409680, 2.523000, -1.00000
c$$$     &     /
c$$$
c$$$      REAL Br43Ca(20,5)
c$$$      DATA Br43Ca /
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! g.s.
c$$$     &     1.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 1
c$$$     &     0.423, 1.0, 0.0, 0.0, 5.0  ! e.x. 2
c$$$     &     1.0, 1.0, 0.1493, 0.0, 0.0, ! e.x. 3
c$$$     &     0.09, 1.0, 0.075, 0.187, 0.0, ! e.x. 4
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 5
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 6
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 7
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 8
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 9
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 10
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 11
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 12
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 13
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 14
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 15
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 16
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 17
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 18
c$$$     &     0.0, 0.0, 0.0, 0.0, 0.0, ! e.x. 19 (empty)
c$$$     &     /
c$$$      REAL St43Ca(20,5)         ! state to decay to
c$$$      DATA St43Ca /
c$$$     &     -1, -1, -1, -1, -1, ! g.s.
c$$$     &     0, -1, -1, -1, -1, ! e.x. 1
c$$$     &     1, 0, -1, -1, -1  ! e.x. 2
c$$$     &     0, 1, 2, -1, -1, ! e.x. 3
c$$$     &     0, 1, 2,  3, -1, ! e.x. 4
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 5
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 6
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 7
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 8
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 9
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 10
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 11
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 12
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 13
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 14
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 15
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 16
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 17
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 18
c$$$     &     -1, -1, -1, -1, -1, ! e.x. 19 (empty)
c$$$     &     /
      
C
C======================================================================C
C                                                                      C
C      Initial beam and reaction information passed via FKIN card      C
C                                                                      C
C      LKINE reaction number                                           C
C                                                                      C
C      ( 1) 13N(p,g)14O                                                C
C      ( 2) 15O(alpha,g)19Ne                                           C
C      ( 3) 25Al(p,g)26Si                                              C
C      ( 4) 17F(p,g)18Ne                                               C
C      ( 5) 18F(p,g)19Ne                                               C
C      ( 6) 19Ne(p,g)20Na                                              C
C      ( 7) 20Na(p,g)21Mg 532 3/2-                                     C
C      ( 8) 21Na(p,g)22Mg(220)                                             C
C      ( 9) 23Mg(p,g)24Al                                              C
C      (10) 26mAl(p,g)27Si                                             C
C      (11) 7Be(p,g)8B                                                 C
C      (12) 21Na(d,n)22Mg
C      (13) 22Na(a,n)25Mg 11.83(gs)                                    C
C      (14) 23Na(p,g)24Mg 2+                                           C
C      (15) 20Ne(p,g)21Na nonres
C      (16) 20Na(p,g)21Mg131                                           C
C      (17) 22Ne(alpha,g)26Mg
C      (18) 21Na(p,g)22Mg(825)
C      (19) 12C(a,g)16O 
C======================================================================C
C
C-- MOD 10/06/03 C.Ruiz
C-- Gamma cascades (isotropic emission) are included for
C-- product nuclei through up to 32 states (including final 
C-- ground state).
C--
C-- Each product energy level is defined as a GEANT particle.
C-- Mass in GeV, lifetime and branching ratios to the lower states
C-- are specified for each level.
C--
C-- Within this routine the array  BRAT and MODE   specify
C-- the decay branches. MODE = part1 + 100*part2
C-- where part1 is the GEANT particle number correspnding to the
C-- energy level or gamma (=1).
C--
C--                   LABEL     IDPART
C--  resonance                   81
C--                              82
C--                              83
C--                              84
C--                              85
C--                              86
C--                              .
C--                              .
C--  ground state              irecoil < 100
C--
C=======================================================================
C
      ilight = 1 ! default to radiative capture, light particle is gamma
      ireaction = abs(lkine)
C
      Select case (ireaction)
C
      case(0)
C
        STOP
C
      case(1)
C
C       ' (1) 13N(p,g)14O '
C
        zbeam =  7.
        abeam = 13.
        atarg =  1.
        aprod = atarg + abeam
C
        resenerg = 0.526
        reswidth = 0.000037
C
        elevel = 0.0
C
        write(6,*)'|**** 13N(p,g)14O reaction NOT implemented yet ****|'
        write(6,*)'Resonance energy ', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(2)
C
C       ' (2) 15O(alpha,g)19Ne '
C
        zbeam =  8.
        abeam = 15.
        atarg =  4.
        zprod = 10.
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
        devmass = 2855.4E-6
        aamass  = abeam*amugev + devmass
        tlif    = 122.2
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'O15',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        resenerg = 0.5036
        reswidth = 1.E-11
        aamass   = aamass + resenerg/1000. + hemass
        tlif     = hbar/reswidth
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Ne19_3/2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Ne19 gamma decays from resonance
C
        devmass = 1751.1E-6
        prodm   = aprod*amugev + devmass
C
        elevel  = 1.536
        aamass  = prodm + elevel/1000.
        tlif    = 2.8E-11
C
        CALL gspart(82,'3_Ne19_3/2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = .275
        aamass = prodm + elevel/1000.
        tlif   = 6.3E-11
C
        CALL gspart(83,'2_Ne19_1/2',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel =.238
        aamass = prodm + elevel/1000.
        tlif   = 2.6E-8
C
        CALL gspart(84,'1_Ne19_5/2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 85
        irecoil = 85
C
        tlif = 1000.
C
        CALL gspart(85,'Ne19_1/2+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1) = 80.
        mode(1) = 1 + 100*85
        brat(2) = 15.
        mode(2) = 1 + 100*83
        brat(3) =  5.
        mode(3) = 1 + 100*82
C
        CALL uzero(brat,4,6)
        CALL uzero(mode,4,6)
C
        CALL gsdk(81,brat,mode)
C
C--     3/2+ state
C
        brat(1) = 95.
        mode(1) = 1 + 100*84
        brat(2) =  5.
        mode(2) = 1 + 100*83
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(82,brat,mode)
C
        brat(1) = 100.
        mode(1) = 1 + 100*95
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(83,brat,mode)
        CALL gsdk(84,brat,mode)
C
        write(6,*)'|**** 15O(alpha,gamma)19Ne reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
      case(3)
C
C       ' (3) 25Al(p,g)26Si '
C
        zbeam = 13.
        abeam = 25.
        atarg =  1.
        zprod = 14.
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
        devmass = -8915.7E-6
        aamass  = abeam*amugev + devmass
        tlif    = 7.18
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Al25',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        resenerg = 0.452   !calculated state...see C. Illiadis Phys
C                           Rev C53(1)(1995)475
        reswidth = 0.00006
        aamass   = aamass + resenerg/1000. + hmass
        tlif     = hbar/reswidth
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Si26_3+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Si26 gamma decays from resonance
C
        devmass = -7145.E-6
        prodm   = aprod*amugev + devmass
C
        elevel = 4.183
        aamass = prodm + elevel/1000.
        tlif    = 150.E-15
C
        CALL gspart(82,'4_Si26_3+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 3.756
        aamass = prodm + elevel/1000.
        tlif   = 700.E-15
C
        CALL gspart(83,'3_Si26_3+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 2.783
        aamass = prodm + elevel/1000.
        tlif   = 210.E-15
C
        CALL gspart(84,'2_Si26_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 1.796
        aamass = prodm + elevel/1000.
        tlif   = 620.E-15
C
        CALL gspart(85,'1_Si26_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 86
        irecoil=86
C
        tlif = 1000.
C
        CALL gspart(96,'Si26_0+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1) = 87.
        mode(1) = 1 + 100*82
        brat(2) = 8.
        mode(2) = 1 + 100*83
        brat(3) = 5.
        mode(3) = 1 + 100*85
C
        CALL uzero(brat,4,6)
        CALL uzero(mode,4,6)
C
        CALL gsdk(81,brat,mode)
C
        brat(1) = 47.
        mode(1) = 1 + 100*84
        brat(2) = 53.
        mode(2) = 1 + 100*85
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(82,brat,mode)
C
        brat(1) = 70.
        mode(1) = 1 + 100*85
        brat(2) = 30.
        mode(2) = 1 + 100*84
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(83,brat,mode)
C
        brat(1) = 69.
        mode(1) = 1 + 100*85
        brat(2) = 31.
        mode(2) = 1 + 100*86
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(84,brat,mode)
C
        brat(1) = 100.
        mode(1) = 1 + 100*86
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(85,brat,mode)
C
        write(6,*)'|**** 25Al(p,g)26Si  reaction ****|'
        write(6,*)'Resonance energy',resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
      case(4)
C
C       ' (4) 17F(p,g)18Ne '
C
        zbeam =  9.
        abeam = 17.
        atarg =  1.
        aprod = abeam +atarg
C
        resenerg = 0.64
        reswidth = 0.
C
        elevel = 2.67
C
        write(6,*)'|**** 17F(p,g)18Ne reaction NOT implemented ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(5)
C
C       ' (5) 18F(p,g)19Ne '
C
        zbeam =  9.
        abeam = 18.
        atarg =  1.
C
        resenerg = 0.3308
        reswidth = 0.0
C
        elevel = 6.742
C
        write(6,*)'|**** 18F(p,g)19Ne reaction NOT implemented ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(6)
C
C       ' (6) 19Ne(p,g)20Na '
C
        zbeam = 10.
        abeam = 19.
        atarg =  1.
        zprod = 11.
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
        devmass = 1751.E-6
        aamass  = abeam*amugev + devmass
        tlif    = 1.
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Ne19',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        resenerg = 0.451
        reswidth = 7.E-9       !w=1 don't know spins
        aamass   = aamass + resenerg/1000. + hmass
        tlif     = hbar/reswidth
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Na20',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Na20 gamma decays from resonance
C
        elevel = 2.660
	aamass = aamass - elevel/1000.
C
C--     ground state --> idpart = 82
C
        irecoil = 82
        tlif = 1.                        ! made up
C
        CALL gspart(82,'gs_Na20',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1) = 100.
	mode(1) = 1 + 100*82
C
        CALL uzero(brat,2,6)
	CALL uzero(mode,2,6)
C
	CALL gsdk(81,brat,mode)
C
        write(6,*)'|**** 19Ne(p,gamma)20Na reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
      case(7)
C
C       ' (7) 20Na(p,g)21Mg536 '
C
        zbeam = 11.
        abeam = 20.
        atarg =  1.
        zprod = 12
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
       	devmass = 6845E-6
       	aamass  = abeam*amugev + devmass
        tlif    = 0.03
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Na20',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        resenerg = 0.536
        reswidth = hbar/tlif
        aamass   = aamass + resenerg/1000. + hmass
        ubuf(1)  = fkine(2)
        tlif   = 4.0E-14
C
        CALL gspart(81,'res_Mg21_536',8,aamass,zprod,tlif,ubuf,nubuf)
C    Define states for Mg21 gamma decays from resonance state
        devmass = 10912E-6
        prodm   = aprod*amugev + devmass

C
        elevel = 0.0
        amass = prodm +elevel/1000.
        tlif    = 1000.
        irecoil = 82
C
        CALL gspart(82,'Mg21_gs',8,aamass,zprod,tlif,ubuf,nubuf)

C--     branch info -- resonance decays
C
        brat(1) = 100.
        mode(1) = 1 + 100*82
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(81,brat,mode)
C
C
C
        write(6,*)'|**** 20Na(p,g)21Mg reaction  ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        
C
      case(8)
C
C       ' (8) 21Na(p,g)22Mg(220) '
C
        zbeam = 11.
        abeam = 21.
        atarg =  1.
        zprod = 12
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
       	devmass = -2184.3E-6
       	aamass  = abeam*amugev + devmass
        tlif    = 0.03
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Na21',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        tlif     = 4.E-14
        resenerg = 0.2124
        reswidth = hbar/tlif
        aamass   = aamass + resenerg/1000. + hmass
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Mg22_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Mg22 gamma decays from resonance
C
        devmass = -396.8E-6
        prodm   = aprod*amugev + devmass
C
        elevel  = 1.246
        aamass  = prodm + elevel/1000.
        tlif    = 3.E-11
C
        CALL gspart(82,'Mg22_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 83
C
        irecoil=83
        tlif   = 1000.
C
        CALL gspart(83,'Mg22_0+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1) = 87.
        mode(1) = 1 + 100*82
        brat(2) = 13.
        mode(2) = 1 + 100*83
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(81,brat,mode)
C
        brat(1) = 100.
        mode(1) = 1 + 100*83
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(82,brat,mode)
C
        write(6,*)'|**** 21Na(p,g)22Mg reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
      case(9)
C
C       ' (9) 23Mg(p,g)24Al '
C
        zbeam = 12.
        abeam = 23.
        atarg =  1.
C
        resenerg = 0.51
        reswidth = 0.0
C
        elevel = 2.38
C
        write(6,*)'|**** 23Mg(p,g)24Al reaction NOT implemented ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(10)
C
C       ' (10) 26mAl(p,g)27Si '
C
        zbeam = 13.
        abeam = 26.
        atarg =  1.
C
        resenerg = 0.201
        reswidth = 0.
C
        elevel = 7.893
C
        write(6,*)'|**** 26mAl(p,g)27Si reaction NOT implemented ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(11)
C
C       ' (11) 7Be(p,g)8B '
C
        zbeam = 4.
        abeam = 7.
        atarg = 1.
C
        resenerg = 0.2
        reswidth = 0.
C
        elevel = 0.338
C
        write(6,*)'|**** 7Be(p,g)8B reaction NOT implemented ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
        STOP
C
      case(12)
C
C       ' (12) 21Na(d,n)22Mg '
C
        zbeam = 11.
        abeam = 21.
        atarg =  2.
        zprod = 12
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
       	devmass = -2184.3E-6
       	aamass  = abeam*amugev + devmass
        tlif    = 0.03
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Na21',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     nonresonant level populated --> idpart = 81
C
        tlif     = 1e-20  !10 KeV width
        resenerg = 0.00
        reswidth = hbar/tlif
        aamass = sqrt((aamass+deutmass)**2 + 2.*deutmass*beamenerg*.001)
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'nonres_Mg23_',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Mg22  from neutron decays from resonance
C
        devmass = -396.8E-6
        prodm   = (aprod-1)*amugev + devmass !gs of 22Mg
C

        elevel  = 4.401

        aamass  = prodm + elevel/1000.
        tlif    = 3.E-11
C
        CALL gspart(82,'Mg22_2+_2',8,aamass,zprod,tlif,ubuf,nubuf)
        elevel  = 1.246

        aamass  = prodm + elevel/1000.
        tlif    = 3.E-11
C
        CALL gspart(83,'Mg22_2+_1',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 84
C
        irecoil=84
        tlif   = 1000.
C
        CALL gspart(84,'Mg22_0+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1)= 100.
        mode(1)= 13+100*82
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
        CALL gsdk(81,brat,mode)
C
C
        brat(1) = 87.
        mode(1) = 1 + 100*83
        brat(2) = 13.
        mode(2) = 1 + 100*84

C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(82,brat,mode)
C
        brat(1) = 100.
        mode(1) = 1 + 100*84
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(83,brat,mode)
C
        write(6,*)'|**** 21Na(d,ng)22Mg reaction ****|'
        write(6,*)' 100% to gs  + neutron'
      case(13)
C
C       ' (13) 22Ne(a,n)25Mg (gs) '
C
        zbeam = 10.
        abeam = 22.
        atarg =  4.
        zprod = 12
        aprod = atarg + abeam
        targmass = hemass
        ilight = 13             !neutron emission
C
C--     create beam particle --> idpart = 80
C
       	devmass = -8024.7202E-06
       	aamass  = abeam*amugev + devmass
        tlif    = 1000.
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Ne22',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        tlif     = hbar/0.931e-3 ! Jaeger 1.1keV lab->0.931 CM
        resenerg = 1.213
        reswidth = hbar/tlif
        aamass   = aamass + resenerg/1000. + hemass
C        aamass = sqrt((aamass+targmass)**2 + 2.*targmass*beamenerg*.001)
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'nonres_Mg26_',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Mg25  from neutron decays from resonance
C
        devmass = -13192.785E-6
        prodm   = (aprod-1)*amugev + devmass !gs of 25Mg
        write(6,*) "25Mg MASS: ", prodm
C

c$$$        elevel  = 0.
c$$$
c$$$        aamass1  = prodm + elevel/1000.
c$$$        if(aamass .lt. aamass1) 
c$$$     &  write(6,*) "Beam energy below 2+ threshold"
c$$$        tlif    = 1000
c$$$C
c$$$        CALL gspart(82,'Mg24_2+_2',8,aamass1,zprod,tlif,ubuf,nubuf)
c$$$C--     ground state --> idpart = 93
c$$$C

        irecoil=82
        tlif   = 1000.
C
        CALL gspart(82,'Mg25_gs',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
C     particle 13 is the neutron
        brat(1)= 100.
        mode(1)= 13+100*82
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
        CALL gsdk(81,brat,mode)
C
C

c$$$

C        brat(1) = 100.
c$$$        mode(1) = 1 + 100*82
c$$$
c$$$C
c$$$        CALL uzero(brat,2,6)
c$$$        CALL uzero(mode,2,6)
c$$$C
c$$$        CALL gsdk(82,brat,mode)
C
C
        write(6,*)'|**** 22Ne(a,n)25Mg reaction ****|'
        write(6,*)' 100% to gs'
C
      case(14)
C
C       ' (14) 23Na(p,g)24Mg '
C
        zbeam = 11.
        abeam = 23.
        atarg =  1.
        zprod = 12
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 90
C
       	devmass = -9530.0E-6
       	aamass  = abeam*amugev + devmass
        tlif    = 0.03
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Na23',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 91
C
        tlif     = 4.E-14
        elevel = 1.368
        devmass = -13930.7E-6
        prodm   = aprod*amugev + devmass
        resenerg = (prodm + elevel/1000. -aamass -hmass)*1000.
        if (resenerg .le.  0.) then
           resenerg = 0.0
           elevel = (aamass + hmass -prodm)*1000.
           write(6,*) " resonance energy negative - set to 0!"
        end if
        aamass = aamass+ hmass + resenerg/1000.
        reswidth = hbar/tlif
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Mg24_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Mg24 gamma decays from resonance
C
C
C
C--     ground state --> idpart = 82
C
        irecoil=82
        tlif   = 1000.
C
        CALL gspart(82,'Mg24_0+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
C
        brat(1) = 100.
        mode(1) = 1 + 100*82
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(81,brat,mode)
C
        write(6,*)'|**** 23Na(p,g)24Mg reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
      case(18)
C
C       ' (18) 21Na(p,g)22Mg(822) '
C
        zbeam = 11.
        abeam = 21.
        atarg =  1.
        zprod = 12
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
       	devmass = -2184.3E-6
       	aamass  = abeam*amugev + devmass
        tlif    = 0.03
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'Na21',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        tlif     = 84.0e-19
        resenerg = 0.8225
        reswidth = hbar/tlif/2 !This is gamma/2 =HWHM for BW
        aamass   = aamass + resenerg/1000. + hmass
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_Mg22_822',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     Define states in Mg22 gamma decays from resonance
C
        devmass = -396.8E-6
        prodm   = aprod*amugev + devmass
C
        elevel = (aamass -prodm)*1000.
        write(6,*)'|**** 21Na(p,g)22Mg(822) reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth*2,' MeV ', 'level ', elevel,' MeV'
C
        elevel  = 4.401
        aamass  = prodm + elevel/1000.
        tlif    = 3.E-14
C
        CALL gspart(82,'Mg22_4401+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C
        elevel  = 3.308
        aamass  = prodm + elevel/1000.
        tlif    = 3.E-14
C
        CALL gspart(83,'Mg22_3308',8,aamass,zprod,tlif,ubuf,nubuf)
C
C
        elevel  = 1.246
        aamass  = prodm + elevel/1000.
        tlif    = 3.E-11
C
        CALL gspart(84,'Mg22_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 85
C
        irecoil=85
        tlif   = 1000.
C
        CALL gspart(85,'Mg22_0+',8,prodm,zprod,tlif,ubuf,nubuf)
C
C--     branch info -- resonance decays
        brat(1) = 50.
        mode(1) = 1+100*84
        brat(2) = 50.
        mode(2) = 1 + 100*85
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(81,brat,mode)

C
        brat(1) = 87.
        mode(1) = 1 + 100*84
        brat(2) = 13.
        mode(2) = 1 + 100*85
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
C
        CALL gsdk(82,brat,mode)
C
        brat(1) = 100.
        mode(1) = 1 + 100*84
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(83,brat,mode)
C        brat(1) = 100.
        mode(1) = 1 + 100*85
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
C
        CALL gsdk(84,brat,mode)
C
C
      case(19)
C
C       ' (19) 12C(alpha,g)16O '
C
        zbeam =  6.
        abeam = 12.
        atarg =  4.
        zprod = 8.
        aprod = atarg + abeam
C
C--     create beam particle --> idpart = 80
C
        devmass = 0.0E-6
        aamass  = abeam*amugev + devmass
        tlif    = 1000.
        ubuf(1) = fkine(1)
C
        CALL gspart(80,'C12gs',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C--     resonant level populated --> idpart = 81
C
        resenerg = 4.358
        reswidth = 7.E-05
        aamass   = aamass + resenerg/1000. + hemass
        print*, aamass
        tlif     = hbar/reswidth
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_O16_2+',8,aamass,zprod,tlif,ubuf,nubuf)

        elevel = 11.520
        print*, 'Resonant mass: ', aamass
        rmass = aamass

        write(6,*)'|**** 12C(alpha,gamma)16O reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ', ' width ',
     &  reswidth,' MeV ', 'level ', elevel,' MeV'
C
C--     Set common block variables for cross-section
C
        m1 = abeam*amugev + devmass
        m2 = hemass
        z1 = zbeam
        z2 = 2.
        er = resenerg
        gp = 70./1000.
        gg = 0.007/1000.
        omg = 5.
        ell = 1.
        ires = 81
C
C--     Define states in O16 gamma decays from resonance
C
        devmass = -4736.998E-6
        prodm   = aprod*amugev + devmass
C
        elevel  = 11.260
        aamass  = prodm + elevel/1000.
        tlif    = hbar/2500.E-06
C
        CALL gspart(82,'12_O16_0+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 11.0967
        aamass = prodm + elevel/1000.
        tlif   = hbar/0.28E-06
C
        CALL gspart(83,'11_O16_4+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 11.080
        aamass = prodm + elevel/1000.
        tlif   = hbar/12.E-06
C
        CALL gspart(84,'10_O16_3+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 10.957
        aamass = prodm + elevel/1000.
        tlif = 5.5E-15
C
        CALL gspart(85,'9_O16_0-',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 10.356
        aamass = prodm + elevel/1000.
        tlif = hbar/26.E-06
C
        CALL gspart(86,'8_O16_4+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 9.8445
        aamass = prodm + elevel/1000.
        tlif = hbar/0.62E-06
C
        CALL gspart(87,'7_O16_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 9.585
        aamass = prodm + elevel/1000.
        tlif = hbar/420.E-06
C
        CALL gspart(88,'6_O16_1-',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 8.8719
        aamass = prodm + elevel/1000.
        tlif = 125.E-15
C
        CALL gspart(89,'5_O16_2-',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 7.11685
        aamass = prodm + elevel/1000.
        tlif = 8.3E-15
C
        CALL gspart(90,'4_O16_1-',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 6.9171
        aamass = prodm + elevel/1000.
        tlif = 4.7E-15
C
        CALL gspart(91,'3_O16_2+',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 6.12989
        aamass = prodm + elevel/1000.
        tlif = 18.4E-12
C
        CALL gspart(92,'2_O16_3-',8,aamass,zprod,tlif,ubuf,nubuf)
C
        elevel = 6.0494
        aamass = prodm + elevel/1000.
        tlif = 67.E-12
C
        CALL gspart(93,'1_O16_0+',8,aamass,zprod,tlif,ubuf,nubuf)
C
C--     ground state --> idpart = 104
        irecoil = 94
C
        tlif = 1000.
C
        CALL gspart(94,'gs_O16_0+',8,prodm,zprod,tlif,ubuf,nubuf)
        print*, 'Ground state mass: ', prodm
C
C--     branch info -- resonance decays
C

        brat(1) = 90.99
        mode(1) = 1 + 100*94
        brat(2) = 4.19
        mode(2) = 1 + 100*93
        brat(3) = 4.00
        mode(3) = 1 + 100*91
        brat(4) = 0.82
        mode(4) = 1 + 100*90
C
        CALL uzero(brat,5,6)
        CALL uzero(mode,5,6)
        CALL gsdk(81,brat,mode)
C
C--     12th ex. state
C
C
        CALL uzero(brat,1,6)
        CALL uzero(mode,1,6)
C
        CALL gsdk(82,brat,mode)
C
C--     11th ex. state
C
        brat(1) = 55.25
        mode(1) = 1 + 100*92
        brat(2) = 44.75
        mode(2) = 1 + 100*91
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
        CALL gsdk(83,brat,mode)
C
C--     10th ex. state
C
        CALL uzero(brat,1,6)
        CALL uzero(mode,1,6)
        CALL gsdk(84,brat,mode)
C
C--     9th ex. state
C
C
        brat(1) = 100.
        mode(1) = 1 + 100*90
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
        CALL gsdk(85,brat,mode)
C
C--     8th ex. state
C
        brat(1) = 1.57
        mode(1) = 1 + 100*92
        brat(2) = 98.43
        mode(2) = 1 + 100*91
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
        CALL gsdk(86,brat,mode)
C
C--     7th ex. state
C
        brat(1) = 60.98
        mode(1) = 1 + 100*94
        brat(2) = 18.29
        mode(2) = 1 + 100*93
        brat(3) = 20.73
        mode(3) = 1 + 100*91
C
        CALL uzero(brat,4,6)
        CALL uzero(mode,4,6)
        CALL gsdk(87,brat,mode)
C
C--     6th ex. state
C
        brat(1) = 89.29
        mode(1) = 1 + 100*94
        brat(2) = 10.71
        mode(2) = 1 + 100*91
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
        CALL gsdk(88,brat,mode)
C
C--     5th ex. state
C
        brat(1) = 7.22
        mode(1) = 1 + 100*94
        brat(2) = 0.12
        mode(2) = 1 + 100*93
        brat(3) = 77.67
        mode(3) = 1 + 100*92
        brat(4) = 3.57
        mode(4) = 1 + 100*91
        brat(5) = 11.42
        mode(5) = 1 + 100*90
C
        CALL uzero(brat,6,6)
        CALL uzero(mode,6,6)
        CALL gsdk(89,brat,mode)
C
C--     4th ex. state
C
        brat(1) = 99.93
        mode(1) = 1 + 100*94
        brat(2) = 0.07
        mode(2) = 1 + 100*92
C
        CALL uzero(brat,3,6)
        CALL uzero(mode,3,6)
        CALL gsdk(90,brat,mode)
C
C--     3rd ex. state
C
        brat(1) = 99.97
        mode(1) = 1 + 100*94
        brat(2) = 0.027
        mode(2) = 1 + 100*93
        brat(3) = 0.008
        mode(3) = 1 + 100*92
C
        CALL uzero(brat,4,6)
        CALL uzero(mode,4,6)
        CALL gsdk(91,brat,mode)
C
C--     2nd ex. state
C
        brat(1) = 100.
        mode(1) = 1 + 100*94
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
        CALL gsdk(92,brat,mode)
C
C--     1st ex. state
C  
        brat(1) = 100.
        mode(1) = 1 + 100*94
C
        CALL uzero(brat,2,6)
        CALL uzero(mode,2,6)
        CALL gsdk(93,brat,mode)



C
C--> Case 20: reaction drawn from user input card
C    only implemented for (a,g) or (p,g) reactions at present
C    CR: 22.07.2003
C    Now implemented for (c12,g) reactions.
C    JS: 15.02.2004
        case(20)
C
        INPUT = 'INPUT'
        CALL getenv(INPUT,CARDNAME)
        IF(CARDNAME.eq.'   ')CARDNAME='c12ag.dat'
        INQUIRE(file=CARDNAME,exist=fexist)
        IF(fexist)then
          open(20,file=CARDNAME,status='old')
        Else
          print*, 'No input reaction card! (LKINE set to 20)'
          stop
        Endif
        read(20,NML=params)
        close(20)

C    Setup beam particle
C
        aprod = atarg + abeam
        devmass = beam_mass_excess
        aamass  = abeam*amugev + devmass
        tlif    = beamlifetime
        ubuf(1) = fkine(1)
C
        CALL gspart(80,beamtyp//'gs',8,aamass,zbeam,tlif,ubuf,nubuf)
C
C    Create resonant particle
C
C    Note: Resonant particle energy based on beam mass, target mass and
C          variable resenerg.  Not based on energy levels specified in INPUT
C          file.  (However, energy of states below resonance are based on
C          eneryg levels specified in INPUT)
C
        If(atarg.eq.4.)then
        targtyp = 'a'
        targmass = hemass
        aamass   = aamass + resenerg/1000. + hemass
        print*, aamass, hemass, resenerg
        Elseif(atarg.eq.1)then
        targtyp = 'p'
        targmass = hmass
        aamass   = aamass + resenerg/1000. + hmass
        Elseif(atarg.eq.12)then
        targtyp = '12C'
        targmass = c12mass
        aamass   = aamass + resenerg/1000. + c12mass 
        Endif
        tlif     = hbar/((part_width+gam_width)/1000.)
        ubuf(1)  = fkine(2)
C
        CALL gspart(81,'res_'//rectyp//'_'//num(rstate)//'',8,
     &              aamass,zprod,tlif,ubuf,nubuf)
        ires = 81
        elevel = level(rstate)
        print*, 'Resonant mass: ', aamass
C        rmass = aamass
        resmass = aamass

        

        write(6,*)'|**** '//beamtyp//'('//targtyp//',gamma)'//rectyp//
     &           ' reaction ****|'
        write(6,*)'Resonance energy', resenerg, ' MeV ',
     &   'level ', elevel,' MeV'  

C
C    Setup excited states below resonance
C
        prodm   = aprod*amugev + recoil_mass_excess
        DO i = 0, rstate-1

        elevel = level(i)
        aamass = prodm + elevel/1000.
        tlif = life(i)
C
        CALL gspart(81+rstate-i,''//rectyp//'_'//num(i)//'',8,
     &              aamass,zprod,tlif,ubuf,nubuf)
        ENDDO
        irecoil = 81+rstate
C
C    Setup decay branching ratios and modes
C
        DO i = 1, rstate
        CALL uzero(brat,1,10)
        CALL uzero(mode,1,10)
         DO j = 1, 10
         if(br(i,j).ne.0)then
          brat(j) = br(i,j)
          mode(j) = 1 + 100*(irecoil-md(i,j))
         endif
         ENDDO
        CALL gsdk(81+rstate-i,brat,mode)
        ENDDO

C
C    Set cross-section variables
C
        gp = part_width
        gg = gam_width
        er = resenerg
        omg = spin_stat_fac
        ell = ell
        m1 = abeam*amugev + beam_mass_excess
        print*, m1
        If(atarg.eq.1.)then
        m2 = hmass
        print*, '++++++', m2
        Elseif(atarg.eq.4.)then
        m2 = hemass
        Elseif(atarg.eq.12.)then
        m2 = c12mass
        Endif
        z1 = zbeam
        z2 = ztarg
        
C     alpha source
        case(21)
        alpha = .true.
        atarg = 1
        aamass = hemass
        ubuf(1) = FKINE(1)
        ubuf(2) = FKINE(2)
        tlif = 1000.
        CALL gspart(80,'Alpha',8,aamass,2.0,tlif,ubuf,nubuf)

C     43K(p,n)43Ca
      case(22)
         zbeam = 19.
         abeam = 43.
         atarg =  1.
         zprod = 20.
         aprod = atarg + abeam
         targmass = hmass
         ilight = 13            !neutron emission
C     
C--   create beam particle --> idpart = 80
C      
         devmass = -3.65753926E-02 ! 43K
         aamass  = abeam*amugev + devmass
         tlif    = 1000.
         ubuf(1) = fkine(1)
C     
         CALL gspart(80,'K43',8,aamass,zbeam,tlif,ubuf,nubuf)
C     
C--   NON resonant level populated --> idpart = 81
         tlif     = hbar/10.e-3 ! 10 keV width (non-resonant)
         resenerg = 0 
         reswidth = hbar/tlif
         print*,'aamass(0),resenerg,reswidth',aamass, resenerg, reswidth
         aamass = sqrt((aamass+targmass)**2 + 
     &        2.*targmass*beamenerg*.001)
         print*,"aamass(1),targmass: ",aamass,targmass
         ubuf(1)  = fkine(2)

         
C     
         CALL gspart(81,'nonres_44Ca',8,aamass,zprod,tlif,ubuf,nubuf)
C     
C--   Define states in 43Ca from neutron decays
C     
         devmass = -3.84088273E-02 ! 43Ca
         prodm   = (aprod-1)*amugev + devmass !gs of 43Ca
         write(6,*) "43Ca MASS, A: ", prodm, aprod-1
C     
C     neutron decay to state in 43Ca
C     particle 13 is the neutron
C     Decide what to do next based on elevel, set from fkine(3)
C         
C     First create particles for the g.s. and first 18 excited
C     states (max # that are accessible @1.5 MeV/u). The g.s.
C     is always the final state, so set that to irecoil
C
C     Ground state (irecoil)
         irecoil = 82
         elevel = 0
         tlif = 1000.
         CALL gspart(irecoil,'43Ca_gs',8,prodm+elevel/1E3,
     &        zprod,tlif,ubuf,nubuf)
C
C     Excited states
         tlif = 1.E-15
         DO i=1, NLEV-1
            elevel = Ex43Ca(i)
            CALL gspart(irecoil+i,'43Ca_ex'//num(i),8,prodm+elevel/1E3,
     &           zprod,tlif,ubuf,nubuf)
            print*, 'EX::',i,elevel,irecoil+i
         ENDDO
         elevel = 0.
         ilevel = NINT(fkine(3))
C     Check bounds: valid entries are 0–18 inclusive
         IF (ilevel .GE. 0 .AND. ilevel .LE. 18) THEN
            elevel = Ex43Ca(ilevel)
            PRINT *, 'Read ELEVEL: ', elevel
        ELSE
            PRINT *, 'Error: fkine(3) =', fkine(3),
     &           ' out of valid range (0–18)'
            STOP
         ENDIF

C
C     First Decay 44Ca --> Chosen Excited State in 43Ca
C
         CALL uzero(brat,1,10)
         CALL uzero(mode,1,10)         

         brat(1)= 100.
         mode(1)= 13 + 100*(irecoil+ilevel)
         CALL gsdk(81,brat,mode)
         write(6,*)'|**** 43K(p,n)43Ca reaction ****|'
         write(6,*)' n decay to level',ilevel,
     +        ' particle:',irecoil+ilevel,' ex:',Ex43Ca(ilevel)
        
         
C     
C     Now set up gamma decays
C     
         DO i = 1, ilevel
            CALL uzero(brat,1,10)
            CALL uzero(mode,1,10)
            DO j = 1, MAXBR
               brat(j) = Br43Ca(i,j)
               mode(j) = 1 + 100*(irecoil+St43Ca(i,j))
            ENDDO
            CALL gsdk(irecoil+i,brat,mode)
         ENDDO

      End select
C
      RETURN
      END
C










