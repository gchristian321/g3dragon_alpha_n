C  RESCOM
       REAL*4  eres,eres0, resenerg, reswidth, abeam, atarg, prodm,
     + entdens, exitdens,z_react,erescm, rmass, resmass, targmass
       Integer*4 irecoil, ires
       LOGICAL resflag, McpHit, alpha
       COMMON / rescom / eres, eres0, resenerg, reswidth,
     + abeam, atarg, prodm, entdens, exitdens,irecoil,resflag,z_react,
     + ires,erescm, rmass, resmass, targmass, alpha
C

      REAL ptorr, ttemp, tdens
C
      COMMON /PTORR/ ptorr, ttemp, tdens


      REAL dedx_scl_b, dedx_scl_r
      COMMON /DEDXSCL/ dedx_scl_b, dedx_scl_r
      