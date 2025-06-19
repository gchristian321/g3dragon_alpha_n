c---- 67---- gamma angular distribution

      REAL FUNCTION angdist(X)

      REAL pi
      parameter (pi = 3.14159265358979)

C  A uniform angular distribution for gammas
      angdist = 1
CC  A dipole angular distribution for gammas
C      angdist = (3./(8.*pi))*(1.-X**2)
CC  A quad. angular distribution for gammas
C      angdist = (15./(8.*pi))*(1.-X**2)*X**2

      END

      REAL FUNCTION angdist1(X)

      REAL pi
      parameter (pi = 3.14159265358979)

C  A uniform angular distribution for gammas
      angdist1 = 1
CC  A dipole angular distribution for gammas
C      angdist1 = (3./(8.*pi))*(1.-X**2)
CC  A quad. angular distribution for gammas
C      angdist1 = (15./(8.*pi))*(1.-X**2)*X**2

      END

      REAL FUNCTION angdist2(X)

      REAL pi
      parameter (pi = 3.14159265358979)

C  A uniform angular distribution for gammas
C      angdist2 = 1
CC  A dipole angular distribution for gammas
      angdist2 = (3./(8.*pi))*(1.-X**2)
CC  A quad. angular distribution for gammas
C      angdist2 = (15./(8.*pi))*(1.-X**2)*X**2

      END

      REAL FUNCTION angdist3(X)

      REAL pi
      parameter (pi = 3.14159265358979)

C  A uniform angular distribution for gammas
C      angdist3 = 1
CC  A dipole angular distribution for gammas
C      angdist3 = (3./(8.*pi))*(1.-X**2)
CC  A quad. angular distribution for gammas
      angdist3 = (15./(8.*pi))*(1.-X**2)*X**2

      END
