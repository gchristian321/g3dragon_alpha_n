# Special fork of DRAGON GEANT3 Simulation
## Implements 22Ne(a,n)25Mg reaction

This code originally forked from the `master` branch of:
https://github.com/DRAGON-Collaboration/G3_DRAGON

Some notes about the special implementation:
 - 22Ne(a,n)26Mg reaction is hard-coded as case(13)
 - Resonance energy is hard-coded for E_x = 11.83 MeV resonance
 - Resonance widths also hard coded
 - Central tune energy calculated from 22Ne(a,n) kinematics (needs checking still)
 - Outgoing recoil vectors should be tracked as normal for rad. capt, and saved in the same E_rec, cost_r, cosp_r hbook outputs
 - Outgoing neutron vectors stored as new hbook branch E_n, cost_n, cosp_n; neutrons not tracked in the simulation (AFAIK)
 - Pressure is set (hard-coded) to 4.87 torr to match 4.87(8) torr pressure in S2230 runs 12939-12946 [`ugmate_trgt.f`]

To run, first source `dsinit-ne22an.dat`, which looks like the following:
```
## Define necessary environment variables for BASH shell
export DSROOT="`pwd`"

export DSDAT="$DSROOT/dat"
export DSLIB="$DSROOT/lib"
export DSSOURCE="$DSROOT/src"
export DSINCLUDE="$DSROOT/include"
export DSBIN="$DSROOT/bin"
export MITRAY="$DSDAT/dragon_2014_IC.dat"
export FFCARD="$DSROOT/dragon_2003_ne22an.ffcards"
#export INPUT="$DSDAT/26alpg.dat"
```
Then ensure that `dragon_2003_ne22an.ffcards` makes sense, specifically
beam energy should be appropriate to reach the resonance in target; the
original version is:
```
BEAM 7.979
```
(note -- this has not been tested with BEAM unset, only with it set).

The ffcards file should also set the reacton case 13 and 4+ charge states
for both beam and recoil:
```
FKIN 13 4.0 4.0
```

## Magnetic & Electric field scaling
### 2025/06/18
A new option has been added to manually scale magnetic and electric fields.
These can be set with the SCAL variable in the ffcards file. For example:
```
SCAL  0.948416 0.8936756
```
The first argument is the scale factor for the B field, and the second for the E
field.

Note that the scaling factors are applied AFTER the scaling factors for the mistune
(MTUN), where the 5th argument is a per-cent deviation from the "standard" tune
(defined by recoils at 90-deg COM angle).  So the process is:
 1. Calculate energies of recoils at 90 DEG COM. Then calculate scaling factors
relative to reference tune (E=1.8885 MeV, A=19, Q=4) to produce fields that 
center these recoils.
 2. Multiply fields in (1) by the scaling factors set by SCAL in the ffcards file.
(If SCAL un-set, nothing happens in this step).

The code now also prints some information about the fields, including an estimate of
the MD1 values and ED1 setpoints. *These need to be double checked to ensure the
reference tune --> scaled tune fields are consistent with the equations used.**
To find these printouts easily do:
```
./bin/dsbatch | grep -A 10 -B 5 Reference
```
**Note that the re-scaling only works if BEAM is used in the ffcards file to set
the beam energy manually.** Hence, simulaitons of (a,n) equations should always use
this option.

## Angular distribution options
### 2025/06/19
This version adds some options to set angular distributions (center of mass)
using the ANGD variable in ffcards options are:
```
ANGD 0 <ignored> <ignored>
```
--> isotropic over full range of center of mass angles

```
ANGD 1 <min> <max>
```
--> isotropic with cos(THETA_CM) limited to <max> -> <min> range

```
ANGD 2 <min> <max>
```
--> dipole with cos(THETA_CM) limited to <max> -> <min> range

```
ANGD 3 <min> <max>
```
--> quadrupole with cos(THETA_CM) limited to <max> -> <min> range

If you want dipole or quadrupole over the full range of angles, then set <min> = -1, <max> = 1, e.g.:
```
ANGD 2 -1. 1.
```

BE SURE TO INCLUDE THE DECIMAL on the ranges to make them floating point!!


# **It is only recommended to use this version for specific 22Ne(a,n) simulations!!**

# Ordinary README for DRAGON GEANT3 proceeds below...


# Geant3

## Prerequisites:

The GEANT3 simulation of DRAGON requires [cernlib2006](http://cernlib.web.cern.ch/cernlib/index.html), [ROOT](https://root.cern.ch), and a g77 compiler.

## Quick Start:

To compile the batch and interactive versions of the simulation, cd to src and type
```
> make dsbatch && make dsinter
```

This will create the binary executables in the bin directory. To run an example simulation with 100 events using the 26Al(p,γ)27Si reaction
(i.e. - using the included 26alpg.dat file in the reaction definition card) type
```
> source dsinit.sh
> ./bin/dsbatch
```

This will produce the file dragon1.hbook. Convert it to a rootfile by typing
```
> h2root dragon1.hbook
```

## Detailed Instructions:

The detection efficiency of DRAGON's BGO γ-ray array as well as the DRAGON separator transmission can be calculated by running the GEANT3 
simulation of the DRAGON separator. To calculate these values for a given reaction at a given observation energy, one must set up an input 
file (e.g 26alpg.dat included in the repository). Up to 15 energy levels are allowed, each with up to 6 decay modes and branching ratios. 
Be sure to set the correct level energies (in MeV), lifetimes (in seconds), branching ratios (expressed in as a percentage of the sum of 
all the intensities in www.nndc.bnl.gov style units, so that all the branching ratios add up to 100%) and the decay modes. Only include 
the levels that actually have cascades from the resonance of interest, to minimize the number of levels. Be sure to set the variable 
'rstate' to the index of the energy level via which the resonance proceeds. Also be sure to set the C.M. resonance energy as Ex - Q. Check 
the file geant/src/angdist to see that the angular distribution is set to either 1 (isotropic) or a function (like the quadrupole one 
listed). Change to what you want and then (if necessary) recompile by typing
```
> make dsbatch && dsinter
```

Edit 'dsinit.sh' so that it lists the correct input file (i.e. - INPUT="$DSROOT/26alpg.dat"). Type 
```
> source dsinit.sh
```

Run about 5,000 events in GEANT (edit the 'TRIG" variable in 'dragon_2003.ffcards' to suit). One can run the simulation in the background 
using the command:
```
> nohup ./dsbatch &
```

When the run is finished, the output file will be called 'dragonXX.hbook', where XX is the number following RUNG defined in the
'dragon_2003.ffcards' file.  Convert this file to a root file of your choice by typing, for example:
```
> h2root dragon01.hbook 26alpg_368keV.root
```

Copy the output root file to the directory 'BGO Efficiency' in dragon@isdaq04:/home/dragon/
Submit the output root file to the DRAGON elog 'BGO Efficiency' in dragon@isdaq04:/home/dragon/
The file containing the thresholds, 'thresholds.root', is already in that directory, all you have to do is run the efficiency.C macro by 
typing:
```
root 26alpg_368keV.root
```

to start a ROOT session, then:
```
root> .L efficiency.C
```

to load the Macro, then:
```
root> efficiency("23napg_646keV.root")
```

You will be prompted to enter a number from 1 to 4 to choose which thresholds were used. Once this has been entered, the macro will spit 
out a number for the efficiency based on the total number of counts above threshold in all BGOs, plus a statistical error based on the 
number of detected counts. To run again for a different threshold, reload the macro and run again.


