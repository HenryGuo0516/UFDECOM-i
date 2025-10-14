*=======================================================================
* include file "DEFS.h"
* This file contains a set of predetermined macro definitions which are 
* inserted into the individual files by C-preprocessor.
*
* By Sylvers Ding, Jan. 18, 2020
*=======================================================================  


* Turn on/off OpenMP compilation
* #define OMP


*=======================================================================
*                          MODULE SWITCHES
*=======================================================================
* Turn on/off salinity simulation
#define MODULE_SAL

* Turn on/off temperature simulation
#define MODULE_TMP

* Turn on/off sediment simulation
* #define MODULE_SED

* Turn on / off WAVE simulation
* #define MODULE_WAVE

* Turn on / off other material simulation
* #define MODULE_MATERIAL

* Turn on / off  lagrange tracking
* #define MODULE_LAG

* Turn on / off SLUICE
* #define MODULE_SLUICE


*=======================================================================
*                         FUNCTIONALITIES
*=======================================================================
* Turn on/off baroclinic pressure gradient term calculation
#define BPG

* Turn on/off the restart file archiving
#define RSTARC

* Choose bottom drag coefficient (BDC) type
* Available choices include:
* Z0BTYP_CONSTANT Z0BTYP_USERDEF
* If defining Z0BTYP_CONSTANT, BDC will be calculated by the logarithmic 
* formula
* If defining Z0BTYP_USERDEF, establish "_lba.dat" file
#define Z0BTYP_USERDEF

* Choose the horizontal viscosity and mixing parameterization scheme
* Available choices include:
* HORZMIX_CLOSURE
#define HORZMIX_CLOSURE

* Choose the vertical viscosity and mixing parameterization scheme
* Available choices include:
* VERTMIX_CLOSURE, VERTMIX_CONSTANT
#define VERTMIX_CLOSURE

* Choose advection scheme for transport process
* Available choices: 
* HSIMT HSIMTMPL KIM MPDATA SUPERBEE VANLEER MINIMOD CENTER UPWIND
#define HSIMT

* Add/remove weir
* #define WEIR

* Turn on/off model warning
* #define MODEL_WARNING

*=======================================================================
*                  SETTINGS OF SEDIMENT SIMULATION
*=======================================================================
* Turn on/off bathymetry update
* #define BED

* Choose tau type of sediment simulation
* Available choices include:
* TAUTYP_CONSTANT, TAUTYP_USERDEF
* If choosing TAUTYP_CONSTANT, tau will be constant in the whole field
* and be calculated using the parameters given in the .nml file
* If choosing TAYTYP_USERDEF, establsh "_d50_tau.dat" file in the input
* directory
* #define TAUTYP_USERDEF


*=======================================================================
*                   OPEN BOUNDARY CONDITIONS
*=======================================================================
* Turn on/off residual water level at open boundary.
* If defined, establish "rel_obc" file in input directory
#define ELB

* Choose the tidal boundary type
* Available choices are :
* TIDE_EL: drive by elevation
* TIDE_FLUX: drive by flux
* TIDE_FLATHER : drive by elevation and flux
# define TIDE_EL

* if consider shelf circulation (on this condition must choose Tide_flux)
* # define shelf_circulation

* Choose lateral boundary conditions for elevation
* If this option is avtivated, elevation at the open boundary will be 
* overwritten despite of the elevation given by tidal constants (EBC)
* and residual water level (ELB).
* Available choices are:
* LBCel_gra: gradient boundary condition
* LBCel_cha: Chapman boundary condition
* LBCel_rad: 2D radiation boundary condition
* #define LBCel_rad

* Choose lateral boundary conditions for velocity
* Available choices include:
* LBCuv_gra: gradient boundary condition
* LBCuv_cha: Chapman boundary condition (implicit)
* LBCuv_rad: 2D radiation boundary condition
#define LBCuv_cha

*=======================================================================
*                            OUTPUT
*=======================================================================
* Turn on/off the field output
#define FPT
#ifdef FPT
* # define FSMPT //turn on/off fsm field output
# define VPT //turn on/off velocity field output
# define EPT //turn on/off elevation field output
# define SPT //turn on/off salinity field output
# define MPT //turn on/off KM and KHfield output
# define TPT //turn on/off temperature field output
* # define SEDPT //turn on/off sediment field output
*#define BEDPT //turn on/off sediment bed field output
*#define MMTPT //turn on/off momentum field output
#endif


*=======================================================================
*                  METEOROLOGICAL CONFIGURATIONS
*=======================================================================
* Choose wind data type
* Available choices include:
* WDTYP_UNIFORM, WDTYP_FIELD
* If defining WDTYP_UNIFORM, establish an ASCII "_wds.dat" file
* If defining WDTYP_FIELD, establsih a binary "_wds.dat" file
#define WDTYP_FIELD

* Turn on / off switch air pressure
#define AIRPRESSURE


* Turn on/off switch of heat flux
* If defining HEATFLUX_SPECI, gives the heat flux directedly(CONNOT USE NOW)
* If defining HEATFLUX_BULK, calculate the heat flux ... 
* via sea / air - temperature, air pressure, relative humidity and wind speed.
#define HEATFLUX_BULK


*=======================================================================
*                   INFO-EXCHANGING OPERATORS
*=======================================================================
* Choose update/interpolation scheme for the infomation exchange process
* For update process, available choices include:
* UD_DR (directly-replacing), 
* UD_IDW (inverse distance weighting interpolation), 
* UD_AVE (area-averaging), 
* UD_SF (9-point Shapiro filtering),
* UD_FWO (full-weighting operator)
#define UD_FWO

* For interpolation precess, available choices include:
* INT_UNI (0th order uniform), 
* INT_IDW (inverse distance weighting interpolation), 
* INT_LINEAR (inverse bilinear interpolation), (not recommended)
* INT_QDT (quadratic interpolation),
* INT_HPI (HSIMT parabolic interpolation),
* INT_UAEI (upwind AEI), 
* INT_HAEI (HSIMT AEI) (not recommended)
* RECOMMENDATIONS:
* If choosing INT_QDT or INT_UAEI,it is recommended to set the refining
* multiple to odd numbers
#define INT_HPI

* Turn on/off Flux Correction Algorithm
*#define NIFluxCon


* Turn on / off MOD_ONLINE_VIEW
* #define MOD_ONLINE_VIEW

* USE IMPLICIT OR EXPLICIT SOLUTION FOR EL
* CHOICE EXP_EL OR IMP_EL
#define IMP_EL
