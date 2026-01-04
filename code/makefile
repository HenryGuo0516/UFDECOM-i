
# Compiler and compilation options
# The parallel strategy of UFDECOM-i is controlled solely
# by compiler flags, without modifying the source code.

# MC-UFDECOM-i Configuration:
FC = nvfortran -Mpreprocess -O4 -stdpar=multicore
# GPU-UFDECOM-i Configuration:
#FC = nvfortran -Mpreprocess -stdpar=gpu
# GPUA-UFDECOM-i Configuration(Platform 2):
#FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=nomanaged
# GPUA-UFDECOM-i Configuration(Platform 2):
#FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=mem:separate

# Executable for multi-core CPU version
EXEC = UFDE_MC
# Alternative executaable name for GPU versions
#EXEC = UFDE_G

# Object files: modules
# These contain global variables, data structures,
# and physical process modules.

MODS = MOD_GLOBAL.o MOD_ADVT_TEST1.o MOD_DATETIME.o \
MOD_WEIR.o MOD_HEAT.o MOD_OUTPUT.o MOD_RESTART.o \
MOD_SED.o MOD_ADVT.o MOD_TG.o mod_solver.o

# Object files: subroutines
# These include core numerical solvers, physical parameterizations, 
# boundary condition handling, and I/O routines.
SUBS = ADVQ.o ADVU_TVD_3RD.o ADVV_TVD_3RD.o \
ALLOC_VARS.o BAROPG5b.o BCDATA.o BCOND.o \
BCOND_TWO.o BROUGH.o BROUGH_USERDEFINED.o CAL_EL.o \
CAL_MATERIAL.o CAL_SAL.o CAL_SED.o CAL_TMP.o DENS.o \
DENSD.o DISTINGUISH.o eltest_sor.o FIELD_CHECK.o \
FindCtrd.o FIND_INTERFACE_AG.o FIRST.o Flux2uv.o \
FLUX_BALANCE.o GEN_N.o INFO_EXCH.o \
GETEL_17_EBC.o GETEL_16_EBC.o TIMESTEP_MC.o \
INFO_EXCH_GPU.o INFO_EXCH3D_GPU.o INFO_EXCH3D.o \
INIT_CBCCOF.o INIT_COF.o INSIDE.o INTERP_POINT.o \
MAXMIN.o OPEN_FILE.o PRE_IE.o PRE_TSR.o PROFQ.o \
PROFT.o PROFU.o PROFU_OLD.o PROFV.o PROFV_OLD.o \
PROF_SED.o READ_NML.o REUV.o SEC.o SETDOM.o \
SINTER.o SMAG.o SOLVER_SOR_TEST.o STEP2_SED2.o \
TANDS3C.o TIMESTEP.o TVD_LIMITERS.o UnFECOM.o \
UVCon.o VARBD.o VARBD_GPU.o VERTVL.o VGA_CHECK.o \
WEIGHT_IDW.o WREAL.o ZERO.o \
INTER_S_VERTICAL.o NMC_DIFF.o MOD_LAG.o \
CAL_TH.o NO_INVERSION.o

# Combined object list
OBJS = $(MODS) $(SUBS)

# Convert module object names to source file names
mods_f = $(MODS:.o=.f)

# Link all object files into the final executable
$(EXEC) : $(OBJS)
	$(FC) -o $(EXEC) $(OBJS)

# Compile rules for object files
$(OBJS) : $(mods_f)

# Ensure MOD_GLOBAL is compiled first
# (it defines global variables used by other modules)
$(MODS) : MOD_GLOBAL.f