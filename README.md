
# This program is the source code of parallel UFDECOM-i based on DO CONCURRENT and OpenACC(GPUA-UFDECOM-i)
#Here are the instructions for use
#Before compiling, You need to load the NVHPC so that you can use the NVFORTRAN compiler
# For example,
#export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/bin:$PATH
#export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:$LD_LIBRARY_PATH
#export CPATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/include:$CPATH



# ------------------ COMPILE ----------------------------
#When compiling GPUA-UFDECOM-i on NVFORTRAN 25.3, set FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=mem:separate;  
#When compiling GPUA-UFDECOM-i on NVFORTRAN 24.7, set FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=nomanaged;  
#When compiling MC-UFDECOM-i, set FC = nvfortran -Mpreprocess -O4 -stdpar=multicore;  
#When compiling GPU-UFDECOM-i, set FC = nvfortran -Mpreprocess -stdpar=gpu.
#Makefile：

#FC = nvfortran -Mpreprocess -O4 -stdpar=multicore
FC = nvfortran -Mpreprocess -stdpar=gpu
#FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=nomanaged
#FC = nvfortran -Mpreprocess -stdpar=gpu -acc=gpu -gpu=mem:separate
#EXEC = UFDE_MC
EXEC = UFDE_G

MODS = MOD_GLOBAL.o MOD_ADVT_TEST1.o MOD_DATETIME.o \
MOD_WEIR.o MOD_HEAT.o MOD_OUTPUT.o MOD_RESTART.o \
MOD_SED.o MOD_ADVT.o MOD_TG.o mod_solver.o \

SUBS = ADVQ.o ADVU_TVD_3RD.o ADVV_TVD_3RD.o \
ALLOC_VARS.o BAROPG5b.o BCDATA.o BCOND.o \
BCOND_TWO.o BROUGH.o BROUGH_USERDEFINED.o CAL_EL.o \
CAL_MATERIAL.o CAL_SAL.o CAL_SED.o CAL_TMP.o DENS.o \
DENSD.o DISTINGUISH.o eltest_sor.o FIELD_CHECK.o \
FindCtrd.o FIND_INTERFACE_AG.o FIRST.o Flux2uv.o \
FLUX_BALANCE.o GEN_N.o INFO_EXCH.o \
GETEL_17_EBC.o GETEL_16_EBC.o TIMESTEP_GPU.o \
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

OBJS = $(MODS) $(SUBS)

mods_f = $(MODS:.o=.f)

$(EXEC) : $(OBJS)
	$(FC) -o $(EXEC) $(OBJS)

$(OBJS) : $(mods_f)
$(MODS) : MOD_GLOBAL.f
# ------------------ run ----------------------------
#run.sh:
#!/bin/sh
###BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -gpu "num=1"
#BSUB -n 1
#BSUB -q gpu
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -J run_gpu
nvidia-smi > out

./UFDE_G
#bsub <run.sh; run.sh is to run the script, you can copy the above required options to it.

