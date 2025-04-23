#Mach Variable defines,can defined in environment variables
# surport: default,JNCSSW2,SW4CC,SW5CC,THGZ,THGZMIC
#for Make mode ;surport nml,nmlc,omp,sw,mic

MACH   ?= NORM
#MACH    ?= HWARM
EXE     ?= nwam
NO_MPI  ?= 1
DEBUG		?= 0
minimake?= 0
defmake ?= irrp
NJ      ?= -j 10
hpctooldir ?= /pacific_fs/HPCKit/latest
hpcthirddir?= /pacific_fs/HPCKit/third
LIBNETCDF  ?=-lnetcdf -lnetcdff
ifeq ($(DEBUG),1)
	CDEBUG=-g -DDEBUG
	FDEBUG=-g -DDEBUG #-fcheck=bounds
	LDEBUG=-g
else
	CDEBUG= -O2 -g -DNODBGINFO
	FDEBUG= -O2 -g -DNODBGINFO #-fcheck=bounds
	LDEBUG= -g
endif
ifeq ($(LPID),imptest)
	NO_MPI=1
endif
NO_MPI := $(strip $(NO_MPI))
ifeq ($(NO_MPI),1)
	DFLAG := $(DFLAG) -DNO_MPI
	PR    := n$(PR)
endif
.PHONY:defmake all zall v1 v2 v3 v4

defmake:$(defmake)
zall:v1 v2 v3 v4

irrp:
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
v1:
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
v2:
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
v3:
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all
v4:
	OD=objs_$(PR)$@; mkdir -p $${OD};cd $${OD};make $(NJ) -f ../src/Makefile EXE=nwam_$(PR)$@ LPID=$@ all

clean_v1:
	rm -rf objs_$(PR)v1 nwam_v1
clean_v2:
	rm -rf objs_$(PR)v2 nwam_v2
clean_v3:
	rm -rf objs_$(PR)v3 nwam_v3
clean_v4:
	rm -rf objs_$(PR)v4 nwam_v4
clean_all:    clean_v1 clean_v2 clean_v3 clean_v4 
clean : clean_$(defmake)

#MODEL_DEF:=  -DC_CALCULATE -DNPROPAGAT #  -DOCHECK
MODEL_DEF :=

#default :
FC = mpif90
CC = mpicc
CXX =g++
CPP = cpp
LD = mpif90
DFLAG := $(DFLAG) -DFORTRANUNDERSCORE
FFLG = -cpp -D_LFORTRAN  -ffree-line-length-0
CFLG =
FMATH=
#OMPFLG=-openmp

MACH := $(strip $(MACH))
LPID := $(strip $(LPID))
MACH_State=default

ifeq ($(MACH),NORM)
 MACH_State=$(MACH)
	NETCDF_BASE ?=/pacific_fs/hpctool/netcdf/netcdf
endif
ifeq ($(MACH),HWARM)
	NETCDF_BASE ?=/pacific_fs/hpctool/netcdf/netcdf
 MACH_State=$(MACH)
 FC=mpif90 -mcpu=hip11
 #export NAN_REPORT_CONFIG=error+detail
 #-mcpu=linxicore6100
 CFLG :=$(CFLG) -mcpu=hip11 -march=armv9.2-a+sve2
 CFLG :=$(CFLG) -I$(hpctooldir)/kml/bisheng/include -I$(hpcthirddir)/include
 #-finstrument-functions
 CC =mpicc
 CXX =$(CC)
 CPP=cpp
 LD = mpif90  -mcpu=hip11 -lfpe --rtlib=compiler-rt
 FFLG := $(FFLG) -g #-fp-detect-exceptions
 SLIBS= -L$(hpcthirddir)/lib64 -lpthread -lksvml -lmemkind -lsdma_dk -lhwloc
 #SLIBS= -L$(hpctooldir)/kml/bisheng/lib/sve512 -lksvml
endif


ifdef MACH
	ifeq ($(MACH_State),default)
	  $(error MACH SET error ,$(MACH) is not surport ,you can change it or add it)
	endif
endif

LPID := $(strip $(LPID))

OBJD=objs_$(LPID)

endif
ifeq ($(LPID),v2)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE
	CINC:= $(CINC) zarm.h
	CSRC   = svars.c
	CPPSRC =  propagat.cpp implsch.cpp output_cal.cpp boundary.cpp
endif
ifeq ($(LPID),v3)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DMTHREAD
	CINC:= $(CINC) zarm.h
	CINC:= $(CINC) mthread.h

	CSRC   = svars.c mthread.c
	CPPSRC =  propagat.cpp implsch.cpp output_cal.cpp boundary.cpp
endif
ifeq ($(LPID),v4)
	MODEL_DEF :=  $(MODEL_DEF) -DC_CALCULATE -DMMTHREAD
	CINC:= $(CINC) zarm.h
	CINC:= $(CINC) mthread.h
	CINC:= $(CINC) gpart.h
	#FSRC = split_mod.F90
	CSRC   = svars.c mthread.c
	CPPSRC =  propagat.cpp implsch.cpp output_cal.cpp boundary.cpp gpart.cpp
endif

DFLAG := $(DFLAG) $(MODEL_DEF)
INCLDIR := -I. -I../src $(INCLDIR)

FFLAGS = $(DFLAG) $(FDEBUG) $(FFLG) $(INCLDIR) $(UOMPFLG)
CFLAGS = $(DFLAG) $(CDEBUG) $(CFLG) $(INCLDIR) $(UOMPFLG)
LFOPTS = $(SLIBS) $(LDEBUG) -L${NETCDF_BASE}/lib $(LIBNETCDF) $(UOMPFLG)

#FILE Structure:
# main.F90
#  |_ waveinit_mod.F90        InitWaveMdl
#  |  |_wav_mpi_mod.F90         1:InitMpi
#  |  |_varcommon_mod.F90       2:InitPara setwave ,use mpi
#  |  |_output_mod.F90          3:ReadTopog 5:InitMonitor 7:settopog 8:InitOutPut
#  |  |_partition_mod.F90       4:partition
#  |  |  |_use irrp_*.F90
#  |  |_platform_init_mod.F90   6:InitAfterPartition 12:Init_Distinct
#  |  |_implsch_mod.F90         9:InitImplsch
#  |  |_propagat_mod.F90        10:InitPropagat
#  |  |_boundary_mod.F90        11:InitSetspec
#  |  |_boundary_mod.F90        12:setspec
#  |  |_output_mod.F90          13:InitCheckOutPut
#  |  |_wav_mpi_mod.F90         14:SetModelTime
#  |  |_varcommon_mod.F90       15:NextTime
#  |
#  |_ wavemdl_mod.F90         RunWaveMdl
#  |  |_
#  |  |_
#  |  |_
#  |  |_
#  |
#  |_ waveinit_mod.F90        EndWaveMdl
#     |_End_Distinct
#     |_InitOutPut(-1)
#     |_IDataIO
#     |_InitMpi


ifneq ($(NO_MPI),1)
	SRCT = wav_mpi_mod.F90 irrp_smpi_mod.F90 irrp_split_mod.F90 irrp_kernal_mod.F90  irrp_package_mod.F90
	CSRCT = mympi.c
endif
	CSRCT := zfile.c $(CSRCT) ctools.c  std.c
	SRCT :=$(SRCT) wav_cal_mod.F90 netcdf_mod.F90
	SRC =  debughlp_mod.F90  varcommon_mod.F90 \
	      partition_mod.F90 dataio_mod.F90 windin_mod.F90 \
	      output_cal_mod.F90 output_loc_mod.F90 output_nc_mod.F90 output_mod.F90 \
	      boundary_mod.F90 implsch_mod.F90  propagat_mod.F90 \
	      restart_mod.F90  $(FSRC) \
	      platform_init_mod.F90 waveinit_mod.F90 wavemdl_mod.F90 main.F90
# msg_mod.F90

	FOBJS = ${SRC:.F90=.o}
	COBJS = $(CPPSRC:.cpp=.o) $(CSRC:.c=.o)
	FOBJST = ${SRCT:.F90=.o}
	COBJST = $(CSRCT:.c=.o)

	OBJST =$(FOBJST) $(COBJST) #$(SSOBJS) $(SCOBJS)
	OBJSM =$(FOBJS) $(COBJS) #$(SSOBJS) $(SCOBJS)
	OBJS =$(OBJST) $(OBJSM) #$(SSOBJS) $(SCOBJS)

#============== Set src file Path =====================
# make
ifdef EDIR
VPATH = ../src/$(LPID)/$(EDIR) ../src/$(LPID) ../src/io ../src/tools ../src
else
VPATH = ../src/$(LPID) ../src/io ../src/tools ../src
endif

#default/dump target:
#check  directory is exist,and make real/dump target in objs directory,
#and copy EXE file to current directory

.SUFFIXES:
.SUFFIXES:.F90 .f90 .o .c .oo .ss .s

%.o:%.F90
	$(FC) $(FFLAGS) -c $<

%.o:%.c
	$(CC) $(CFLAGS) -c $<

%.o:%.cpp
	$(CXX) $(CFLAGS) -c $<

$(OBJSM) : ../src/wavedef.h


CPP1 = $(CPP) -P -w  -traditional-cpp # -U__FILE__ -U__LINE__
CPPE = $(CPP) -P -w  -traditional-cpp
# Function: -P : Inhibit generation of linemarkers.
#           -w : Avoid warning due to undefine __FILE__ and __LINE__.
#           -U__FILE__ -U__LINE__ : undefine __FILE__ and __LINE__ ,for use in future.

SEDP = sed -e 's/\\$$/~\\/g'
# Function: Add char '~' before  '\' at end of line
#    where: \\ used as \ at command line.
#           $$ treated as $ in makefile.
#           $ means end of line at command of sed.

SEDN = sed -e " s/ *;~/\n/g ; s/~/\&\n/g ; s/^ *&//g;s/\#\#\s*//g"
SEDE = sed -e "/^\s*$$/{ N;  /^\s*$$/D }"
# Function: (1) Replace ";~" and any spaces into a newline
#           (2) Replace "~" into "&\n"
#           (3) Delete spaces before & for the case that there
#               is none valid characters before this &.

irrp_package_mod.o:irrp_package_mod.F90
	$(SEDP) $< | $(CPP1) -D__FILE__=\"$*.F90\" | $(SEDN) | $(SEDE)  >$*.f90
	$(FC) $(FFLAGS) -c $*.f90

netcdf_mod.o:netcdf_mod.F90
	$(SEDP) $< | $(CPP1) -D__FILE__=\"$*.F90\" | $(SEDN) | $(SEDE) |$(CPPE) -D__FILE__=\"$*.f90\" >$*.f90
	$(FC) $(FFLAGS) -I$(NETCDF_BASE)/include -c $*.f90

ifneq ($(minimake),1)
ifneq ($(NO_MPI),1)
irrp_split_mod.o :irrp_smpi_mod.o
irrp_kernal_mod.o:irrp_smpi_mod.o irrp_split_mod.o
irrp_package_mod.o:irrp_kernal_mod.o irrp_smpi_mod.o
partition_mod.o  :wav_mpi_mod.o irrp_split_mod.o irrp_smpi_mod.o irrp_kernal_mod.o  irrp_package_mod.o
debughlp_mod.o   : wav_mpi_mod.o
wav_mpi_mod.o    :
endif
gpart.o 				 :$(CINC)
svars.o 				 :$(CINC)
mthread.o				 :$(CINC)
propagat.o 			 :$(CINC)
implsch.o 			 :$(CINC)
output_cal.o 		 :$(CINC)
boundary.o			 :$(CINC)
wav_cal_mod.o    :
debughlp_mod.o   :
varcommon_mod.o  : defcheck debughlp_mod.o  netcdf_mod.o wav_cal_mod.o
boundary_mod.o   : varcommon_mod.o
implsch_mod.o    : varcommon_mod.o partition_mod.o
propagat_mod.o   : varcommon_mod.o partition_mod.o

output_cal_mod.o : varcommon_mod.o
output_nc_mod.o  : varcommon_mod.o partition_mod.o wav_cal_mod.o netcdf_mod.o
output_loc_mod.o : varcommon_mod.o partition_mod.o wav_cal_mod.o netcdf_mod.o
output_mod.o     : varcommon_mod.o output_cal_mod.o partition_mod.o output_nc_mod.o output_loc_mod.o   netcdf_mod.o

partition_mod.o  : varcommon_mod.o
dataio_mod.o     : varcommon_mod.o partition_mod.o netcdf_mod.o wav_cal_mod.o
windin_mod.o     : varcommon_mod.o dataio_mod.o wav_cal_mod.o

restart_mod.o   : varcommon_mod.o partition_mod.o windin_mod.o
platform_init_mod.o:varcommon_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o

waveinit_mod.o  : varcommon_mod.o platform_init_mod.o boundary_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o
wavemdl_mod.o   : varcommon_mod.o waveinit_mod.o platform_init_mod.o boundary_mod.o implsch_mod.o propagat_mod.o windin_mod.o output_mod.o partition_mod.o  restart_mod.o
main.o          : debughlp_mod.o wavemdl_mod.o

else
propagat_mod.o   : varcommon_mod.o
endif

$(OBJSM) :defcheck

#real/dump targets
all : $(EXE)

#check if defines are changed?
defcheck: ../src/Makefile
	@echo LPID= $(LPID) MACH=$(MACH)
	@echo $(FFLAGS) $(CFLAGS) >defcheck.new;cmp -s defcheck defcheck.new || mv defcheck.new defcheck

#real target,depend defcheck & files $(OBJS)
$(EXE) : defcheck $(OBJS)  ../src/Makefile
	$(LD) -o  $@ $(OBJS)  $(LFOPTS)
	cp $(EXE) ..
bak:
	sdate=`date +"%y%m%d_%H%M"`; \
	mkdir -p bak/src$${sdate} ;\
	cp -rp src/* bak/src$${sdate};)\
	(cd ../bak;zip -r src$${sdate} src$${sdate} )
