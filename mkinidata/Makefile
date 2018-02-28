# Makefile for CLM initial data processing

include ../include/Makeoptions

# name of executable
EXECUTABLE = initial.x

####################################################################

OBJS =  precision.o              \
	MOD_TimeInvariants.o     \
	timemanager.o            \
	PhysicalConstants.o      \
	MOD_TimeVariables.o      \
	CLMINI.o                 \
	initialize.o             \
	ALBEDO.o                 \
	IniTimeConst.o           \
	IniTimeVar.o             \
	lai_empirical.o          \
	LAI_readin.o             \
	lakedepth_readin.o       \
	orb_coszen.o             \
	snowfraction.o           \
	soil_color_refl.o        \
	soil_parameters_readin.o

####################################################################

all: ${EXECUTABLE}

${EXECUTABLE} : ${OBJS} 
	${FF} ${FOPTS} ${OBJS} -o $@
	@echo 'making CLM initial data completed!'

$(OBJS):%.o:%.F90
	${FF} -c ${FOPTS} $(INCLUDE_DIR) -o $@ $< 

clean :
	${RM} -f *.o *.mod
