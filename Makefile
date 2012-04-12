.SUFFIXES : .o .f90 .f .F90

F90	:= mpif90
F90FLAGS := -g -W -fbounds-check -Jmod
LINK_FLAGS :=
LIBS :=
CPP_DIRECTIVES :=
INC_DIR :=

OBJ_FILES = $(wildcard obj/*.o)

.PHONY: depend clean
xMHD3D: src/mhd3d.f90
	$(F90) ${F90FLAGS} $< -o $@ \
			$(OBJ_FILES) \
			$(LIBS) \
			${LINK_FLAGS} \
			${CPP_DIRECTIVES} \
			${INC_DIR}

obj/%.o: src/%.f90
	${F90} -c ${F90FLAGS} $< -o $@ ${CPP_DIRECTIVES} ${INC_DIR}

obj/%.o: src/%.F90
	${F90} -c ${F90FLAGS} $< -o $@ ${CPP_DIRECTIVES} ${INC_DIR}

include Makefile.deps

clean:
	rm -f xMHD3D mod/*.mod obj/*.o
