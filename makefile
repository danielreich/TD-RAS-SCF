F90 = gfortran
FFLAGS = -O3 -fopenmp -ffree-line-length-0 -fbounds-check

PATH_BLAS=BLAS
PATH_LAPACK=LAPACK
PATH_OBJ = obj_dicp

LIBBLAS=~/Programas/software/lapack-3.1.1/blas_LINUX.a
LIBLAPACK=~/Programas/software/lapack-3.1.1/lapack_LINUX.a 

.SUFFIXES : .o .f90 .f
.f90.o: ; $(F90) $(FFLAGS) -c $< 
.f.o:  ; $(F90) $(FFLAGS) -c $<

tdrasscf2014= sub.o  module_global.o log_information.o readin.o  module_auxiliary.o  module_fedvr3d_basis_set.o  module_index_radial.o  module_operator_radial.o module_operator_fedvr3d.o module_twoe_basis_set.o module_wf_juan.o module_operator_spatial_space_omp.o module_density.o module_solver.o module_analysis.o prop_a_phi.o  derivate.o combin_p_q.o space_a.o space_p.o space_q_juan.o orb_eqs_p_zet.o laser.o absorb.o main.o   

LIB = LINEAR
obj = $(PATH_OBJ)/lalib.o $(PATH_OBJ)/lineq.o $(PATH_OBJ)/mmlib.o \
      $(PATH_OBJ)/mtlib.o $(PATH_OBJ)/op1lib.o \
      $(PATH_OBJ)/op2lib.o $(PATH_OBJ)/rmlib.o \
      $(PATH_OBJ)/schmidtortho.o $(PATH_OBJ)/sdlib.o \
      $(PATH_OBJ)/xvlib.o $(PATH_OBJ)/rs.o

tdrasscf : $(tdrasscf2014) $(obj)  
	$(F90) $(FFLAGS)  -o  $@ $(tdrasscf2014) $(obj) $(LIBLAPACK) $(LIBBLAS)

#-------------------------------------------------------------
#Librery of MCTDHB
#-------------------------------------------------------------

$(PATH_OBJ)/schmidtortho.o : $(LIB)/schmidtortho.f
	$(F90) $(FFLAGS) -c $(LIB)/schmidtortho.f
	mv schmidtortho.o  $(PATH_OBJ)/

 $(PATH_OBJ)/lalib.o : $(LIB)/lalib.f
	$(F90) $(FFLAGS) -c $(LIB)/lalib.f 
	mv lalib.o $(PATH_OBJ)/
 $(PATH_OBJ)/lineq.o :$(LIB)/lineq.f
	$(F90) $(FFLAGS) -c $(LIB)/lineq.f
	mv lineq.o $(PATH_OBJ)/
 $(PATH_OBJ)/mmlib.o : $(LIB)/mmlib.f
	$(F90) $(FFLAGS) -c $(LIB)/mmlib.f
	mv mmlib.o $(PATH_OBJ)/
 $(PATH_OBJ)/mtlib.o : $(LIB)/mtlib.f
	$(F90) $(FFLAGS) -c $(LIB)/mtlib.f
	mv mtlib.o $(PATH_OBJ)/
 $(PATH_OBJ)/op1lib.o : $(LIB)/op1lib.f
	$(F90) $(FFLAGS) -c $(LIB)/op1lib.f
	mv op1lib.o $(PATH_OBJ)/
 $(PATH_OBJ)/op2lib.o : $(LIB)/op2lib.f
	$(F90) $(FFLAGS) -c $(LIB)/op2lib.f
	mv op2lib.o $(PATH_OBJ)/
 $(PATH_OBJ)/rmlib.o : $(LIB)/rmlib.f
	$(F90) $(FFLAGS) -c $(LIB)/rmlib.f
	mv rmlib.o $(PATH_OBJ)/
 $(PATH_OBJ)/sdlib.o :$(LIB)/sdlib.f
	$(F90) $(FFLAGS) -c $(LIB)/sdlib.f
	mv sdlib.o $(PATH_OBJ)/
 $(PATH_OBJ)/xvlib.o :$(LIB)/xvlib.f
	$(F90) $(FFLAGS) -c $(LIB)/xvlib.f
	mv xvlib.o $(PATH_OBJ)/
 $(PATH_OBJ)/rs.o :$(LIB)/rs.f
	$(F90) $(FFLAGS) -c $(LIB)/rs.f
	mv rs.o $(PATH_OBJ)/

#--------------------------------------------------

.PHONY:clean
clean:
	-rm tdrasscf
	-rm *.o
	-rm *.mod
#	-rm obj_dicp/*.o 