SHELL=	/bin/sh

CMD= nsst.x

OBJS =  read_tf_sfcanl.o get_mask.o get_tim_wgt.o get_tf_clm.o get_seaice.o get_sal_clm.o lalo_to_cubed.o nc_check.o togrib.o tfreez.o sor_lapl.o iterate.o putsi.o setup.o nsst.o

build: $(CMD)

$(CMD):	$(OBJS)
	$(FCOMP) $(FFLAGS) -I$(NEMSIO_INC) $(NETCDF_INCLUDE) -o $(CMD) $(OBJS) $(NETCDF_LDFLAGS_F) $(NEMSIO_LIB) $(BACIO_LIB4) $(W3NCO_LIBd)

read_tf_sfcanl.o:  read_tf_sfcanl.f90
	$(FCOMP) $(FFLAGS) $(NETCDF_INCLUDE) -I$(NEMSIO_INC) -c read_tf_sfcanl.f90
get_mask.o:  get_mask.f90
	$(FCOMP) $(FFLAGS) -c get_mask.f90
get_tim_wgt.o:  get_tim_wgt.f90
	$(FCOMP) $(FFLAGS) -c get_tim_wgt.f90
get_tf_clm.o:  get_tf_clm.f90
	$(FCOMP) $(FFLAGS) -c get_tf_clm.f90
get_seaice.o:  get_seaice.f90
	$(FCOMP) $(FFLAGS) -c get_seaice.f90
get_sal_clm.o:  get_sal_clm.f90
	$(FCOMP) $(FFLAGS)  $(NETCDF_INCLUDE) -c get_sal_clm.f90
lalo_to_cubed.o:  lalo_to_cubed.f90
	$(FCOMP) $(FFLAGS) -c lalo_to_cubed.f90
nc_check.o:  nc_check.f90
	$(FCOMP) $(FFLAGS) $(NETCDF_INCLUDE) -c nc_check.f90
togrib.o:  togrib.f90
	$(FCOMP) $(FFLAGS) -c togrib.f90
tfreez.o:  tfreez.f90
	$(FCOMP) $(FFLAGS) -c tfreez.f90
sor_lapl.o:  sor_lapl.f90
	$(FCOMP) $(FFLAGS) -c sor_lapl.f90
iterate.o:  iterate.f90
	$(FCOMP) $(FFLAGS) -c iterate.f90
putsi.o:  putsi.f90
	$(FCOMP) $(FFLAGS) -c putsi.f90
setup.o:  setup.f90
	$(FCOMP) $(FFLAGS) -c setup.f90
nsst.o:  nsst.f90
	$(FCOMP) $(FFLAGS) -I$(NEMSIO_INC) -c nsst.f90

install:
	cp ${CMD} ../../exec/

clean:
	rm -f *.o *.mod ${CMD} 

test:
	@echo NO TESTS YET

