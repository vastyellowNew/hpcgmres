include ../../Makefile.rule

OBJS =
OBJS +=

# double precision
OBJS += d_zero_horizon_ocp.o 
OBJS += d_inexact_newton_for_zero_horizon_ocp.o 
OBJS += d_mfgmres_for_cgmres_initializer.o 
OBJS += d_cgmres_initializer.o 

# single precision
OBJS += s_zero_horizon_ocp.o 
OBJS += s_inexact_newton_for_zero_horizon_ocp.o 
OBJS += s_mfgmres_for_cgmres_initializer.o 
OBJS += s_cgmres_initializer.o 

obj: $(OBJS)

clean:
    rm -f *.o
    rm -f *.s

d_zero_horizon_ocp.o: d_zero_horizon_ocp.c x_zero_horizon_ocp.c
s_zero_horizon_ocp.o: s_zero_horizon_ocp.c x_zero_horizon_ocp.c
d_inexact_newton_for_zero_horizon_ocp.o: d_inexact_newton_for_zero_horizon_ocp.c x_inexact_newton_for_zero_horizon_ocp.c
s_inexact_newton_for_zero_horizon_ocp.o: s_inexact_newton_for_zero_horizon_ocp.c x_inexact_newton_for_zero_horizon_ocp.c
d_mfgmres_for_cgmres_initializer.o: d_mfgmres_for_cgmres_initializer.c ../common/d_mfgmres.c ../common/x_mfgmres.c
s_mfgmres_for_cgmres_initializer.o: s_mfgmres_for_cgmres_initializer.c ../common/s_mfgmres.c ../common/x_mfgmres.c
d_cgmres_initializer.o: d_cgmres_initializer.c x_cgmres_initializer.c
s_cgmres_initializer.o: s_cgmres_initializer.c x_cgmres_initializer.c