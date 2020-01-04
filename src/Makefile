include ../Makefile.rule

OBJS =

ifeq ($(TARGET), GENERIC)
OBJS +=
endif

# double
OBJS += d_zero_horizon_ocp.o 
OBJS += d_inexact_newton_for_zero_horizon_ocp.o 
OBJS += d_cgmres_initializer.o 
OBJS += d_time_varying_smooth_horizon.o 
OBJS += d_single_shooting_ocp.o 
OBJS += d_single_shooting_continuation.o 
OBJS += d_single_shooting_cgmres.o 

# single
OBJS += s_zero_horizon_ocp.o 
OBJS += s_inexact_newton_for_zero_horizon_ocp.o 
OBJS += s_cgmres_initializer.o 
OBJS += s_time_varying_smooth_horizon.o 
OBJS += s_single_shooting_ocp.o 
OBJS += s_single_shooting_continuation.o 
OBJS += s_single_shooting_cgmres.o 

obj: $(OBJS)

clean:
    rm -f *.o
    rm -f *.s

d_zero_horizon_ocp.o: d_zero_horizon_ocp.c x_zero_horizon_ocp.c 
s_zero_horizon_ocp.o: s_zero_horizon_ocp.c x_zero_horizon_ocp.c
d_inexact_newton_for_zero_horizon_ocp.o: d_inexact_newton_for_zero_horizon_ocp.c x_inexact_newton_for_zero_horizon_ocp.c
s_inexact_newton_for_zero_horizon_ocp.o: s_inexact_newton_for_zero_horizon_ocp.c x_inexact_newton_for_zero_horizon_ocp.c
d_cgmres_initializer.o: d_cgmres_initializer.c x_cgmres_initializer.c d_mfgmres.c x_mfgmres.c
s_cgmres_initializer.o: s_cgmres_initializer.c x_cgmres_initializer.c s_mfgmres.c x_mfgmres.c
d_time_varying_smooth_horizon.o: d_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
s_time_varying_smooth_horizon.o: s_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
d_single_shooting_ocp.o: d_single_shooting_ocp.c x_single_shooting_ocp.c
s_single_shooting_ocp.o: s_single_shooting_ocp.c x_single_shooting_ocp.c
d_single_shooting_continuation.o: d_single_shooting_continuation.c x_single_shooting_continuation.c
s_single_shooting_continuation.o: s_single_shooting_continuation.c x_single_shooting_continuation.c
d_single_shooting_cgmres.o: d_single_shooting_cgmres.c x_single_shooting_cgmres.c d_mfgmres.c x_mfgmres.c
s_single_shooting_cgmres.o: s_single_shooting_cgmres.c x_single_shooting_cgmres.c s_mfgmres.c x_mfgmres.c