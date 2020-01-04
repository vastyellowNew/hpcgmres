include ../../Makefile.rule

OBJS =

ifeq ($(TARGET), GENERIC)
OBJS +=
endif

# double precision
OBJS += d_time_varying_smooth_horizon.o 
OBJS += d_givens_rotation.o 

# single precision
OBJS += s_time_varying_smooth_horizon.o 
OBJS += s_givens_rotation.o 

obj: $(OBJS)

clean:
    rm -f *.o
    rm -f *.s

d_time_varying_smooth_horizon.o: d_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
s_time_varying_smooth_horizon.o: s_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
d_givens_rotation.o: d_givens_rotation.c x_givens_rotation.c
s_givens_rotation.o: s_givens_rotation.c x_givens_rotation.c