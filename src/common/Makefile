include ../../Makefile.rule

OBJS =
OBJS +=

# double precision
OBJS += d_time_varying_smooth_horizon.o 
OBJS += d_givens_rotation.o 
OBJS += d_linear_algebra.o 
OBJS += d_memory_manager.o

# single precision
OBJS += s_time_varying_smooth_horizon.o 
OBJS += s_givens_rotation.o 
OBJS += s_linear_algebra.o 
OBJS += s_memory_manager.o

obj: $(OBJS)

clean:
    rm -f *.o
    rm -f *.s

d_time_varying_smooth_horizon.o: d_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
s_time_varying_smooth_horizon.o: s_time_varying_smooth_horizon.c x_time_varying_smooth_horizon.c
d_givens_rotation.o: d_givens_rotation.c x_givens_rotation.c
s_givens_rotation.o: s_givens_rotation.c x_givens_rotation.c
d_linear_algebra.o: d_linear_algebra.c 
s_linear_algebra.o: s_linear_algebra.c 
d_memory_manager.o: d_memory_manager.c x_memory_manager.c
s_memory_manager.o: s_memory_manager.c x_memory_manager.c