include ../../Makefile.rule

OBJS =
OBJS +=

# double precision
OBJS += d_single_shooting_ocp.o 
OBJS += d_single_shooting_continuation.o 
OBJS += d_mfgmres_for_single_shooting_cgmres.o 
OBJS += d_single_shooting_cgmres.o 

# single precision
OBJS += s_single_shooting_ocp.o 
OBJS += s_single_shooting_continuation.o 
OBJS += s_mfgmres_for_single_shooting_cgmres.o 
OBJS += s_single_shooting_cgmres.o 

obj: $(OBJS)

clean:
    rm -f *.o
    rm -f *.s

d_single_shooting_ocp.o: d_single_shooting_ocp.c x_single_shooting_ocp.c
s_single_shooting_ocp.o: s_single_shooting_ocp.c x_single_shooting_ocp.c
d_single_shooting_continuation.o: d_single_shooting_continuation.c x_single_shooting_continuation.c
s_single_shooting_continuation.o: s_single_shooting_continuation.c x_single_shooting_continuation.c
d_mfgmres_for_single_shooting_cgmres.o: d_mfgmres_for_single_shooting_cgmres.c x_mfgmres_for_single_shooting_cgmres.c
s_mfgmres_for_single_shooting_cgmres.o: s_mfgmres_for_single_shooting_cgmres.c x_mfgmres_for_single_shooting_cgmres.c
d_single_shooting_cgmres.o: d_single_shooting_cgmres.c x_single_shooting_cgmres.c
s_single_shooting_cgmres.o: s_single_shooting_cgmres.c x_single_shooting_cgmres.c