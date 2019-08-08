# LosGPU Makefile
# Possible targets: test, cpu, gpu
# Requires: PGI compiler, HDF5 and OpenMPI libraries
# 30/07/19

EXEC := ../LosGPU
OBJS := main.o useful.o output.o send_data.o read_snapshot.o 
INCL := global_vars.h prototypes.h

.DEFAULT_GOAL := test

# Set up environment variables #################################################

CC           := h5c++
HDF5_CXX     := mpic++
HDF5_CLINKER := mpic++

.PHONY: test
test: CFLAGS := -fast -ta=host -Mlarge_arrays -mcmodel=medium -std=c++11
test: TARGET := test

.PHONY: gpu
gpu: CFLAGS := -fast -ta=tesla:cc60 -Mlarge_arrays -mcmodel=medium  \
	               -Minfo=accel -std=c++11
gpu: TARGET := gpu

.PHONY: cpu
cpu: CFLAGS := -fast -ta=multicore -Minfo=accel -Mlarge_arrays  \
	             -mcmodel=medium -std=c++11
cpu: TARGET := cpu

# Build executable #############################################################

test:    $(EXEC)
gpu: $(EXEC)
cpu:   $(EXEC)

%.o: %.cpp $(INCL)
	$(CC) $(CFLAGS) -c $<

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	@echo Compiled for: $(TARGET)

# Clean ########################################################################
#
.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
