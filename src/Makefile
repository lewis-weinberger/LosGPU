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
wilkes2: CFLAGS := -fast -ta=tesla:cc60 -Mlarge_arrays -mcmodel=medium  \
	               -Minfo=accel -std=c++11
wilkes2: TARGET := wilkes2

.PHONY: cpu
peta4: CFLAGS := -fast -ta=multicore -Minfo=accel -Mlarge_arrays  \
	             -mcmodel=medium -std=c++11
peta4: TARGET := peta4

# Build executable #############################################################

test:    $(EXEC)
ioa:     $(EXEC)
wilkes2: $(EXEC)
peta4:   $(EXEC)

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
