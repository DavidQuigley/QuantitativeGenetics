BOOST_DIR=/home/dquigley/boost_1_54_0
CC=g++
CFLAGS=-I$(BOOST_DIR) -Isrc -O3 -Wall -c -pthread
ODIR=obj
RM := rm -rf

SRC += \
../src/Attributes.cpp \
../src/ClassMinerOptions.cpp \
../src/DataStructures.cpp \
../src/Dataset.cpp \
../src/Discretize.cpp \
../src/Parser.cpp \
../src/Rawdata.cpp \
../src/graph.cpp \
../src/spear.cpp 

DEPS += \
../src/Attributes.h \
../src/ClassMinerOptions.h \
../src/DataStructures.h \
../src/Dataset.h \
../src/Discretize.h \
../src/ParseOptions.h \
../src/Parser.h \
../src/Rawdata.h \
../src/graph.h \
../src/spear.h 

OBJ += \
obj/Attributes.o \
obj/ClassMinerOptions.o \
obj/DataStructures.o \
obj/Dataset.o \
obj/Discretize.o \
obj/ParseOptions.o \
obj/Parser.o \
obj/Rawdata.o \
obj/graph.o \
obj/spear.o 
 
$(ODIR)/%.o: $(SRC) $(DEPS)
	$(CC) -c -o"$@" "$<" $(CFLAGS)

libCarmen: $(OBJ)
	@echo 'Building target: $@'
	ar -r "libCarmen.a" $(OBJ)
	@echo 'Finished building target: $@'
	@echo ' '

all: libCarmen

clean:
	-$(RM) $(OBJ) libCarmen.a
	-@echo ' '
