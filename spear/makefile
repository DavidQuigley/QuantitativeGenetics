ODIR=obj
BINDIR=bin
$(shell mkdir -p $(ODIR)) 
$(shell mkdir -p $(BINDIR)) 
SRC_DIR=src
BOOST_DIR=/home/david/software/boost_1_57_0
CC=g++
CFLAGS=-I$(BOOST_DIR) -I$(SRC_DIR) -O2 -Wall -pthread
RM := rm -rf

SRC += \
$(SRC_DIR)/Attributes.cpp \
$(SRC_DIR)/ClassMinerOptions.cpp \
$(SRC_DIR)/DataStructures.cpp \
$(SRC_DIR)/Dataset.cpp \
$(SRC_DIR)/Discretize.cpp \
$(SRC_DIR)/Parser.cpp \
$(SRC_DIR)/ParseOptions.cpp \
$(SRC_DIR)/Rawdata.cpp \
$(SRC_DIR)/graph.cpp \
$(SRC_DIR)/spear.cpp

DEPS += \
$(SRC_DIR)/Attributes.h \
$(SRC_DIR)/ClassMinerOptions.h \
$(SRC_DIR)/DataStructures.h \
$(SRC_DIR)/Dataset.h \
$(SRC_DIR)/Discretize.h \
$(SRC_DIR)/Parser.h \
$(SRC_DIR)/ParseOptions.h \
$(SRC_DIR)/Rawdata.h \
$(SRC_DIR)/graph.h \
$(SRC_DIR)/spear.h

OBJ += \
$(ODIR)/Attributes.o \
$(ODIR)/ClassMinerOptions.o \
$(ODIR)/DataStructures.o \
$(ODIR)/Dataset.o \
$(ODIR)/Discretize.o \
$(ODIR)/Parser.o \
$(ODIR)/ParseOptions.o \
$(ODIR)/Rawdata.o \
$(ODIR)/graph.o \
$(ODIR)/spear.o 

LIBS= -lCarmen -lboost_system -lboost_regex -lboost_thread
LIBS_DIR := -L bin -L /usr/local/lib 

all: carmen spear

$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -o"$@" "$<" $(CFLAGS) -c

carmen: $(OBJ)
	@echo 'Building target: $@'
	chmod 744 obj/*
	ar -r "bin/libCarmen.a" $(OBJ)
	chmod 744 bin/libCarmen.a
	@echo 'Finished building target: $@'
	@echo ' '

spear: carmen
	@echo 'Building target: spear'
	$(CC) $(CFLAGS) src/SpearMain.cpp $(LIBS_DIR) $(LIBS) -obin/spear
	chmod 744 bin/spear
	@echo 'Finished building target: spear'
	@echo ' '

clean:
	-$(RM) $(OBJ) bin/libCarmen.a bin/spear
	-@echo ' '
