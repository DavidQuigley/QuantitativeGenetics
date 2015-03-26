# This should be converted to a proper makefile at a later date

DIR_CARMEN=/notebook/code/src/carmen
DIR_EQTL=/notebook/code/src/eqtl
g++ -o $DIR_CARMEN/build/obj/Attributes.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/Attributes.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/ClassMinerOptions.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/ClassMinerOptions.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/Dataset.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/Dataset.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/DataStructures.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/DataStructures.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/Discretize.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/Discretize.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/graph.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/graph.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/ParseOptions.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/ParseOptions.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/Parser.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/Parser.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/Rawdata.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/Rawdata.cpp -O3 -Wall -pthread 
g++ -o $DIR_CARMEN/build/obj/spear.o -I/usr/local/include -I$DIR_CARMEN/src -c $DIR_CARMEN/src/spear.cpp -O3 -Wall -pthread 
chmod 755 $DIR_CARMEN/build/obj/*.o

g++ -o $DIR_EQTL/obj/QTL_Calculator.o \
 -I/usr/local/include -I$DIR_CARMEN/src -I$DIR_EQTL/src  \
 -c $DIR_EQTL/src/QTL_Calculator.cpp \
 -O3 -Wall -pthread 

g++ -o $DIR_EQTL/obj/eQTL.o \
 -I/usr/local/include -I$DIR_CARMEN/src -I$DIR_EQTL/src \
 -c $DIR_EQTL/src/eQTL.cpp \
 -O3 -Wall -pthread 

g++ -o "eqtl" \
  -pthread \
  $DIR_CARMEN/build/obj/Attributes.o \
  $DIR_CARMEN/build/obj/ClassMinerOptions.o \
  $DIR_CARMEN/build/obj/Dataset.o \
  $DIR_CARMEN/build/obj/DataStructures.o \
  $DIR_CARMEN/build/obj/Discretize.o \
  $DIR_CARMEN/build/obj/graph.o \
  $DIR_CARMEN/build/obj/ParseOptions.o \
  $DIR_CARMEN/build/obj/Parser.o \
  $DIR_CARMEN/build/obj/Rawdata.o \
  $DIR_CARMEN/build/obj/spear.o \
  $DIR_EQTL/obj/QTL_Calculator.o $DIR_EQTL/obj/eQTL.o \
  -L/usr/local/lib -lRmath -lboost_system -lboost_regex -lboost_thread -lboost_filesystem 
