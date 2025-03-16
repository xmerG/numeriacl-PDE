CC = g++
CFLAG = -std=c++11 -O2
LDFLAGS = -llapacke -lblas
SRC = .
OD = ..
OBJ = $(SRC)/Circle.o	#	object files of all common files
CIRCLEOBJ = $(SRC)/test1.o			#	object file for pdadd

all: $(OD)/pdadd $(OD)/pdlist $(OD)/pdremove $(OD)/pdshow

	

.cpp.o:
	$(CC) $(CFLAG) -c $< -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(OBJ) $(ADDOBJ) $(LISTOBJ) $(REMOVEOBJ) $(SHOWOBJ) $(OD)/pdadd $(OD)/pdlist $(OD)/pdremove $(OD)/pdshow
