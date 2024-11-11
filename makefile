SLATEC 	= $(wildcard ./*.cpp)
SLATEC_OBJ = $(addprefix ./,$(notdir $(SLATEC:.cpp=.o)))

FLAGS = -O1
FCOMP = g++

all: main

main: $(SLATEC_OBJ)
	$(FCOMP) $(FLAGS) -o main *.o

%.o: %.cpp
	$(FCOMP) $(FLAGS) -c $<

clean:
	rm -f *.o