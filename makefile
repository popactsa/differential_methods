TARGET = main

.PHONY: all clean print

CPP = g++

S_DIR = ./source

H_DIR = ./headers

O_DIR = ./object_files

SRCS = $(shell ls $(S_DIR) | grep ".cpp") 

OBJS = $(SRCS:%.cpp=$(O_DIR)/%.o)

INC_FLAGS = $(addprefix -I,$(H_DIR))

CFLAGS = $(INC_FLAGS) -O1

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CPP) $(OBJS) -o $@

$(O_DIR)/%.o: $(S_DIR)/%.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(O_DIR)/* $(TARGET)

print:
	echo $(SRCS)
