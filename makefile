TARGET = main

.PHONY: all clean print

CPP = g++

S_DIR = ./source

H_DIR = ./headers

O_DIR = ./object_files

SRCS = $(shell find $(S_DIR) | grep ".cpp" | head -n -1)

OBJS = $(SRCS:$(S_DIR)%.cpp=$(O_DIR)%.o)

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
