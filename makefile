TARGET = main

.PHONY: all install clean print

CPP = g++

S_DIR = ./source

H_DIR = ./headers

O_DIR = ./object_files

D_DIR = ./data

.SILENT: print install $(O_DIR) $(D_DIR)

SRCS = $(shell ls $(S_DIR) | grep ".cpp") 

OBJS = $(SRCS:%.cpp=$(O_DIR)/%.o)

INC_FLAGS = $(addprefix -I,$(H_DIR))

CFLAGS = $(INC_FLAGS) -O0 -Wall

all: $(O_DIR) $(D_DIR) $(TARGET) 

$(TARGET): $(OBJS)
	$(CPP) $(OBJS) -o $@

$(O_DIR)/%.o: $(S_DIR)/%.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

$(O_DIR):
	if [ -d "$(O_DIR)" ]; then \
		echo "$(O_DIR) exists"; \
	else \
		mkdir $(O_DIR); \
		echo "$(O_DIR) created"; \
	fi

$(D_DIR):
	if [ -d "$(D_DIR)" ]; then \
		echo "$(D_DIR) exists";\
	else \
		mkdir $(D_DIR); \
		echo "$(D_DIR) created"; \
	fi

install:
	if [ -d "$(O_DIR)" ]; then \
		echo "$(O_DIR) exists"; \
	else \
		mkdir $(O_DIR); \
		echo "$(O_DIR) created"; \
	fi
	if [ -d "$(D_DIR)" ]; then \
		echo "$(D_DIR) exists";\
	else \
		mkdir $(D_DIR); \
		echo "$(D_DIR) created"; \
	fi

clean:
	rm -rf $(O_DIR)/* $(TARGET) $(D_DIR)/*

print:
	echo $(SRCS)
