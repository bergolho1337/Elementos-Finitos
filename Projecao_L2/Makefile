INCLUDE_DIR = include
SOURCE_DIR = src
OBJECT_DIR = obj
CC = gcc
CFLAGS = -Wall -O2 -I$(INCLUDE_DIR)
CC_SOURCES = $(wildcard src/*)
STRING_OBJ_AUX = $(CC_SOURCES:.c=.o)
STRING_OBJ = $(subst src/,,$(STRING_OBJ_AUX))
CC_OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(STRING_OBJ))
PROGRAM_NAME = projection_L2

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(CC_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $< -lm

clean:
	rm $(OBJECT_DIR)/*.o *.dat $(PROGRAM_NAME)

print-%  : ; @echo $* = $($*)
