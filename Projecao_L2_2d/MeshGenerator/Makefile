CC = gcc
CFLAGS = -Wall -O2 
PROGRAM_NAME = meshGenerator

all: $(PROGRAM_NAME)

$(PROGRAM_NAME): $(PROGRAM_NAME).o
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(PROGRAM_NAME).o: $(PROGRAM_NAME).c
	$(CC) $(CFLAGS) -o $@ -c $< -lm

clean:
	rm *.o $(PROGRAM_NAME)

print-%  : ; @echo $* = $($*)
