EXEC = projet
SRC = projet.c
OBJ = $(SRC:.c=.o)
CC = gcc
CFLAGS = -Wall -Wextra -O2

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

run: $(EXEC)
	@echo "Running the program with redirection to solution.mtx..."
	@./$(EXEC) inputMatrix.mtx inputRhs.mtx 1e-6 > solution.mtx

clean:
	rm -f $(OBJ) $(EXEC)

mrproper: clean
	rm -f solution.mtx
