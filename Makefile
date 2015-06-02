CC = gcc
INC_DIR = include
SRC_DIR = src
LIB_DIR = lib
TEST_DIR = tests
EXE_DIR = executables

CC_ARGS = -std=c99 $(INC)
GMP = -lgmp
OMP = -fopenmp
USED_LIBS = $(GMP) $(OMP)
INC = -I$(INC_DIR)

$(EXE_DIR)/gaussian_elimination-con-main: $(LIB_DIR)/sieve.o $(LIB_DIR)/vector.o $(LIB_DIR)/matrix.o $(TEST_DIR)/gaussian_elimination-con-main.c
	$(CC) $^ $(CC_ARGS) $(USED_LIBS) -o $@

$(LIB_DIR)/eratostene.o: $(SRC_DIR)/eratostene.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/factor_base.o: $(SRC_DIR)/factor_base.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/gaussian_elimination.o: $(SRC_DIR)/gaussian_elimination.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/legendre.o: $(SRC_DIR)/legendre.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/matrix.o: $(SRC_DIR)/matrix.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/sieve.o: $(SRC_DIR)/sieve.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/vector.o: $(SRC_DIR)/vector.c
	$(CC) $^ $(CC_ARGS) -c -o $@

clean:
	-rm lib/*