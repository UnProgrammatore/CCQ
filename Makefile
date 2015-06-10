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

#main1: $(EXE_DIR)/qs

$(EXE_DIR)/qs: $(LIB_DIR)/trivial_fact.o $(LIB_DIR)/eratostene.o $(LIB_DIR)/smart_sieve.o $(LIB_DIR)/vector.o $(LIB_DIR)/matrix.o $(LIB_DIR)/linear_algebra.o $(LIB_DIR)/base_fattori.o $(LIB_DIR)/quadratic_sieve.o $(SRC_DIR)/main-qs.c
	$(CC) $^ $(CC_ARGS) $(USED_LIBS) -o $@

$(EXE_DIR)/main-prova: $(LIB_DIR)/eratostene.o $(LIB_DIR)/sieve.o $(LIB_DIR)/vector.o $(LIB_DIR)/matrix.o $(LIB_DIR)/linear_algebra.o $(LIB_DIR)/base_fattori.o $(SRC_DIR)/main-prova.c
	$(CC) $^ $(CC_ARGS) $(USED_LIBS) -o $@

$(EXE_DIR)/gaussian_elimination-con-main: $(LIB_DIR)/sieve.o $(LIB_DIR)/vector.o $(LIB_DIR)/matrix.o $(TEST_DIR)/gaussian_elimination-con-main.c
	$(CC) $^ $(CC_ARGS) $(USED_LIBS) -o $@

$(EXE_DIR)/prova_le_altre_matrix: $(LIB_DIR)/matrix.o $(SRC_DIR)/prova_le_altre_matrix.c
	$(CC) $^ $(CC_ARGS) $(USED_LIBS) -o $@

$(LIB_DIR)/trivial_fact.o: $(SRC_DIR)/trivial_fact.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/quadratic_sieve.o: $(SRC_DIR)/quadratic_sieve.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/eratostene.o: $(SRC_DIR)/eratostene.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/factor_base.o: $(SRC_DIR)/factor_base.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/base_fattori.o: $(SRC_DIR)/base_fattori.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/linear_algebra.o: $(SRC_DIR)/linear_algebra.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/legendre.o: $(SRC_DIR)/legendre.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/matrix.o: $(SRC_DIR)/matrix.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/smart_sieve.o: $(SRC_DIR)/smart_sieve.c
	$(CC) $^ $(CC_ARGS) -c -o $@

$(LIB_DIR)/vector.o: $(SRC_DIR)/vector.c
	$(CC) $^ $(CC_ARGS) -c -o $@

clean:
	-rm lib/*
