CC = gcc -std=gnu99
CFLAGS = -Wall -Wextra --pedantic -O3
BUILD_DIR = build
LIB_DIR = lib

all : dirs $(LIB_DIR)/libOpenSimplex2.a

$(LIB_DIR)/libOpenSimplex2.a : $(BUILD_DIR)/OpenSimplex2.o $(BUILD_DIR)/OpenSimplex2S.o
	ar rcs \
		$(LIB_DIR)/libOpenSimplex2.a \
		$(BUILD_DIR)/OpenSimplex2.o \
		$(BUILD_DIR)/OpenSimplex2S.o

$(BUILD_DIR)/OpenSimplex2.o : OpenSimplex2.c _OpenSimplex2.h
	$(CC) $(CFLAGS) -o $(BUILD_DIR)/OpenSimplex2.o -c OpenSimplex2.c

$(BUILD_DIR)/OpenSimplex2S.o : OpenSimplex2S.c _OpenSimplex2S.h
	$(CC) $(CFLAGS) -o $(BUILD_DIR)/OpenSimplex2S.o -c OpenSimplex2S.c

utils : dirs $(BUILD_DIR)/MakeGradients
	$(BUILD_DIR)/MakeGradients

$(BUILD_DIR)/MakeGradients : $(BUILD_DIR)/MakeGradients.o
	$(CC) $(CFLAGS) -o $(BUILD_DIR)/MakeGradients $(BUILD_DIR)/MakeGradients.o

$(BUILD_DIR)/MakeGradients.o : MakeGradients.c _OpenSimplex2.h _OpenSimplex2S.h
	$(CC) $(CFLAGS) -o $(BUILD_DIR)/MakeGradients.o -c MakeGradients.c

dirs :
	mkdir -p $(BUILD_DIR) $(LIB_DIR) gradients/fast gradients/smooth

clean :
	rm -rf $(BUILD_DIR) $(LIB_DIR)
