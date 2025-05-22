CC = mpicc
CFLAGS = -O3 -Wall
LDFLAGS = -lm

TARGET = shearsort
SRC = main.c source_files/shearsort.c
HEADERS = header_files/shearsort.h

all: $(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)
