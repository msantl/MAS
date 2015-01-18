PROJECT = matija_santl_dz1.exe
SOURCE = matija_santl_dz1.c
HEADERS = 

CC = gcc
CFLAGS = -Wall -g 
LIBS = -lm
OBJECTS = ${SOURCE:.c=.o}

$(PROJECT): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(PROJECT) $(LIBS)

$(OBJECTS): $(HEADERS)

clean:
	-rm -f $(PROJECT) $(OBJECTS) *.core
