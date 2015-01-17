PROJECT = matija_santl_dz1.exe
SOURCE = matija_santl_dz1.c
HEADERS = 

CC = gcc
CFLAGS = -Wall -g 
OBJECTS = ${SOURCE:.c=.o}

$(PROJECT): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(PROJECT)

$(OBJECTS): $(HEADERS)

clean:
	-rm -f $(PROJECT) $(OBJECTS) *.core
