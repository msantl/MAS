PROJECT = jpeg.exe
SOURCE = jpeg.c
HEADERS = 

CC = gcc
CFLAGS = -Wall -g 
OBJECTS = ${SOURCE:.c=.o}

$(PROJECT): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(PROJECT)

$(OBJECTS): $(HEADERS)

clean:
	-rm -f $(PROJECT) $(OBJECTS) *.core
