include ../Machine

all: moduledep.x dep

moduledep.x : moduledep.c
	$(CC) -o moduledep.x -O2 moduledep.c

dep: moduledep.x
	( cd ../ ; ./shdep )
