raytracer: rt05w.c
	x86_64-w64-mingw32-gcc rt05w.c -lSDL2 -o rt05.exe -lm 

glitch: glitch.c
	gcc glitch.c -lSDL2 -o glitch -lm


