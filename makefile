PRG = Test
OBJ = Test.o, Fundamentals.o, Lsbox.o
#MCU_TARGET = atmega165
MCU_TARGET = atmega3250
target:
	avr-gcc -g -Os -mmcu=atmega3250 -c Test.c
	avr-gcc -g -Os -mmcu=atmega3250 -c Fundamentals.c
	avr-gcc -g -Os -mmcu=atmega3250 -c LSbox.c
	avr-gcc -g -mmcu=atmega3250 -o Test.elf Test.o LSbox.o Fundamentals.o
run:
	
clean:
	rm *.o
	rm *.elf
	rm *.map

