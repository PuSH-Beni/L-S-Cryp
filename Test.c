//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//	This file is a test file, to test the functions implemented in the project.



#include "LSbox.h"
#include <stdio.h>
#include <avr/io.h>
#include <avr/interrupt.h>

/* <time.h> must be declared explicitly in C99 */
//#include <time.h>


/*=========================================================*/
/* MARK:  Toggle Of Test Options     */
/*=========================================================*/
#define CONSTRUCT_MAT_TEST				 1
#define FILE_IO_TEST					 0
#define PIN_IO							 1
#define ENCRYPT_TEST					 1
#define TIMES							 2

#define PERIOD_TIME						 1 // ( 1 / FCPU * TCCR0  ) us , where FCPU should be some Mhz

static volatile
unsigned int totOverflows;

#define uc unsigned char

static
const uc digits[]={0xc0,0xf9,0xa4,0xb0,0x99,0x92,0x82,0xf8,0x80,0x90};
//uc digits[10]={0x3F,0x06,0x5B,0x4F,0x66,0x6D,0x7D,0x07,0x7F,0x6F};

static
uc ledBuff[4] = {
	0xFF, 0xFF, 0xFF, 0xFF
};

void delayus(int us)
{
	  int i;
	  for( i=0;i<us;i++);
}

void delayms(int ms)
{
	while(ms--){
		delayus(1000);
	}
}

void show(unsigned int toShow)
{
	DDRG |= 0x0F;
	DDRD  = 0XFF;
	
	int ledIndex = 0;
	ledBuff[0] = digits[toShow%10];
	ledBuff[1] = digits[toShow/10%10];
	ledBuff[2] = digits[toShow/100%10];
	ledBuff[3] = digits[toShow/1000%10];
	while(1){
		delayms(1);
		switch (ledIndex){
			case 0: PORTG |= 0x01; PORTD=ledBuff[0]; break;
			case 1: PORTG |= 0x02; PORTD=ledBuff[1]; break;
			case 2: PORTG |= 0x03; PORTD=ledBuff[2]; break;
			case 3: PORTG |= 0x04; PORTD=ledBuff[3]; break;
			default:break;
		}
		++ledIndex;
		if(ledIndex > 3) ledIndex = 0;
	}
	
}

/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
int main(){

	
	/*=====================================================*/
	/* =============    TEST BEGINS  =================== */
	/*=====================================================*/
	
	#if CONSTRUCT_MAT_TEST
	Mat *matX = newMat(DIM_S, DIM_L, NULL, 0x00);
	Mat *matY = newMat(DIM_S, DIM_L, NULL, 0x00);
	#endif

	#if FILE_IO_TEST
	
	FILE *fin = (FILE*)fopen("in.txt", "r");
	FILE *fout = (FILE*)fopen("out.txt", "a");
	srand((WORD)time(NULL));
	
	/* Dectect the input file existed or not */
	if (fin == NULL || fout == NULL) {
		printf( "File Doesnt Exist\n");
		return 1;
	}
	else{
		printf(  "File Read Seccessfully\n");
	}

	/* Read the first matrix */
	int i;
	for (i = 0; i < DIM_L; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matX->vect + i) = (BASE)tem;
	}
	/* Show it */
	printf( "\n ==>Mat X :\n");
	for (i = 0; i < DIM_L; ++i) {
		printf( "%x ", *(matX->vect + i));
	}
	printf( "\n");

	/* Read the second matrix */
	for (i = 0; i < DIM_L; ++i) {
		int tem;
		fscanf(fin, "%x ", &tem);
		*(matY->vect + i) = (BASE)tem;
	}
	/* Show it */
	printf( "\n ==>Mat Y  :\n");
	for (i = 0; i < DIM_L; ++i) {
		printf( "%x ", *(matY->vect + i));
	}
	printf( "\n");
	#endif

	#if PIN_IO
	DDRA = 0x00;
	DDRB = 0xFF;
	int IOIndex;
	for (IOIndex = 0; IOIndex < DIM_L; ++IOIndex) {
		*(matX->vect + IOIndex) = (BASE)0x1b;
	}
	for (IOIndex = 0; IOIndex < DIM_L; ++IOIndex) {
		*(matY->vect + IOIndex) = (BASE)0x2c;
	}
	/* Show it */
	for (IOIndex = 0; IOIndex < DIM_L; ++IOIndex) {
		PORTB =  *(matX->vect + IOIndex);
	}
	for (IOIndex = 0; IOIndex < DIM_L; ++IOIndex) {
		PORTB =  *(matY->vect + IOIndex);
	}
	#endif


	#if ENCRYPT_TEST
	/* MARK: ENCRYPTO */

	/* L-S-Model Encryption */

	Mat *cipher;
	//printf( "\n ---> Encrytion: \n");

	#if DIM_A
	setupEn();
	#endif
	newPreCal();

	int repeatIndex;
	//double time_Start = (double)clock();

	
	// Pre-scaler = FCPU/8(0b010) , 64(0b011), 256(0b100), 1024(0b101)
	TCCR1B |= (1<<CS11)|(0<<CS10);
	//Enable Overflow Interrupt Enable
	TIMSK |=(1<<TOIE1);
	//Enable Global Interrupts
	sei();
	totOverflows = 0;
	TCNT1 = 0;
	show(0x0001);
	
	for (repeatIndex = 0; repeatIndex != TIMES; ++repeatIndex){
		cipher = encrypto(matX, matY);

		//if (repeatIndex < 2){
			////printf("\n==>LSout:\n");
			//int byteIndex;
			//for (byteIndex = 0; byteIndex < DIM_L; ++byteIndex) {
				//PORTB =  *(cipher->vect + byteIndex);
			//}
		//}
		//deMat(cipher);
		
	}
	
	unsigned int timeCost = TCNT1 * PERIOD_TIME;
	show(timeCost);
	
	//double time_End = (double)clock();
	dePostCal();
	//printf("\n ---> Time Cost: \n");

	//printf("%.fms\n", (time_End - time_Start));
	#endif


	/*=====================================================*/
	/* MARK: Free the allocated memories */
	/*=====================================================*/
	#if CONSTRUCT_MAT_TEST
	deMat(matX);
	deMat(matY);
	#endif

	return 0;
}
ISR (TIMER0_OVF_vect)  // timer0 overflow interrupt
{
	++totOverflows;
}
