//
//  Test.c
//  L-S-model
//
//  Created by benny on 16/1/5.
//
//	This file is a test file, to test the functions implemented in the project.



#include "LSbox.h"
#include <stdio.h>

//#include <util/delay.h>

/* <time.h> must be declared explicitly in C99 */
//#include <time.h>


/*=========================================================*/
/* MARK:  Toggle Of Test Options     */
/*=========================================================*/
#define CONSTRUCT_MAT_TEST				 0
#define FILE_IO_TEST					 0
#define PIN_IO							 0
#define ENCRYPT_TEST					 1
#define TIMES							 2

#define PERIOD_TIME						 1 // ( 1 / FCPU * TCCR0  ) us , where FCPU should be some Mhz

static volatile
unsigned int totOverflows;

#define uc unsigned char
//#define F_CPU 16000000UL

//static
//const uc digits[]={0xc0,0xf9,0xa4,0xb0,0x99,0x92,0x82,0xf8,0x80,0x90};
uc tab[10]={0x3F,0x06,0x5B,0x4F,0x66,0x6D,0x7D,0x07,0x7F,0x6F};


const int K = 577;

#define SEG_PORT PORTD

#define SEG4_ON   PORTG &= 0XFE;
#define SEG4_OFF  PORTG |= 0X01;

#define SEG3_ON   PORTG &= 0XFD;
#define SEG3_OFF  PORTG |= 0X02;

#define SEG2_ON   PORTG &= 0XF7;
#define SEG2_OFF  PORTG |= 0X08;

#define SEG1_ON   PORTG &= 0XEF;
#define SEG1_OFF  PORTG |= 0X10;

void delay_us(unsigned int microsecond)
{
	do
	{
		microsecond--;
	}
	while (microsecond>1);
}

void delay_ms(unsigned int millisecond)
{
	while (millisecond--)
	{
		delay_us(999);
	}
}

void delay(void)
{
	int i;
	for(i=0;i<2000;i++){
		delay_ms(1000);
	}
}
void port_init(void)
{

	PORTD = 0x00;
	DDRD  = 0xFF;
	
	PORTG = 0x1B;
	DDRG  = 0x1B;
	
	//DDRA = 0x00;
	//DDRB = 0xFF;
	//PORTB = 0x00;
}

void init_devices(void)
{
	//stop errant interrupts until set up
	cli(); //disable all interrupts
	XDIV  = 0x00; //xtal divider
	XMCRA = 0x00; //external memory
	port_init();

	MCUCR = 0x00;
	EICRA = 0x00; //extended ext ints
	EICRB = 0x00; //extended ext ints
	EIMSK = 0x00;
	TIMSK = 0x00; //timer interrupt sources
	sei(); //re-enable interrupts
	//all peripherals are now initialized
}

void show(unsigned int n)
{
	SEG_PORT = tab[n/1000];
	SEG1_ON;
	delay_ms(K);
	SEG1_OFF;
	SEG_PORT = tab[(n%1000)/100];
	SEG2_ON;
	delay_ms(K);
	SEG2_OFF;
	SEG_PORT = tab[(n%100)/10];
	SEG3_ON;
	delay_ms(K);
	SEG3_OFF;
	SEG_PORT = tab[n%10];
	SEG4_ON;
	delay_ms(K);
	SEG4_OFF;
}



/*=========================================================*/
/* MARK: Main Function Begins    */
/*=========================================================*/
int main(){

	
	/*=====================================================*/
	/* =============    TEST BEGINS  =================== */
	/*=====================================================*/
	init_devices();
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

	//#if DIM_A
	//setupEn();
	//#endif
	//newPreCal();

	//int repeatIndex;
	//double time_Start = (double)clock();

	while (1)
	{
		// Pre-scaler = FCPU/8(0b010) , 64(0b011), 256(0b100), 1024(0b101)
		//TCCR1B |= (1<<CS11)|(0<<CS10);
		//Enable Overflow Interrupt Enable
		//TIMSK |=(1<<TOIE1);
		//Enable Global Interrupts
		//totOverflows = 0;
		//TCNT1 = 0;
		
		show(1234);
		
		int m;
		for(m = 0; m < 2000; ++m){
			show(1234);
			delay_ms(1000);
		}
		for(m = 0; m < 2000; ++m){
			show(5678);
			delay_ms(1000);
		}
		/*for (repeatIndex = 0; repeatIndex != TIMES; ++repeatIndex){
		cipher = encrypto(matX, matY);

		//if (repeatIndex < 2){
		////printf("\n==>LSout:\n");
		//int byteIndex;
		//for (byteIndex = 0; byteIndex < DIM_L; ++byteIndex) {
		//PORTB =  *(cipher->vect + byteIndex);
		//}
		//}
		//deMat(cipher);
		
		}*/
		
		//unsigned int timeCost = TCNT1 * PERIOD_TIME;
		//show(5678);
	}
	
	//double time_End = (double)clock();
	//dePostCal();
	//printf("\n ---> Time Cost: \n");
	//printf("%.fms\n", (time_End - time_Start));
	#endif


	/*=====================================================*/
	/* MARK: Free the allocated memories */
	/*=====================================================*/
	//#if CONSTRUCT_MAT_TEST
	//deMat(matX);
	//deMat(matY);
	//#endif
	return 0;
}
//ISR (TIMER0_OVF_vect)  // timer0 overflow interrupt
//{
//++totOverflows;
//}
