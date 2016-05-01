//WIN-GCC application builder : 2015/12/13 20:52:33
/*******************************************************************************
* 版权:     涛哥工作室
*
* 单片机:   ATMAGE128A AU 1036
* 晶振:     外部8MHz
* 编译器:   AVRStduio
*
* 文件名:   main.c
* 版本:     1.0
* 功能描述: 在8M晶振下,实现数码管动态显示操作
*******************************************************************************/
//#include "avr/io.h"
//#include "iom128.h"
#include <string.h>
#include <stdio.h>
#include <avr/io.h>
#include <avr/iom128.h>
//#include <avr/signal.h>     //中断信号头文件
#include <avr/interrupt.h>

#define LED PORTF
#define SEG_PORT PORTD     //数码管数据口

#define SEG4_ON   PORTG &= 0XFE;
#define SEG4_OFF  PORTG |= 0X01;

#define SEG3_ON   PORTG &= 0XFD;
#define SEG3_OFF  PORTG |= 0X02;

#define SEG2_ON   PORTG &= 0XF7;
#define SEG2_OFF  PORTG |= 0X08;

#define SEG1_ON   PORTG &= 0XEF;
#define SEG1_OFF  PORTG |= 0X10;

unsigned char key_state=0;
const unsigned char tab[]={0X3F,0X06,0X5B,0X4F,0X66,0X6D,0X7D,0X07,0X7F,0X6F};//数码管
unsigned int i = 0,K = 577,m = 0;
void delay_us(unsigned int microsecond)
{
	do
	{
		microsecond--;
	}
	while (microsecond>1);
}

void delay(void)
{
	int i;
	for(i=0;i<200;i++);
}
void port_init(void)
{

	PORTD = 0x00;
	DDRD  = 0xFF;
	PORTG = 0x1B;
	DDRG  = 0x1B;
}
void SEG(unsigned int n)
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
//call this routine to initialize all peripherals
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

void delay_ms(unsigned int millisecond)
{
	while (millisecond--)
	{
		delay_us(999);
	}
}
int main(void)
{
	init_devices();

	while(1)
	{
		for(i=0;i<14;i++)
		{
			for(m=0;m<2000;m++)
			{
				if(i%2)SEG(1234);
				else SEG(5678);
				delay_ms(1000);
			}
			m=0;
		}
	}
}

