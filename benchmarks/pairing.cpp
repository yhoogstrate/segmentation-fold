#include <iostream>
#include <stdio.h>

// benchmark to see which kind of Pairing initiazation is fastest

// double switch - 13.19s
char init_old(char i, char j)
{
	switch(i)
	{
		case 0:
			switch(j)
			{
				case 3:
					return 1;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case 1:
			switch(j)
			{
				case 2:
					return 2;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case 2:
			switch(j)
			{
				case 3:
					return 3;
				case 1:
					return 4;
				default:														// Makes the compiler happy..
					break;
			}
			break;
			
		case 3:
			switch(j)
			{
				case 0:
					return 5;
				case 2:
					return 6;
				default:														// Makes the compiler happy..
					break;
			}
			break;
	}
	
	return 0;
}

// multiply - 18.00s
char init_mp(char i, char j)
{
	char mp = (i*4) + j;
	switch(mp)
	{
		case 3:
			return 1;
		case 6:
			return 2;
		case 9:
			return 4;
		case 11:
			return 3;
		case 12:
			return 5;
		case 14:
			return 6;
		default:
			return 0;
	}
}

// bitshift - 18.86
char init_bs(char i, char j)
{
	char bs = (i << 2) + j;
	switch(bs)
	{
		case 3:
			return 1;
		case 6:
			return 2;
		case 9:
			return 4;
		case 11:
			return 3;
		case 12:
			return 5;
		case 14:
			return 6;
		default:
			return 0;
	}
}

int main()
{
	unsigned int n_tests = 100000000;
	
	char product;
	
	for(unsigned int k = 0; k < n_tests; k++)
	{
		for(char i = 0; i< 4; i++)
		{
			for(char j = 0; j < 4; j++)
			{
				//product = init_old(i, j);
				//product = init_mp(i, j);
				product = init_bs(i, j);
				
				if(product > 6)
				{
					//printf("%i,%i:  %i  %i  %i \n",i,j, product_a, product_b, product_c);
					throw std::invalid_argument("incorrect product");
				}
				/*
				if(product_b != product_c)
				{
					//printf("%i,%i:  %i  %i  %i \n",i,j, product_a, product_b, product_c);
					throw std::invalid_argument("incorrect product c");
				}
				*/
			}
		}
	}
}