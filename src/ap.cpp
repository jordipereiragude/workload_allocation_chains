/*
 * Copyright (c) 2023 Jordi Pereira, Marcus Ritt
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "ap.hpp"

using namespace std;

//Assignment problem code translated from pascal code from Jonkers and Volgenat
// ``A shortest augmenting path algorithm for dense and sparse linear assignment problems''
// Computing volume 38, pages 325–340 (1987)
// The original code is available at http://www.assignmentproblems.com/LAPJVBOT.htm

namespace ap { //assignment problem code

int inf;
//arrays and vectors numbered from 1 to n 
int** c;
int* x;
int* y;
 
void freeAP(int n) {
	for(int i=0;i<=n;i++)
	{
		free(c[i]);
	}
	free(c);
  free(x);
  free(y);
  return;
}

void initAP(int n) {
  
  inf=RAND_MAX;
	c=(int**)malloc((n+1)*sizeof(int *));
	for(int i=0;i<=n;i++)
	{
		c[i]=(int *) malloc((n+1)*sizeof(int));
	}
	x=(int *) malloc((n+1)*sizeof(int));
	y=(int *) malloc((n+1)*sizeof(int));
  return;
}

int bottleneckAP(int n) {
  int h,i,j,k,i0=0,j0,dj,nn,min,sol,low,up;
  int col[n+1];
  int lab[n+1];
  int d[n+1];

	sol=0-inf;
	nn=0;
	for(i=1;i<=n;i++) x[i]=0;
	for(j=n;j>=1;j--)
	{
		col[j]=j;
		min=inf;
		i=1;
		do
		{
			dj=c[i][j];
			if(dj<=sol)
			{
				min=sol;
				i0=i;
				i=n+1;
			}
			else
			{
				if(dj<min)
				{
					min=c[i][j];
					i0=i;
				}
				i=i+1;
			}
		}while(i<=n);
		if(x[i0]==0)
		{
			nn=nn+1;
			x[i0]=j;
			y[j]=i0;
		}
		else
		{
			y[j]=0;
		}
		if(min>sol) sol=min;
	}
	for(i0=1;i0<=n;i0++)
	{
		if(x[i0]==0)
		{
			up=1;
			low=1;
			for(k=1;k<=n;k++)
			{
				j=col[k];
				d[j]=c[i0][j];
				lab[j]=i0;
				if(d[j]<=sol)
				{
					if(y[j]==0) goto label1;
					else
					{
						col[k]=col[up];
						col[up]=j;
						up=up+1;
					}
				}
			}
			do
			{
				if(up==low)
				{
					min=inf;
					for(k=up;k<=n;k++)
					{
						j=col[k];
						dj=d[j];
						if(dj<=min)
						{
							if(dj<min)
							{
								up=low;
								min=dj;
							}
							col[k]=col[up];
							col[up]=j;
							up++;
						}
  				}
  				sol=min;
  				for(h=low;h<=(up-1);h++)
  				{
  					j=col[h];
  					if(y[j]==0) goto label1;
  				}
  			}
  			j0=col[low];
  			low++;
  			i=y[j0];
  			for(k=up;k<=n;k++)
  			{
  				j=col[k];
  				dj=c[i][j];
  				if(dj<d[j])
  				{
  					lab[j]=i;
  					if(dj<=sol)
  					{
  						if(y[j]==0) goto label1;
  						else
  						{
  							col[k]=col[up];
  							col[up]=j;
  							up++;
  						}
  					}
  					d[j]=dj;
  				}
  			}
  		}while(1);
label1:
			i=lab[j];
			y[j]=i;
			k=j;
			j=x[i];
			x[i]=k;
			for(;i!=i0;)
			{
				i=lab[j];
				y[j]=i;
				k=j;
				j=x[i];
				x[i]=k;
			}
		}
	}
	return(sol);
}

}
