#include<iostream>
#include<cmath>
using namespace std;
void printmatrix(double *A,int row,int col)
{
	for(int i=0;i<row;i++) 
	{
 		cout<<"[ ";
 		for(int j=0;j<col;j++)
 			{
 				cout<<" "<<A[i*col+j];
 			}
 		cout<<" ]"<<endl;
 	}
}
void scalarmul(double* r,double *b,double s,int row)
{
 	int y=0; 
 	for(int i=0;i<row;i++)
 	{
  		r[y++] = b[i]*s;
 	} 
}
void sub(double* a,double* b,double* c,int row)
{
 	for(int i=0;i<row;i++)
 	{
  		c[i]=a[i]-b[i];
 	}
}
void add(double* a,double* b,double* c,int row)
{
 	for(int i=0;i<row;i++)
 	{
  		c[i]=a[i]+b[i];
 	}
}
void extraction(double *A,double *b,int row,int col,int p)
{
 	int t=0;
 	for(int i=0;i<row;i++)
 	{
  		b[t++]=*(A+i*col+p);
 
 	}
}
double dotproduct(double* v1,double* v2,int row)
{
 	double x=0;
 	for(int i=0;i<row;i++)
	{
  		x+=v1[i]*v2[i];
	}
 	return x;
}
double norm(double *b,int row)
{
 	double n;
	n=sqrt(dotproduct(b,b,row));
 	return n;
}
void vectorcap(double* v,double *b,double n,int row)
{
 	int s=0; 
 	for(int i=0;i<row;i++)
 	{
 	 	v[s++] = b[i]/n;
 	}
}
void transpose(double *A,double *B,int row,int col)
{
	for(int j=0;j<col;j++)
 	{ 
  		for(int i=0;i<row;i++)
  		{
   			B[j*row+i]=A[i*col+j];
  		}
  
 	}
}
void multiply(double *C,double *A,double *B,int row1,int col1,int row2,int col2)
{
   if(col1==row2)
  {
   for(int i=0;i<row1;i++)
     {
	 for(int j=0;j<col2;j++)
	  {
	    C[i*col2+j]=0;
	    for(int k=0;k<col1;k++)
		{
		 C[i*col2+j]+=A[i*col1+k]*B[k*col2+j];
		}
	  }
      }
  }
}
void QR(double *Q,double *R,double *A, int row,int col)
{
 double n;
 double *b;
 double *v;
 double *q;
 double *r;
 double *c;
 double *z;
 b=new double[row];
 v=new double[row];
 q=new double[row];
 r=new double[row];
 c=new double[row];
 z=new double[row];
 extraction(A,b,row,col,0);
 n=norm(b,row);
 vectorcap(v,b,n,row);

//################### Q formation
   for(int i=0;i<row;i++)
    {
    	for(int j=0;j<col;j++)
        {
 		if(j==0)
 		 {
  			 Q[i*col+j]=v[i];
  		 }
  		else
		{
   			 Q[i*col+j]=0;
 		}
 	}
    }
for(int i=1;i<col;i++)
{
    double *t;
    t=new double[row];
        for (int w=0;w<row;w++)
        {
 		t[w]=0;
	}

 	for(int j=0;j<i;j++)
 	{
 		 extraction(A,b,row,col,i);
 		 extraction(Q,q,row,col,j);
 		 double x=dotproduct(b,q,row);
  		 scalarmul(r,q,x,row);
  		 add(t,r,t,row);	
 	}
  extraction(A,b,row,col,i);
  sub(b,t,c,row);
  double e=norm(c,row);
  vectorcap(z,c,e,row);

 	for(int k=0;k<row;k++)
	{
 		 Q[k*col+i] = z[k];
	}
}
//################### R formation;

 double *B;
 B=new double[col*row];
 transpose(Q,B,row,col);
 multiply(R,B,A,col,row,row,col); 
}


