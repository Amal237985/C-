#include<iostream>
#include<cmath>
#include "project.hpp"
#include<fstream>
#include <sstream>
using namespace std;

int main()
{
 
double *L;
int row,col;
cin>>row>>col;

int rows, cols;

rows=416; //change according to your dimension
cols=235; //change according to your dimesion
int rm=rows*cols;


/////////////////////////////////////////////// 
ifstream file("R.csv"); //file name
    if (!file.is_open()) {
        cerr << "Failed to open the file." << endl;
        return 1;
    }
    L=new  double[rm];
    
int numRows = 0;
string line;
while (getline(file, line) && numRows < rows) {
        stringstream lineStream(line);
        string cell;
        int numCols = 0;

        while (getline(lineStream, cell, ',') && numCols < cols) {
            try {
                L[numRows*col+numCols] = stod(cell); // Convert the cell to double
            } catch (const invalid_argument& e) {
                // Handle conversion errors
                L[numRows*col+numCols] = 0.0; // Default value
            }
            numCols++;
        }

        numRows++;
    }

    file.close();





double *u;
u = new double[row*col];
double *v;
v = new double[col*row];

double *A;
A = new double[row*row];
double *Q;
Q = new double[row*row];
double *R;
R = new double[row*row];
double *Q0;
Q0 = new double[row*row];
double *G;
G = new double[row*row];
double *B;
B = new double[col*col];
double *Q1;
Q1 = new double[col*col];
double *R1;
R1 = new double[col*col];
double *Q01;
Q01 = new double[col*col];
double *G1;
G1 = new double[col*col];

//int s=0;
for(int i=0;i<row;i++)
{ 
	for(int j=0;j<col;j++)
	{
        	u[i*col+j]=L[i*col+j];
   	}
}
printmatrix(u,row,col);
transpose(u,v,row,col);
//printmatrix(v,col,row);
multiply(A,u,v,row,col,col,row);
//printmatrix(A,row,row);

QR(Q,R,A,row,row);
for(int i=0;i<row;i++)
{ 
	for(int j=0;j<row;j++)
	{
        	Q0[i*row+j] = Q[i*row+j];
   	}
}
int i=0;
while(i<1000)
{
	multiply(A,R,Q,row,row,row,row);
	QR(Q,R,A,row,row);
	multiply(G,Q0,Q,row,row,row,row);
	for(int i=0;i<row;i++)
	{ 
		for(int j=0;j<row;j++)
		{
        		Q0[i*row+j] = G[i*row+j];
   		}
	}
	i++;
}
cout<<"_________________________________________________to find U"<<endl;
cout<<"eigen values on diagonal...."<<endl;
printmatrix(A,row,row); 
cout<<"________________________________________________________________________"<<endl;
cout<<"correseponding eigen vectors"<<endl;
printmatrix(G,row,row);

multiply(B,v,u,col,row,row,col);
//printmatrix(B,col,col);
QR(Q1,R1,B,col,col) ;
for(int i=0;i<col;i++)
{ 
	for(int j=0;j<col;j++)
	{
        	Q01[i*col+j] = Q1[i*col+j];
   	}
}
int j=0;
while(j<1000000)
{
	multiply(B,R1,Q1,col,col,col,col);
	QR(Q1,R1,B,col,col);
	multiply(G1,Q01,Q1,col,col,col,col);
	for(int i=0;i<col;i++)
	{ 
		for(int j=0;j<col;j++)
		{
        		Q01[i*col+j] = G1[i*col+j];
   		}
	}
	j++;
}
cout<<"____________________________________to find V"<<endl;
//cout<<"eigen values on diagonal......."<<endl;
//printmatrix(B,col,col); 
cout<<"________________________________________________________________________"<<endl;
cout<<" eigen vectors.......V"<<endl;
printmatrix(G1,col,col);
double *G2;
G2 = new double[col*col];

double *S;
S = new double[col*col];
double *N;
N = new double[col*col];
double *P;
P = new double[col*col];
transpose(G1,G2,col,col);

for(int i=0;i<col;i++)
    {
    	for(int j=0;j<col;j++)
        {
 		if(i==j)
 		 {
  			 S[i*col+j]=sqrt(B[i*col+j]);
  		 }
  		else
		{
   			 S[i*col+j]=0;
 		}
 	}
    }
cout<<"________________________________________________________________________"<<endl;
cout<<"sigma matrix...."<<endl;
printmatrix(S,col,col);
cout<<"________________________________________________________________________"<<endl;
multiply(N,G,S,col,col,col,col);
//printmatrix(N,col,col);
multiply(P,N,G2,col,col,col,col);
cout<<"________________________________________________________________________"<<endl;
cout<<"product of U sigma V^T........"<<endl;
//printmatrix(P,col,col);
cout<<"________________________________________________________________________"<<endl;
//printmatrix(u,col,col);



return 0;
}