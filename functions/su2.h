void su2_rand(double*A){//generate random su2 matrices
for(int i=0;i<4;i++)
    A[i]=Random()-0.5;
double n=sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3]);
for(int i=0;i<4;i++)
    A[i]=A[i]/n;
}
void su2_id(double*A){//generate the identity
A[0]=1;
for(int i=1;i<4;i++)
    A[i]=0;
}

double su2_det(double*A){// determinant of as SU2 matrix given in sigma basis
return A[0]*A[0]+A[1]*A[1]+A[2]*A[2]+A[3]*A[3];
}
void su2_normalize(double*A){//normalization of a SU2 matrix given in sigma basis
double d=sqrt(su2_det(A));
for (int i = 0; i < 4; i++)
    A[i]=A[i]/d;
}
void su2_count(double*A,double a){//generate unitary with a parameter for checks

for(int i=0;i<4;i++)
    A[i]=0.05*(a+i);
su2_normalize(A);
}
void su2_zero(double*A){//zero vector
for(int i=0;i<4;i++)
    A[i]=0;
}
void su2_product(double* A,double* B, double* C){ //define the product of 2 SU2 matrices in sigma basis
C[0]=A[0]*B[0]-A[1]*B[1]-A[2]*B[2]-A[3]*B[3];
C[1]=A[1]*B[0]+A[0]*B[1]+A[2]*B[3]-A[3]*B[2];
C[2]=A[2]*B[0]+A[3]*B[1]+A[0]*B[2]-A[1]*B[3];
C[3]=A[3]*B[0]-A[2]*B[1]+A[0]*B[3]+A[1]*B[2];
}
void su2_product_tc(double* A,double* B, double* C){ //define the product of 2 SU2 matrices in sigma basis with second CT
C[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3];
C[1]=A[1]*B[0]-A[0]*B[1]-A[2]*B[3]+A[3]*B[2];
C[2]=A[2]*B[0]-A[3]*B[1]-A[0]*B[2]+A[1]*B[3];
C[3]=A[3]*B[0]+A[2]*B[1]-A[0]*B[3]-A[1]*B[2];
}
void su2_product_tc2(double* A,double* B, double* C){ //define the product of 2 SU2 matrices in sigma basis both CT
C[0]=A[0]*B[0]-A[1]*B[1]-A[2]*B[2]-A[3]*B[3];
C[1]=-A[1]*B[0]-A[0]*B[1]+A[2]*B[3]-A[3]*B[2];
C[2]=-A[2]*B[0]+A[3]*B[1]-A[0]*B[2]-A[1]*B[3];
C[3]=-A[3]*B[0]-A[2]*B[1]-A[0]*B[3]+A[1]*B[2];
}
double su2_prod_trace(double* A,double* B){// return the real part of 1/2 of the trace of the product
return A[0]*B[0]-A[1]*B[1]-A[2]*B[2]-A[3]*B[3];
}
void su2_sum(double*A,double*B){
for(int i=0;i<4;i++)
    A[i]=A[i]+B[i];
}
void su2_sub(double*A,double*B,double*C){
for(int i=0;i<4;i++)
    C[i]=A[i]-B[i];
}

void su2_update(double*A,double *B,double ep){//update by multiplying by a random SU2 matrix, ep close to identity
double r[4];
for(int i=1;i<4;i++)
    r[i]=Random()-0.5;
double n=sqrt(r[1]*r[1]+r[2]*r[2]+r[3]*r[3]);
for(int i=0;i<3;i++)
    r[i+1]=ep*r[i+1]/n;
//if(r[0]>0)
 r[0]=sqrt(1-ep*ep);
//else
 //r[0]=-sqrt(1-ep*ep);

su2_product(r,A,B);
}
void su2_print(double*A){//print a SU2 matrix given in sigma basis
printf("%lf + i %lf , %lf +i %lf\n%lf + i %lf , %lf +i %lf\n",A[0],A[1],A[2],A[3],-A[2],A[3],A[0],-A[1]);
}
void su2_print_vec(double*A){//print th 4 real values representing U
printf("%lf,%lf,%lf,%lf  ",A[0],A[1],A[2],A[3]);
}
void su2_copy(double*A, double *B){
for(int i=0;i<4;i++)
    B[i]=A[i];
}



