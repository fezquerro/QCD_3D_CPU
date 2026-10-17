
#define pi 3.14159265358979324

static inline int heat_bath(double*A,double* phi, double beta,uint64_t* R){//Kennedy-Pendleton without trigonometric calls
double a=sqrt(su2_det(A)),l,u,v,s,c2,aux,n,X[4],B[4];
int c=0;
while (c==0){
do{//cos^2 of a uniform angle as u^2/(u^2+v^2) with (u,v) uniform in the unit disk
    u=2*Random_xor(R)-1;
    v=2*Random_xor(R)-1;
    s=u*u+v*v;
}while(s>1||s==0);
c2=u*u/s;
l=-(log(1-Random_xor(R))+c2*log(1-Random_xor(R)))/(2*a*beta);
aux=Random_xor(R);
if(aux*aux<1-l)
    c=1;
}
X[0]=1-2*l;n=sqrt(1-X[0]*X[0]);
do{//Marsaglia: uniform direction on the sphere
    u=2*Random_xor(R)-1;
    v=2*Random_xor(R)-1;
    s=u*u+v*v;
}while(s>=1);
aux=2*n*sqrt(1-s);
X[1]=n*(1-2*s);
X[2]=aux*u;
X[3]=aux*v;
aux=1/a;
for(int i=0;i<4;i++)
    B[i]=A[i]*aux;
su2_product_tc(X,B,phi);
return 1;
}

//Previous heat bath (Kennedy-Pendleton with trigonometric calls), kept for reference.
//To recover it, uncomment it and comment out the heat_bath above.
//int heat_bath(double*A,double* phi, double beta,uint64_t* R){
//double a=sqrt(su2_det(A)),r[3],l,aux,n,X[4],B[4],C[4],psi,ctheta;
//int c=0,i;
//while (c==0){
//for(i=0;i<3;i++)
//    r[i]=1-Random_xor(R);
//aux=cos(2*pi*r[1]);
//l=-(log(r[0])+aux*aux*log(r[2]))/(2*a*beta);
//aux=Random_xor(R);
//if(aux*aux<1-l)
//    c=1;
//}
//X[0]=1-2*l;n=sqrt(1-X[0]*X[0]);
//
///*c=0;
//while (c==0){
//for(i=0;i<3;i++)
//    r[i]=2*Random()-1;
//aux=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
//if(aux<=1)
//    c=1;
//
//}
//for(i=0;i<3;i++){
//    X[i+1]=r[i]*n/sqrt(aux);
//}
//*/
//psi=2*pi*Random_xor(R);
//ctheta=2*Random_xor(R)-1;
//aux = n * sqrt(1 - ctheta * ctheta);//n*arcsin(theta)
//X[1]=n*ctheta;
//X[2]=aux*cos(psi);
//X[3]=aux*sin(psi);
//
//su2_copy(A,B);
//su2_normalize(B);
//su2_product_tc(X,B,phi);
//
//return 1;
//}
static inline int sweep_one_point(double*phi,int pos_z,int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm, struct lattice l,uint64_t* R){
int pos=l.d*l.dg*(pos_z+pos_y*l.z+l.z*l.y*pos_x+l.z*l.y*l.x*pos_t),aux=0,i;
double A[4];
//printf("\n%d\n",pos/l.d/l.dg);
compute_inter_z(A,phi,pos,pos_z,pos_y,pos_x,pos_t,tp,tm,xp,xm,yp,ym,zp,zm,l.dg);
//su2_print_vec(A);
//printf("%d %d %d %d %lf %lf %lf %lf %d\n",pos_z,pos_y,pos_x,pos_t,A[0],A[1],A[2],A[3],j);
aux+=heat_bath(A,phi+pos,l.beta,R);
compute_inter_y(A,phi,pos,pos_z,pos_y,pos_x,pos_t,tp,tm,xp,xm,yp,ym,zp,zm,l.dg);
//printf("%d %d %d %d %lf %lf %lf %lf\n",pos_z,pos_y,pos_x,pos_t,A[0],A[1],A[2],A[3]);
//su2_print_vec(A);
aux+=heat_bath(A,phi+pos+l.dg,l.beta,R);
compute_inter_x(A,phi,pos,pos_z,pos_y,pos_x,pos_t,tp,tm,xp,xm,yp,ym,zp,zm,l.dg);
//printf("%d %d %d %d %lf %lf %lf %lf\n",pos_z,pos_y,pos_x,pos_t,A[0],A[1],A[2],A[3]);
//su2_print_vec(A);
aux+=heat_bath(A,phi+pos+2*l.dg,l.beta,R);
compute_inter_t(A,phi,pos,pos_z,pos_y,pos_x,pos_t,tp,tm,xp,xm,yp,ym,zp,zm,l.dg);
//printf("%d %d %d %d %lf %lf %lf %lf\n",pos_z,pos_y,pos_x,pos_t,A[0],A[1],A[2],A[3]);
//su2_print_vec(A);
aux+=heat_bath(A,phi+pos+3*l.dg,l.beta,R);
return aux;
}

static inline double MC_sweep(double* phi,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,struct lattice l,uint64_t* R,int id,int tthread){//update the full lattice
int accept2=0,t,x,y,z,taux,t_ini;


//t_ini = clock();





       for(taux=0;taux<tthread/2;taux++){
        t=2*taux+id*tthread;
        for(x=0; x<l.x; x++)
          {
              for (y=0; y<l.y;y++)
              {
                   for (z=0; z<l.z;z++){
                   accept2+=sweep_one_point(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);

                   }
              }
          }
        }

          #pragma omp barrier




          for(taux=0;taux<tthread/2;taux++){
              t=2*taux+id*tthread+1;
          for(x=0; x<l.x; x++)
          {
              for (y=0; y<l.y;y++)
              {
                   for (z=0; z<l.z;z++)
                   accept2+=sweep_one_point(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);
              }
          }
          }


          #pragma omp barrier

//printf("%d %d\n",id,clock()-t_ini);
//printf("%lf\n",(double) accept/l.V/l.d/l.updates);
      return (double)accept2/l.V/l.d/l.updates;
}


