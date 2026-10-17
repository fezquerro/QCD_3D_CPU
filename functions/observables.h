
static inline double plaquette_zy(double* phi,int pos_z, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int *zp,int dg){
double A[4],B[4];
su2_product(phi,phi+zp[pos_z]+dg,A);
su2_product_tc(A,phi+yp[pos_y],B);
su2_product_tc(B,phi+dg,A);
return 1-A[0];
}
static inline double plaquette_zx(double* phi,int pos_z, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int *zp,int dg){
double A[4],B[4];
su2_product(phi,phi+zp[pos_z]+2*dg,A);
su2_product_tc(A,phi+xp[pos_x],B);
su2_product_tc(B,phi+2*dg,A);
return 1-A[0];
}
static inline double plaquette_zt(double* phi,int pos_z, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int *zp,int dg){
double A[4],B[4];
su2_product(phi,phi+zp[pos_z]+3*dg,A);
su2_product_tc(A,phi+tp[pos_t],B);
su2_product_tc(B,phi+3*dg,A);
return 1-A[0];
}

static inline double plaquette_yx(double* phi, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int dg){
double A[4],B[4];
su2_product(phi+dg,phi+yp[pos_y]+2*dg,A);
su2_product_tc(A,phi+xp[pos_x]+dg,B);
su2_product_tc(B,phi+2*dg,A);
return 1-A[0];
}
static inline double plaquette_yt(double* phi, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int dg){
double A[4],B[4];
su2_product(phi+dg,phi+yp[pos_y]+3*dg,A);
su2_product_tc(A,phi+tp[pos_t]+dg,B);
su2_product_tc(B,phi+3*dg,A);
return 1-A[0];
}
static inline double plaquette_xt(double* phi, int pos_y,int pos_x, int pos_t,int *tp,int *xp,int *yp,int dg){
double A[4],B[4];
su2_product(phi+2*dg,phi+xp[pos_x]+3*dg,A);
su2_product_tc(A,phi+tp[pos_t]+2*dg,B);
su2_product_tc(B,phi+3*dg,A);
return 1-A[0];
}



static inline void energy_av_par(double*phi,int *tp,int *xp,int *yp,int *zp,struct lattice l,double* energ,double* plaq,double* prof,int id, int tthread){
int t,taux,x,y,z,pos,t_ini;double s=0,s1=0,a,b,c,d,e,f,pz[l.z],pnz[l.z];//pz: planes with z (centred at z+1/2), pnz: planes without z (at z)
for(z=0;z<l.z;z++){
    pz[z]=0;
    pnz[z]=0;
}


//t_ini = clock();



    for(taux=0;taux<tthread;taux++){
         t=taux+id*tthread;
        // printf("%d %d\n",t,id);
         for(x=0; x<l.x; x++)
          {
              for (y=0; y<l.y;y++)
              {
                  for (z=0; z<l.z;z++)
                  {
                   pos=l.d*l.dg*(z+y*l.z+l.z*l.y*x+l.z*l.y*l.x*t);
                   a=plaquette_zy(phi+pos,z,y,x,t,tp,xp,yp,zp,l.dg);
                   b=plaquette_zx(phi+pos,z,y,x,t,tp,xp,yp,zp,l.dg);
                   c=plaquette_zt(phi+pos,z,y,x,t,tp,xp,yp,zp,l.dg);
                   d=plaquette_yx(phi+pos,y,x,t,tp,xp,yp,l.dg);
                   e=plaquette_yt(phi+pos,y,x,t,tp,xp,yp,l.dg);
                   f=plaquette_xt(phi+pos,y,x,t,tp,xp,yp,l.dg);
                   s+=(6-a-b-c-d-e-f)/6;
                   s1+=(a + b + c) / 3 - (d + e + f) / 3;
                   pz[z]+=(a + b + c) / 3;
                   pnz[z]+=(d + e + f) / 3;
                   //s1+=a-b-c;
                   //s1+=a+b-c+d-e-f;
                   //printf("%d %d %d %d %lf %lf\n", z, y, x, t, (a + b + c) / 3 - (d + e + f) / 3, (6 - a - b - c - d - e - f) / 6);
                  }

              }
          }
        }



#pragma omp atomic update
*energ+= s1/l.V;
#pragma omp atomic update
   *plaq+=s/l.V;
for(z=0;z<l.z;z++){//site-centred T_zz: average the z planes at z-1/2 and z+1/2
#pragma omp atomic update
   prof[z]+=((pz[(z-1+l.z)%l.z]+pz[z])/2-pnz[z])*l.z/l.V;
}
//printf("%d %d\n",id,clock()-t_ini);




#pragma omp barrier
//printf("measure %d\n",id);
     // printf("%lf %lf\n",s1/l.V,s/l.V);
}
