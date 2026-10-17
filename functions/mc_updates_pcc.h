


static inline int sweep_one_point_pcc(double*phi,int pos_z,int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int*zp,int*zm, struct lattice l,uint64_t* R){
int pos=l.d*l.dg*(pos_z+pos_y*l.z+l.z*l.y*pos_x+l.z*l.y*l.x*pos_t),aux=0,i;
double A[4];
//printf("\n%d\n",pos/l.d/l.dg);
compute_inter_z(A,phi,pos,pos_z,pos_y,pos_x,pos_t,tp,tm,xp,xm,yp,ym,zp,zm,l.dg);
//printf("%d %d %d %d %lf %lf %lf %lf\n",pos_z,pos_y,pos_x,pos_t,A[0],A[1],A[2],A[3]);
//su2_print_vec(A);
aux+=heat_bath(A,phi+pos,l.beta,R);
return aux;
}

static inline double MC_sweep_pcc(double* phi,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int*zm,struct lattice l,uint64_t* R,int id,int tthread){//update the full lattice
int accept2=0,t,x,y,z,taux,t_ini;
    for(taux=0;taux<tthread/2;taux++){
        t=2*taux+id*tthread;
        for(x=0; x<l.x; x++)
          {
              for (y=0; y<l.y;y++)
              {

                  for (z=0;z<l.z;z++){
                    if (!(z==0))
                   accept2+=sweep_one_point(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);
                   else
                   accept2+=sweep_one_point_pcc(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);

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

                  for (z=0;z<l.z;z++){
                    if (!(z==0))
                   accept2+=sweep_one_point(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);
                   else
                   accept2+=sweep_one_point_pcc(phi,z,y,x,t,tp,tm,xp,xm,yp,ym,zp,zm,l,R);

              }
          }
      }
      }

      #pragma omp barrier
      return (double)accept2/l.V/l.d/l.updates;
}


