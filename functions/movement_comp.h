void movement_computing(int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,struct lattice l){
 int i;//compute the movements for the derivatives along the lattice
 for(i=0;i<(l.z-1);i++)//the positive derivative in z
  zp[i]=l.dt;//jump for the boundary
 zp[l.z-1]=-(l.z-1)*l.dt;
 for(i=1;i<l.z;i++)//the negative derivative in z
  zm[i]=-l.dt;//jump for the boundary
 zm[0]=(l.z-1)*l.dt;
 for(i=0;i<(l.y-1);i++)//the positive derivative in y
  yp[i]=l.z*l.dt;//jump for the boundary
 yp[l.y-1]=-l.z*(l.y-1)*l.dt;
 for(i=1;i<l.y;i++)//the negative derivative in y
  ym[i]=-l.z*l.dt;//jump for the boundary
 ym[0]=l.z*(l.y-1)*l.dt;
 for(i=0;i<(l.x-1);i++)//the positive derivative in x
  xp[i]=l.z*l.y*l.dt;//jump for the boundary
 xp[l.x-1]=-l.z*l.y*(l.x-1)*l.dt;
 for(i=1;i<l.x;i++)//the negative derivative in x
  xm[i]=-l.z*l.y*l.dt;//jump for the boundary
xm[0]=l.z*l.y*(l.x-1)*l.dt;
for(i=0;i<(l.t-1);i++)//the positive derivative in t
 tp[i]=l.z*l.y*l.x*l.dt;//jump for the boundary
tp[l.t-1]=-l.z*l.y*l.x*(l.t-1)*l.dt;
for(i=1;i<l.t;i++)//the negative derivative in t
 tm[i]=-l.z*l.y*l.x*l.dt;//jump for the boundary
tm[0]=l.z*l.y*l.x*(l.t-1)*l.dt;
}

static inline void compute_inter_z(double* A,double* phi,int pos,int pos_z, int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,int dg){//compute the term that multiplies the U plaquete
double aux[4],aux2[4];
su2_zero(A);
su2_product_tc(phi+pos+zp[pos_z]+dg,phi+pos+yp[pos_y],aux);//y positive plaquette
su2_product_tc(aux,phi+pos+dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+zp[pos_z]+dg+ym[pos_y],phi+pos+ym[pos_y],aux);//y negative plaquette
su2_product(aux,phi+pos+ym[pos_y]+dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);
su2_product_tc(phi+pos+zp[pos_z]+2*dg,phi+pos+xp[pos_x],aux);//x positive plaquette
su2_product_tc(aux,phi+pos+2*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+zp[pos_z]+2*dg+xm[pos_x],phi+pos+xm[pos_x],aux);//x negative plaquette
su2_product(aux,phi+pos+xm[pos_x]+2*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc(phi+pos+zp[pos_z]+3*dg,phi+pos+tp[pos_t],aux);//t positive plaquette
su2_product_tc(aux,phi+pos+3*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+zp[pos_z]+3*dg+tm[pos_t],phi+pos+tm[pos_t],aux);//t negative plaquette
su2_product(aux,phi+pos+tm[pos_t]+3*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

}

static inline void compute_inter_y(double* A,double* phi,int pos,int pos_z, int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,int dg){//compute the term that multiplies the U plaquete
double aux[4],aux2[4];
su2_zero(A);
su2_product_tc(phi+pos+yp[pos_y],phi+pos+zp[pos_z]+dg,aux);//z positive plaquette
su2_product_tc(aux,phi+pos,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+yp[pos_y]+zm[pos_z],phi+pos+zm[pos_z]+dg,aux);//z negative plaquette
su2_product(aux,phi+pos+zm[pos_z],aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);
su2_product_tc(phi+pos+yp[pos_y]+2*dg,phi+pos+xp[pos_x]+dg,aux);//x positive plaquette
su2_product_tc(aux,phi+pos+2*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+yp[pos_y]+2*dg+xm[pos_x],phi+pos+xm[pos_x]+dg,aux);//x negative plaquette
su2_product(aux,phi+pos+xm[pos_x]+2*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc(phi+pos+yp[pos_y]+3*dg,phi+pos+tp[pos_t]+dg,aux);//t positive plaquette
su2_product_tc(aux,phi+pos+3*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

su2_product_tc2(phi+pos+yp[pos_y]+3*dg+tm[pos_t],phi+pos+tm[pos_t]+dg,aux);//t negative plaquette
su2_product(aux,phi+pos+tm[pos_t]+3*dg,aux2);
//printf("%d %d %d ", pos_y, pos_x, pos_t); su2_print_vec(aux2); printf("\n");
su2_sum(A,aux2);

}

static inline void compute_inter_x(double*A,double*phi,int pos,int pos_z,int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,int dg){//compute the term that multiplies the U plaquete
double aux[4],aux2[4];
su2_zero(A);
su2_product_tc(phi+pos+xp[pos_x],phi+pos+zp[pos_z]+2*dg,aux);//z positive plaquette
su2_product_tc(aux,phi+pos,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+xp[pos_x]+zm[pos_z],phi+pos+zm[pos_z]+2*dg,aux);//z negative plaquette
su2_product(aux,phi+pos+zm[pos_z],aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);
su2_product_tc(phi+pos+xp[pos_x]+dg,phi+pos+yp[pos_y]+2*dg,aux);//y positive plaquette
su2_product_tc(aux,phi+pos+dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+xp[pos_x]+ym[pos_y]+dg,phi+pos+ym[pos_y]+2*dg,aux);//y negative plaquette
su2_product(aux,phi+pos+ym[pos_y]+dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc(phi+pos+xp[pos_x]+3*dg,phi+pos+tp[pos_t]+2*dg,aux);//t positive plaquette
su2_product_tc(aux,phi+pos+3*dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+xp[pos_x]+3*dg+tm[pos_t],phi+pos+tm[pos_t]+2*dg,aux);//t negative plaquette
su2_product(aux,phi+pos+tm[pos_t]+3*dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);
}

static inline void compute_inter_t(double*A,double*phi,int pos,int pos_z, int pos_y,int pos_x,int pos_t,int *tp,int *tm,int *xp,int *xm,int *yp,int *ym,int *zp,int *zm,int dg){//compute the term that multiplies the U plaquete
double aux[4],aux2[4];
su2_zero(A);
su2_product_tc(phi+pos+tp[pos_t],phi+pos+zp[pos_z]+3*dg,aux);//z positive plaquette
su2_product_tc(aux,phi+pos,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+tp[pos_t]+zm[pos_z],phi+pos+zm[pos_z]+3*dg,aux);//z negative plaquette
su2_product(aux,phi+pos+zm[pos_z],aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);
su2_product_tc(phi+pos+tp[pos_t]+dg,phi+pos+yp[pos_y]+3*dg,aux);//y positive plaquette
su2_product_tc(aux,phi+pos+dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+tp[pos_t]+ym[pos_y]+dg,phi+pos+ym[pos_y]+3*dg,aux);//y negative plaquette
su2_product(aux,phi+pos+ym[pos_y]+dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc(phi+pos+tp[pos_t]+2*dg,phi+pos+xp[pos_x]+3*dg,aux);//x positive plaquette
su2_product_tc(aux,phi+pos+2*dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);

su2_product_tc2(phi+pos+tp[pos_t]+2*dg+xm[pos_x],phi+pos+xm[pos_x]+3*dg,aux);//x negative plaquette
su2_product(aux,phi+pos+xm[pos_x]+2*dg,aux2);
//su2_print_vec(aux2);
su2_sum(A,aux2);
}
