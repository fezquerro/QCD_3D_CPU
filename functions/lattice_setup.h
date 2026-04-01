
struct lattice {
    int z;
	int y;
	int x;
	int t;
	int V;
	int d;
	int dg;
	int dt;
	float beta;
	float epsilon;
	int updates;

};

double* rand_init (struct lattice l){//random su2 matrices
    double*  phi = calloc(l.d*l.dg*l.z*l.y*l.x*l.t,sizeof(double));/*main field array*/
    for(int i=0; i<l.V; i++)
          {
            su2_rand(phi+i*l.d*l.dg);//z direction U
            su2_rand(phi+i*l.d*l.dg+l.dg);//y direction U
            su2_rand(phi+i*l.d*l.dg+2*l.dg);//x direction U
            su2_rand(phi+i*l.d*l.dg+3*l.dg);//t direction U

            //printf("%lf %lf\n",phi[l],2*Random() - 1);
          }
    return phi;
}
double* identity_init (struct lattice l){// all indentity
    double*  phi = calloc(l.d*l.dg*l.z*l.y*l.x*l.t,sizeof(double));/*main field array*/
    for(int i=0; i<l.V; i++)
          {
            su2_id(phi+i*l.d*l.dg);//z direction U
            su2_id(phi+i*l.d*l.dg+l.dg);//y direction U
            su2_id(phi+i*l.d*l.dg+2*l.dg);//x direction U
            su2_id(phi+i*l.d*l.dg+3*l.dg);//t direction U
          }
    return phi;
}

double* rand_init_pcc(struct lattice l) {//dirichlet initial conditions
    double*  phi = calloc(l.d*l.dg*l.z*l.y*l.x*l.t,sizeof(double));
    for(int t=0;t<l.t;t++){
        for(int x=0;x<l.x;x++){
            for(int y=0;y<l.y;y++){
                for(int z=0;z<l.z;z++){
                if(z==0){
                su2_rand(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z));//z direction U
                su2_id(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+l.dg);//y direction U
                su2_id(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+2*l.dg);//x direction U
                su2_id(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+3*l.dg);//t direction U
                }
                else{
                su2_rand(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z));//z direction U
                su2_rand(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+l.dg);//y direction U
                su2_rand(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+2*l.dg);//x direction U
                su2_rand(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+3*l.dg);//t direction U
                }
                }
            }
        }
    }
return phi;
}

void read_conf(double* phi, struct lattice l,int* Nterm, int* initt, char bound[12], char nsave[2])
{
	char filename[256]; FILE* f;
	snprintf(filename, sizeof(filename), "Config_%d_%d_%d_%d_B_%.2lf_%s_%s.txt", l.t, l.x, l.y, l.z, l.beta, bound, nsave);
	if ((f = fopen(filename, "rb")) == NULL) {
		printf("Error creating the file");
		exit(1);
	}
	fread(phi, sizeof(phi[0]), l.V * l.dt, f);
	fclose(f);
	snprintf(filename, sizeof(filename), "Config_%d_%d_%d_%d_B_%.2lf_%s_%s_data.txt", l.t, l.x, l.y, l.z, l.beta, bound, nsave);
	if ((f = fopen(filename, "r")) == NULL) {
		printf("Error creating the file");
		exit(1);
	}
	fscanf(f, "N_termalization %d Type_termalization %d", Nterm, initt);
	fclose(f);
}


double* check_init (struct lattice l){// all indentity
    double*  phi = calloc(l.d*l.dg*l.y*l.x*l.t,sizeof(double));/*main field array*/
    for(int t=0;t<l.t;t++){
        for(int x=0;x<l.x;x++){
            for(int y=0;y<l.y;y++){
                su2_count(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x),y+1.33*x+1.55*t+0.111);//y direction U
                su2_count(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+l.dg,y+1.33*x+1.55*t+0.111*2);//x direction U
                su2_count(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+2*l.dg,y+1.33*x+1.55*t+0.111*3);//t direction U
                }

            }
        }
    return phi;
}
void printlattice_matrix(double* phi, struct lattice l)//prints the matrices of the lattice
{
	for (int t = 0; t < l.t; t++)
	{
		for (int x = 0; x < l.x; x++)
		{
			for (int y = 0; y < l.y; y++){
			su2_print(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x));
			printf("\n");
			su2_print(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+l.dg);
			printf("\n");
			su2_print(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+2*l.dg);
			printf("\n");
			}
			printf("\n\n");

		}
		printf("\n\n\n\n");
	}

}
void printlattice(double* phi, struct lattice l)//prints the 4 real values of each U of the lattice
{
	for (int t = 0; t < l.t; t++)
	{
		for (int x = 0; x < l.x; x++)
		{
			for (int y = 0; y < l.y; y++)
			{
              for (int z = 0; z < l.z; z++)
              {
                  su2_print_vec(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z));
			      su2_print_vec(phi+l.d*l.dg*(z+y*l.z+x*l.y*l.z+t*l.y*l.x*l.z)+l.dg);
			      su2_print_vec(phi+l.d*l.dg*(z+l.z*y+x*l.y*l.z+t*l.y*l.x*l.z)+2*l.dg);
			      su2_print_vec(phi+l.d*l.dg*(z+l.z*y+x*l.y*l.z+t*l.y*l.x*l.z)+3*l.dg);
			      printf("\n");


              }
              printf("\n\n");
			}
			printf("\n\n\n");

		}
		printf("\n\n\n\n");
	}

}
void check_unit(double* phi, struct lattice l){//prints the value of the determinant of each U
for (int t = 0; t < l.t; t++)
	{
		for (int x = 0; x < l.x; x++)
		{
			for (int y = 0; y < l.y; y++){
            printf("%.16lf %.16lf %.16lf",su2_det(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)),su2_det(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+l.dg),su2_det(phi+l.d*l.dg*(y+x*l.y+t*l.y*l.x)+2*l.dg));
			printf("\n");
			}
			printf("\n\n");

		}
		printf("\n\n\n\n");
	}


}

void print_det(double* phi, struct lattice l){
    int i;double det=0;
    for(i=0;i<l.V*l.d;i++){
            //printf("%.16lf\n",su2_det(phi+i*l.dg));
            det+=su2_det(phi+i*l.dg);
		}
    printf("Av det %.16lf\n",det/l.V/l.d);
}
void normalize(double* phi, struct lattice l){
 for(int i=0;i<l.V*l.d;i++)
        su2_normalize(phi+i*l.dg);

}

void normalize_pf(double* phi, struct lattice l,int id,int tthread){

 for(int i=0;i<tthread*l.x*l.y*l.z*l.d;i++)
        su2_normalize(phi+i*l.dg+id*tthread*l.x*l.y*l.z*l.dt);

#pragma omp barrier
}
