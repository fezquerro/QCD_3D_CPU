#include <stdlib.h>//Standard C Library
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>
#include <time.h>


#include "functions/rand_gen_par_f.h"
#include "functions/rand_gen_par_xor.h"
#include "functions/rand_gen.h"
#include "functions/su2.h"
#include "functions/lattice_setup.h"
#include "functions/movement_comp.h"
#include "functions/mc_updates.h"
#include "functions/mc_updates_pcc.h"
#include "functions/observables.h"
#include "functions/statistics.h"
#include "functions/output.h"


#define nthreads 4 //number of threads for parallel computing, MUST BE A DIVISOR OF d_t
#define d_z 8 //z direction size
#define d_y 8 //y direction size
#define d_x 8 //x direction size
#define d_t 8 //t direction size

#define Boundary 0 //define boundary condition, 0 for Periodic and 1 for Perfect colour conductor
#define readconf 0 //read initial configuration if 1 and thermalize, read conf and not thermalize if 2
#define init 0 //0 hot start, 1 cold start
#define sconf 1 // save configuration after thermalization if 0
#define save 0 //save data if 0

#define betal 2.6 // coupling value
#define epsilonl 0.3 //range of update
#define upd 1

#define N 1000 //number of measures
#define N_term 200 //number of thermalization
#define N_block 5 // number of times it prints the result while running
#define Nsave 5 //number of times it saves the results while running
#define Nblock_term 5 // number of times it prints the result while thermalizing
#define N_cor 3 //MC steps between measures

#define W 80
int main(){
char nsave[3]="1";
char nconfr[2] = "1";
char nconfs[2] = "1";
int Nterms=0,inits;

struct lattice l = { d_z,d_y,d_x,d_t, d_z*d_y * d_x * d_t,4,4,16,betal,epsilonl,upd };// dimz,dimy,dimx,dimt,V,dim, dim group rep,beta,delta
printf("N %d Nterm %d Lt %d Lx %d Ly %d Lz %d Beta %.6lf Noutput %s nthreads %d\n",N,N_term,l.t,l.x,l.y,l.z,l.beta,nsave,nthreads);
unsigned int try1[16],try2[16];
int thread_id;
//ini_ran(123);
ini_ran(time(NULL));


omp_set_num_threads(nthreads);//set number of threads
printf("N proc %d\n",nthreads);

#if Boundary==0
printf("PERIODIC BOUNDARY CONDITIONS\n\n");
char bound[12] = "per";
 #if readconf==0
 inits = init;
   #if init==0
      double*  phi=rand_init(l);//hot start
   #else
      double*  phi=identity_init(l);//cold start
   #endif
 #else
  double* phi = (double*)malloc(l.V * l.dt * sizeof(double));
  read_conf(phi, l,&Nterms,&inits, bound, nconfr);//read an initial configurations
 #endif
#else
 printf("PERFECT COLOUR CONDUCTOR BOUNDARY CONDITIONS\n\n");
 char bound[12] = "pcc";
 #if readconf==0
 inits = init;
  #if init==0
   double*  phi=rand_init_pcc(l);//hot start
  #else
   double*  phi=identity_init(l);//cold start
  #endif
 #else
  double* phi = (double*)malloc(l.V * l.dt * sizeof(double));
  read_conf(phi, l,&Nterms,&inits, bound, nconfr);//read an initial configurations
 #endif
#endif
fflush(stdout);

//double* phi=check_init(l);
//printlattice(phi,l); //prints the 4 values representing the su2
//check_unit(phi,l);//prints the determinant
print_det(phi,l);


int tp[l.t],tm[l.t],xp[l.x],xm[l.x],yp[l.y],ym[l.y],zp[l.z],zm[l.z],i,j; movement_computing(tp,tm,xp,xm,yp,ym,zp,zm,l);//Compute the movements along the lattice for the derivatives
double* measure = (double*)malloc(N * sizeof(double)), accept=0;//define the array where the measures will be saved
double* measure2 = (double*)malloc(N * sizeof(double));


for(i=0;i<N;i++){
    measure[i]=0;
    measure2[i]=0;
}
int t_ini, t_fin, t_n,block_size_t = N_term / Nblock_term, block_size = N / N_block,block_s_size=N/Nsave;//blocks for print time and saving
t_ini = time(NULL); t_n = time(NULL);
double energy[2],energy2[2];
double autoc[W], autoc2[W];
uint64_t seed=time(NULL);

#pragma omp parallel private(i,j) reduction(+:accept)
{

int id=omp_get_thread_num();
int tthread=l.t/nthreads;// number of time sites each thread doess



uint64_t sr[4];
//uint32_t sr[4];
sr[0]=Rand_init(seed);
//sr[0]=seed;
for(i=1;i<4;i++)
sr[i]=Rand_init(sr[i-1]);
//sr[i]=seed+i*9821764;
for(i=0;i<id;i++)
long_jump(sr);
printf("Hilo %d\n",id);
fflush(stdout);
int t_inis = clock();

#if readconf<2
for (i=1;i<N_term+1;i++){//loop for thermalization
     #if Boundary==0
     accept+=MC_sweep(phi,tp,tm,xp,xm,yp,ym,zp,zm,l,sr,id,tthread);
      #else
     accept+=MC_sweep_pcc(phi,tp,tm,xp,xm,yp,ym,zp,zm,l,sr,id,tthread);
     #endif

    // print_det(phi,l);
   if(i%100==0)
       normalize_pf(phi,l,id,tthread);//project over su2 by normalization

    #pragma omp master
    {


    if (i % block_size_t == 0)  {
				printf("%lf Time %d Left(min)%lf\n", (double)i / (N_term+1), time(NULL) - t_n, (double)(time(NULL) - t_ini) * (N_term+1 - i) / i / 60);
				fflush(stdout);
				t_n = time(NULL);
			}
    }
}


#pragma omp master
{

printf("Accept term %lf\n",accept/(N_term+1));
Nterms = N_term + Nterms;
#if sconf == 0//save configuration
			normalize(phi,l);
			save_conf(phi, l,Nterms,inits, bound, nconfs);
#endif
print_det(phi,l);
}

accept=0;
#endif




		for (i = 1; i < N + 1; i++) {//loop for measuring
			for (j = 0; j < N_cor; j++) {
				 #if Boundary==0
                 accept+=MC_sweep(phi,tp,tm,xp,xm,yp,ym,zp,zm,l,sr,id,tthread);
                 #else
                 accept+=MC_sweep_pcc(phi,tp,tm,xp,xm,yp,ym,zp,zm,l,sr,id,tthread);
                 #endif
			}
			energy_av_par(phi,tp,xp,yp,zp,l,measure+i-1,measure2+i-1,id,tthread);// measure in parallel the energy and plaquette average

			// print_det(phi,l);
            if(i%100==0)
               normalize_pf(phi,l,id,tthread);//project over su2 by normalization

            #pragma omp master
            {


			if (i % block_size == 0) {//print time left
				printf("%lf Time %d Left(hour)%lf\n", (double)i / N, time(NULL) - t_n, (double)(time(NULL) - t_ini) * (N - i) / i / 60/60);
				fflush(stdout);
				t_n = time(NULL);
			}

			if (i % block_s_size == 0){//compute average and save them
                mean_variance(measure, i, energy, true);
		        mean_variance(measure2, i, energy2, true);
                autocorrelation(measure, i, W, energy[0], energy[1], autoc); autocorrelation(measure2, i, W, energy2[0], energy2[1], autoc2);
                fflush(stdout);
            #if save==0
                save_energy(energy[0], energy[1], Nterms, i, N_cor, accept, l,bound,nsave);
		        save_energy_av(energy[0], energy[1], energy2[0], energy2[1], Nterms, i, N_cor, accept, l,bound,nsave);
                save_autocorrelation(autoc, W, i, N_cor, Nterms, l,bound,nsave); save_autocorrelation_av(autoc2, W, i, N_cor, Nterms, l,bound,nsave);
                save_en_seq(measure2,measure,Nterms,i,N_cor,1,l,bound,nsave);
            #endif
			}
            }
		}
		//printf("End thread %d %d %lf\n",i,id,accept);
}
		t_fin = time(NULL);
		printf("%d\n", t_fin - t_ini);
        print_det(phi,l);
        #if Boundary==1
		accept = accept / ((double)(l.z-1)/l.z+1/(double)l.z/4);
        #endif
		printf("Accept %lf\n",accept/N/N_cor);


		mean_variance(measure, N, energy, true);
		mean_variance(measure2, N, energy2, true);
        autocorrelation(measure, N, W, energy[0], energy[1], autoc); autocorrelation(measure2, N, W, energy2[0], energy2[1], autoc2);

        #if save==0

        save_energy(energy[0], energy[1], Nterms, N, N_cor, accept, l,bound,nsave);
		save_energy_av(energy[0], energy[1], energy2[0], energy2[1], Nterms, N, N_cor, accept, l,bound,nsave);

		save_autocorrelation(autoc, W, N, N_cor, Nterms, l,bound,nsave); save_autocorrelation_av(autoc2, W, N, N_cor, Nterms, l,bound,nsave);
        save_en_seq(measure2,measure,Nterms,N,N_cor,1,l,bound,nsave);
        #endif

free(phi);free(measure);free(measure2);
}
