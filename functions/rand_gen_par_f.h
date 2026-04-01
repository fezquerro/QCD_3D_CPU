#define NormRANu (2.3283063671E-10F)
#if rand_gen==0



struct Rand_data {
unsigned int irr_p[256];
unsigned int ir1_p;
unsigned char ind_ran_p;
unsigned char ig1_p;
unsigned char ig2_p;
unsigned char ig3_p;

};



void initialize_Random(int seed, struct Rand_data* R){
int INI,FACTOR,SUM,i;//initialize the rng for each thread

    INI=seed;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        (*R).irr_p[i]=INI;
    }
    (*R).ind_ran_p=(*R).ig1_p=(*R).ig2_p=(*R).ig3_p=0;

}

double Random_p(struct Rand_data* R)
{
    double r;
   (*R).ig1_p=(*R).ind_ran_p-24;
   (*R).ig2_p=(*R).ind_ran_p-55;
   (*R).ig3_p=(*R).ind_ran_p-61;
   (*R).irr_p[(*R).ind_ran_p]=(*R).irr_p[(*R).ig1_p]+(*R).irr_p[(*R).ig2_p];
   (*R).ir1_p=(*R).irr_p[(*R).ind_ran_p]^(*R).irr_p[(*R).ig3_p];
   (*R).ind_ran_p++;
   r=(*R).ir1_p*NormRANu;
    return r;
}


#endif



