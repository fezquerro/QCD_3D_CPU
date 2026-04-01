#define NormRANu (2.3283063671E-10F)
unsigned int irr[256],ir1;
unsigned char ind_ran,ig1,ig2,ig3;

void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;;

    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;

    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}


double Random(void)
{
    double r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}
