/* double drand64(void);  */
#undef MERSENNE_TWISTER                   /* random number, Mersenne or lcg */
#undef NR_RANDOM_NUMBER              /* #define use "Numerical Recipes" ran2 */
#define LCG64                              /* use 64-bit linear congurential */

#include <assert.h>
#include <stdio.h>
                                                                  /* globals */
#ifdef NR_RANDOM_NUMBER   /* then use "Numerical Recipes" random number ran2 */
                                     /* from "NR in C" on page 282, modified */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-13
#define RNMX (1.0-EPS)

static int idum = 1;

double drand64(void)               /* this is "Numerical Recipes" ran2(idum) */
{ 
    int j;
    int k;
    static int idum2=123456789;
    static int iy=0;
    static int iv[NTAB];
    double temp;
    
    if( idum <= 0 )
    {
        if( ((-1)*idum ) < 1 )
            idum = 1;
        else
            idum  = (-1)*idum;
        idum2 = idum;
        for( j = NTAB+7 ; j>=0 ; j-- )
        {
            k = idum/IQ1;
            idum = IA1 * (idum-k*IQ1)- k*IR1;
            if( idum < 0 )
                idum += IM1;
            if( j< NTAB )
                iv[ j ] = idum;
        }
        iy = iv[0];
    }
    
    k = idum/IQ1;
    idum = IA1 * (idum-k*IQ1)- k*IR1;
    if( idum < 0 )
        idum += IM1;
    k = idum2/IQ2;
    idum2 = IA2 * (idum2-k*IQ2)- k*IR2;
    if( idum2 < 0 )
        idum2 += IM2;
    j = iy/NDIV;
    iy = iv[ j ]-idum2;
    iv[ j ] = idum;
    if( iy < 1 )
        iy += IMM1;
    if( (temp = AM*iy ) > RNMX )
        return RNMX;
    else
        return temp; 
}  

void srand64(int seed, FILE *fp)                           /* initialization */
{
   assert(sizeof(int) == 4);                    /* work only with 32-bit int */
   idum = -seed;
   if(idum > 0) idum = - idum;
   drand64();                       /* call with negative idum to initialize */
   fprintf(fp, "NR ran2, initial seed x = %d\n", seed);
}

 
#elif defined(LCG64)                       /* 64-bit random number generator */
                 /* from Knuth "Art of Computer Programming Vol 2", page 108 */
static unsigned long long int x = 1;

double drand64(void) 
{
   x = 6364136223846793005ll * x + (long long int) 1;
   return (double) x * 5.4210108624275218e-20; 
}

void srand64(int seed, FILE *fp)
{
   assert(sizeof(long long int) == 8);
   x = seed;
   fprintf(fp, "drand 64-bit, initial seed x = %d\n", seed);
}
#else 

            /* from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html */
                                             /* with interface modfification */
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed). 

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Any feedback is very welcome.
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned int mt[N];  /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void srand64(int s, FILE *fp)
{
    assert(sizeof(int)==4);

    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
   fprintf(fp, "Mersenne Twister, initial seed x = %d\n", (int) s);
}


/* generates a random number on [0,0xffffffff]-interval */
unsigned int genrand_int32(void)
{
    unsigned int y;
    static unsigned int mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

/****
    assert(mti != N+1);
        if (mti == N+1)  
            init_genrand(5489UL); 
****/

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}


/* generates a random number on [0,1) with 53-bit resolution*/
double drand64(void) 
{ 
    unsigned int a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

#endif
