#include <cstdio>
#include <iostream>
#include <cmath>
#include "rndmCPP.h"
#include "mathconst.h"

RNDM::RNDM(unsigned long seed)
{
  rndmst0(1);
  rndmst1(1);
  rndmst2(1);
  //  rndmst3(1);
  rndmst5(12,34,56,78);
  rndmst6(1);
  rndmst7(seed);
  rndmst8(1);
  iset_ = 0;
  v1_ = 0.0, v2_ = 0.0;
  counterRndm7_ = 0;
}

void  RNDM::rndmst0(int i)
{
  rndm0_store.i=i;
}

float RNDM::rndm0()
{
  const int m=2147483647,a=16807,q=127773,r=2836;
  const float m_inv=(1.0-RNDM_EPS)/(float)m;
  int k;

  k=rndm0_store.i/q;
  if((rndm0_store.i=a*(rndm0_store.i-k*q)-r*k)<0) rndm0_store.i+=m;
  return m_inv*rndm0_store.i;
}


void RNDM::rndmst1(int i)
{
  //    const long m=2147483647,m1=m-1,a=16807,q=127773,r=2836,n=32,nd=1+m1/n;
  //   const float m_inv=(1.0-RNDM_EPS)/(float)m;
  const long m=2147483647, a=16807,q=127773,r=2836,n=32;
  int k,j;

  for (j=n+7;j>=0;j--)
    {
      k=i/q;
      if((i=a*(i-k*q)-r*k)<0) i+=m;
      if (j<n) rndm1_store.is[j]=i;
    }
  rndm1_store.i=i;
  rndm1_store.p=rndm1_store.is[0];
}

float RNDM::rndm1()
{
  const long m=2147483647,m1=m-1,a=16807,q=127773,r=2836,n=32,nd=1+m1/n;
  const float m_inv=(1.0-RNDM_EPS)/(float)m;
  int k,j;

  k=rndm1_store.i/q;
  if((rndm1_store.i=a*(rndm1_store.i-k*q)-r*k)<0) rndm1_store.i+=m;
  j=rndm1_store.p/nd;
  rndm1_store.p=rndm1_store.is[j];
  rndm1_store.is[j]=rndm1_store.i;
  return rndm1_store.p<m1?m_inv*rndm1_store.p:m_inv*m1;
/*  return (m_inv*min(rndm1_store.p,m1));*/
}


void RNDM::rndmst2(int i)
{
  //   const int m1=2147483563,m2=2147483399,m1_1=m1-1;
  //    const int a1=40014,a2=40692,q1=53668,q2=52774,r1=12211,r2=3791,n=32, nd=1+m1_1/n;
  const int m1=2147483563;
  const int a1=40014,q1=53668,r1=12211,n=32;
/*  const float m1_inv=(1.0-RNDM_EPS)/(float)m1,m2_inv=(1.0-RNDM_EPS)/(float)m2;*/
  int k,j;

  for (j=n+7;j>=0;j--)
    {
      k=i/q1;
      if((i=a1*(i-k*q1)-r1*k)<0) i+=m1;
      if (j<n) rndm2_store.is[j]=i;
    }
  rndm2_store.i1=i;
  rndm2_store.i2=12345678;
  rndm2_store.p=rndm2_store.is[0];
}

float RNDM::rndm2()
{
  const int m1=2147483563,m2=2147483399,m1_1=m1-1;
  const int a1=40014,a2=40692,q1=53668,q2=52774,r1=12211,r2=3791,n=32,
            nd=1+m1_1/n;
  const float m_inv=(1.0-RNDM_EPS)/(float)m1;
  int k,j;

  k=rndm2_store.i1/q1;
  if((rndm2_store.i1=a1*(rndm2_store.i1-k*q1)-r1*k)<0) rndm2_store.i1+=m1;
  k=rndm2_store.i2;
  if((rndm2_store.i2=a2*(rndm2_store.i2-k*q2)-r2*k)<0) rndm2_store.i2+=m2;
  j=rndm2_store.p/nd;
  if((rndm2_store.p=rndm2_store.is[j]-rndm2_store.i2)<1) rndm2_store.p+=m1_1;
  rndm2_store.is[j]=rndm2_store.i2;
  return rndm2_store.p*m_inv;
/*  return rndm2_store.p<m1_1?m_inv*rndm2_store.p:m_inv*m1_1;*/
/*  return m_inv*min(rndm2_store.p,m1_1);*/
}

void RNDM::rndmst5(int na1,int na2,int na3, int nb1)
{
  int i,j,nat;
  float s,t;
  rndm5_store.i=96;
  rndm5_store.j=32;
  for (i=0;i<97;i++)
    {
      s=0.0;
      t=0.5;
      for (j=0;j<24;j++)
	{
	  nat=(((na1*na2) % 179)*na3) % 179;
	  na1=na2;
	  na2=na3;
	  na3=nat;
	  nb1=(53*nb1+1) % 169;
	  if ((nb1*nat) % 64 >= 32)
	    {
	      s+=t;
	    }
	  t*=0.5;
	}
      rndm5_store.u[i]=s;
    }
  rndm5_store.c=    362436.0/16777216.0;
  rndm5_store.cd=  7654321.0/16777216.0;
  rndm5_store.cm= 16777213.0/16777216.0;
}

float RNDM::rndm5()
{
  float temp;

/* for (;;){*/
  temp=rndm5_store.u[rndm5_store.i]-rndm5_store.u[rndm5_store.j];
  if (temp<0.0)
    {
      temp+=1.0;
    }
  rndm5_store.u[rndm5_store.i]=temp;
  if (--rndm5_store.i<0) rndm5_store.i=96;
  if (--rndm5_store.j<0) rndm5_store.j=96;
  rndm5_store.c-=rndm5_store.cd;
  if (rndm5_store.c<0.0)
    {
      rndm5_store.c+=rndm5_store.cm;
    }
  temp-=rndm5_store.c;
  if (temp<0.0)
    {
      return temp+1.0;
    }
  else
    {
/*      if (temp>0.0) */
        return temp;
    }
/*}*/
}


void RNDM::rndmst6(int i)
{
  rndm6_store.i=i;
}

float RNDM::rndm6()
{

  if((rndm6_store.i*=48828125)<0)
    {
      rndm6_store.i-=2147483647;
      rndm6_store.i--;
    }
  return (float)rndm6_store.i*0.46566129e-9*(1.0-RNDM_EPS);
}


void RNDM::rndmst7(unsigned long i)
{
  //  std::cout << " test of the random generator rndm7 " << std::endl;

  rndm7_store.i= 0L;
  unsigned long a_coeff;
  unsigned long c_coeff;
  unsigned long test7[11];

  // 32-bits machine 
#ifndef COMPUTER64b

  test7[0] = 0x3c6ef35f;
  test7[1] = 0x47502932;
  test7[2] = 0xd1ccf6e9;
  test7[3] = 0xaaf95334;
  test7[4] = 0x6252e503;
  test7[5] = 0x9f2ec686;
  test7[6] = 0x57fe6c2d;
  test7[7] = 0xa3d95fa8;
  test7[8] = 0x81fdbee7;
  test7[9] = 0x94f0af1a;
  test7[10] = 0xcbf633b1;

  a_coeff = 1664525L;
  c_coeff = 1013904223L;
  if ( rndm_test7(a_coeff , c_coeff, test7) )
    {
      std::cout << " random generator successful test : 32-bits computer " << std::endl;      
      // following factor is 1/(2.exp(32))
      rndm7_modulo_dividing_factor_ = 0.232830643654e-9;
    }
  else
    {
      std::cerr << " RNDM::rndmst7 : test for random generator failed, are we really working on a  32-bits computer ? " << std::endl;
      exit(0);
    }
#else
  test7[0] = 0x1;
  test7[1] = 0x5851f42d4c957f2e;
  test7[2] = 0xc0b18ccf4e252d17;
  test7[3] = 0xcbb5f646404a560c;
  test7[4] = 0xc7033129d2bd141d;
  test7[5] = 0x30705b042917ec1a;
  test7[6] = 0x20fd5db43a776693;
  test7[7] = 0x9a8b7f78da6ef4d8;
  test7[8] = 0x502959d812b031f9;
  test7[9] = 0xab894868b3b04fc6;
  test7[10]= 0x6c0356a743cf3fcf;


  a_coeff = 6364136223846793005L;
  c_coeff = 1L;
  if ( rndm_test7(a_coeff , c_coeff, test7) )
    {
      std::cout << " random generator successful test : 64-bits computer " << std::endl;
      // following factor is 1/(2.exp(64))
      rndm7_modulo_dividing_factor_ = 0.54210108624e-19;
    }
  else
    {
      std::cerr << " RNDM::rndmst7 : test for random generator failed, are we really working on a  64-bits computer ? maybe you should desactive the define COMPUTER64b (CMT requirements) " << std::endl;
      exit(0);
    }
#endif

  rndm7_a_coeff_ = a_coeff;
  rndm7_c_coeff_ = c_coeff;
  rndm7_store.i=i;
}

bool RNDM::rndm_test7(unsigned long a_coeff, unsigned long c_coeff, const unsigned long* test_values)
{
  int k;
  unsigned long test_i = 0;
  bool test = true;
  for (k=0; k< 11; k++)
    {
       test_i =a_coeff* test_i + c_coeff;
       if (test_i != test_values[k])
	 {
	   test = false;
	   printf(" k= %d random value= %lx should be %lx \n", k, test_i, test_values[k]);
	   break;
	 }
    }
  return test;
}

float RNDM::rndm7()
{
  counterRndm7_++;
  rndm7_store.i=rndm7_a_coeff_*rndm7_store.i+rndm7_c_coeff_;
/* additional factor to ensure rndm7!=1.0 if float */
  return  (float)rndm7_store.i*rndm7_modulo_dividing_factor_*(1.0-RNDM_EPS);
}


void RNDM::rndmst8(int idummy)
{

  m1_ =259200;
  ia1_ =7141;
  ic1_ =54773;
  rm1_ =1.0/259200.0;
  m2_ =134456;
  ia2_ =8121;
  ic2_ =28411;
  rm2_ =1.0/134456.0;
  m3_ =243000;
  ia3_ =4561;
  ic3_ =51349;
  iff_ =0;
  rm3_ =1.0/243000.0;


  if (idummy>0) idummy=-idummy;
  ix1_=(ic1_-idummy) % m1_;
  ix1_=(ia1_*ix1_+ic1_) % m1_;
  ix2_=ix1_ % m2_;
  ix1_=(ia1_*ix1_+ic1_) % m1_;
  ix3_=ix1_ % m3_;
  int k;
  for (k=0;k<97;k++)
    {
      ix1_=(ia1_*ix1_+ic1_) % m1_;
      ix2_=(ia2_*ix2_+ic2_) % m2_;
      rndm8_store.r[k]=(ix1_+ix2_*rm2_)*rm1_;
    }
  rndm8_store.ix1=ix1_;
  rndm8_store.ix2=ix2_;
  rndm8_store.ix3=ix3_;
}

float RNDM::rndm8()
{
  float help;

  rndm8_store.ix1=(ia1_*rndm8_store.ix1+ic1_) % m1_;
  rndm8_store.ix2=(ia2_*rndm8_store.ix2+ic2_) % m2_;
  rndm8_store.ix3=(ia3_*rndm8_store.ix3+ic3_) % m3_;
  int k =(int)((float)(97*rndm8_store.ix3)*rm3_);
  help=rndm8_store.r[k];
  rndm8_store.r[k]=(rndm8_store.ix1+rndm8_store.ix2*rm2_)*rm1_;
  return help;
}

float RNDM::expdev()
{
  return -log(1.0-rndm());
}

float RNDM::exp_dev()
{
    return -log(1.0-rndm2());
}

float RNDM::gasdev()
{
  float r;
  if (iset_ == 0)
    {
      for (;;)
	{
          v1_ = 2.0*rndm()-1.0;
          v2_ = 2.0*rndm()-1.0;
          r=v1_*v1_+v2_*v2_;
          if ((r<=1.0) && (r!=0.0))
	    {
	      break;
	    }
        }
      iset_ = 1;

               r=sqrt(-2.0*log((double)r)/r);
      v1_ *= r;
      v2_ *= r;
      return v1_;
    }
  else
    {
      iset_ = 0;
      return v2_;
    }
}

float RNDM::rndm_sincos(float *c)
{
    const float twopi=2.0*PI;
    float tmp;
    tmp=rndm();
    *c=cos(twopi*tmp);
    if (tmp>0.5)
	return sqrt(1.0- *c * *c);
    else
	return -sqrt(1.0- *c * *c);
}

// get a vector containing the sequence of integers from 1 to maxInt 
// randomly shuffled 
void RNDM::getShuffledIntegerSequence(int maxInt, std::vector<unsigned long int>& vec)
{
  int k, j;
  unsigned long int itemp;
  float u;
  vec.clear();
  vec.resize(maxInt);
  for (k = 0; k < maxInt; k++) vec[k] = k + 1;
  j = maxInt;
  while (j > 1)
    {
      u = rndm();
      k = (int)floor(j*u) +1;
      itemp = vec[k-1];
      vec[k-1] = vec[j-1];
      vec[j-1] = itemp;
      j--;
    }
}

int RNDM::rndm_save()
{
    FILE *file;
    file=fopen("rndm.save","w");
    fwrite(&rndm0_store,sizeof(rndm0_store),1,file);
    fwrite(&rndm1_store,sizeof(rndm1_store),1,file);
    fwrite(&rndm2_store,sizeof(rndm2_store),1,file);
    fwrite(&rndm5_store,sizeof(rndm5_store),1,file);
    fwrite(&rndm6_store,sizeof(rndm6_store),1,file);
    fwrite(&rndm7_store,sizeof(rndm7_store),1,file);
    fwrite(&rndm8_store,sizeof(rndm8_store),1,file);
    fclose(file);
    return 0;
}

int RNDM::rndm_load()
{
    FILE *file;	
    if ( (file = fopen("rndm.save","r")) )
      {
        if( !fread(&rndm0_store,sizeof(rndm0_store),1,file)) {fclose(file); return 0;}
        fread(&rndm1_store,sizeof(rndm1_store),1,file);
        fread(&rndm2_store,sizeof(rndm2_store),1,file);
        fread(&rndm5_store,sizeof(rndm5_store),1,file);
        fread(&rndm6_store,sizeof(rndm6_store),1,file);
        fread(&rndm7_store,sizeof(rndm7_store),1,file);
        fread(&rndm8_store,sizeof(rndm8_store),1,file);
        fclose(file);
	std::cout << "  RNDM::rndm_load : the seed is now : " << rndm7_store.i << std::endl;
	return 1;
      }
    std::cout << "  RNDM::rndm_load : no rndm_load file, the seed is : " << rndm7_store.i << std::endl;
    return 0;
}


