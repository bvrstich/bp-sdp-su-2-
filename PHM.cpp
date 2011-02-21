#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;

#include "include.h"

int PHM::counter = 0;

int **PHM::ph2s;
int **PHM::s2ph;

double **PHM::_6j;

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1 of dimension M*M/4.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PHM::PHM(int M,int N) : BlockMatrix(2) {
   
   this->N = N;
   this->M = M;

   //set the dimension of the blocks
   this->setMatrixDim(0,M*M/4,1);
   this->setMatrixDim(1,M*M/4,3);

   if(counter == 0)
      constr_lists();

   ++counter;

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks of dimension M*M/4 and copies the content of phm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){

   this->N = phm_c.gN();
   this->M = phm_c.gM();

   if(counter == 0)
      constr_lists();

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){

   if(counter == 1){

      delete [] s2ph[0];
      delete [] s2ph;

      for(int i = 0;i < M*M/4;++i)
         delete [] ph2s[i];

      delete [] ph2s;

      for(int S = 0;S < 2;++S)
         delete [] _6j[S];

      delete [] _6j;

   }

   --counter;

}

/**
 * Allocate and fill the lists needed in this class
 */
void PHM::constr_lists(){

   s2ph = new int * [M/2];
   s2ph[0] = new int [M*M/4];

   for(int i = 1;i < M/2;++i)
      s2ph[i] = s2ph[i - 1] + M/2;

   //allocation of ph2s
   ph2s = new int * [M*M/4];

   for(int i = 0;i < M*M/4;++i)
      ph2s[i] = new int [2];

   //initialisation of the two arrays
   int teller = 0;

   for(int a = 0;a < M/2;++a)
      for(int b = 0;b < M/2;++b){

         s2ph[a][b] = teller;

         ph2s[teller][0] = a;
         ph2s[teller][1] = b;

         ++teller;

      }

   //allocate
   _6j = new double * [2];

   for(int S = 0;S < 2;++S)
      _6j[S] = new double [2];

   //initialize
   _6j[0][0] = -0.5;
   _6j[0][1] = 0.5;
   _6j[1][0] = 0.5;
   _6j[1][1] = 1.0/6.0;

}

/**
 * access the elements of the matrix in sp mode, 
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(S,i,j)
 */
double &PHM::operator()(int S,int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

/**
 * access the elements of the matrix in sp mode: read only mode
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(S,i,j)
 */
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   for(int S = 0;S < phm_p.gnr();++S){

      output << S << "\t" << phm_p.gdim(S) << "\t" << phm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(S);++i)
         for(int j = 0;j < phm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

               << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * @return number of particles
 */
int PHM::gN() const{

   return N;

}

/**
 * @return number of single particle oribals
 */
int PHM::gM() const{

   return M;

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm(1.0/(N - 1.0),tpm);

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < this->gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            //tp part
            (*this)(S,i,j) = -_6j[S][0]*tpm(0,a,d,c,b) - 3.0*_6j[S][1]*tpm(1,a,d,c,b);

            //norm
            if(a == d)
               (*this)(S,i,j) *= std::sqrt(2.0);

            if(c == b)
               (*this)(S,i,j) *= std::sqrt(2.0);

            //sp part
            if(b == d)
               (*this)(S,i,j) += spm(a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void PHM::bar(const PPHM &pphm){

   int a,b,c,d;

   double ward,hard;

   for(int S = 0;S < 2;++S){//loop over spinblocks of PHM

      for(int i = 0;i < this->gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            (*this)(S,i,j) = 0.0;

            //first the S = 1/2 block of the PPHM matrix
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * _6j[S][S_ab] * _6j[S][S_de];

                  for(int l = 0;l < M/2;++l){

                     hard = ward * pphm(0,S_ab,l,a,b,S_de,l,c,d);

                     //norms
                     if(l == a)
                        hard *= std::sqrt(2.0);

                     if(l == c)
                        hard *= std::sqrt(2.0);

                     (*this)(S,i,j) += hard;

                  }

               }

            //then the S = 3/2 block
            if(S == 1)
               for(int l = 0;l < M/2;++l)
                  (*this)(S,i,j) += 4.0/3.0 * pphm(1,1,l,a,b,1,l,c,d);

         }
      }

   }

   this->symmetrize();

}

/**
 * Will print out a spin uncoupled version of the PHM (*this) in uncoupled sp coordinates to the file with name filename
 */
void PHM::uncouple(const char *filename){

   ofstream output(filename);

   output.precision(10);

   int a,b,c,d;
   int s_a,s_b,s_c,s_d;

   double ward_0;
   double ward_1;

   //the uncoupled row vector i
   for(int alpha = 0;alpha < M;++alpha)
      for(int beta = 0;beta < M;++beta){

         //watch it now, 0 is down, 1 is up.
         a = alpha/2;
         s_a = alpha%2;

         b = beta/2;
         s_b = beta%2;

         //the collom vector j
         for(int gamma = alpha;gamma < M;++gamma)
            for(int delta = 0;delta < M;++delta){

               c = gamma/2;
               s_c = gamma%2;

               d = delta/2;
               s_d = delta%2;

               //eerst S = 0 stuk
               if(s_a == s_b && s_c == s_d)
                  ward_0 = 0.5 * (*this)(0,a,b,c,d);
               else
                  ward_0 = 0.0;

               //dan S = 1 stuk
               if( (s_a != s_b) && (s_c != s_d) && (s_a == s_c) )
                  ward_1 = (*this)(1,a,b,c,d);
               else if(s_a == s_b && s_c == s_d){

                  if(s_a == s_c)
                     ward_1 = 0.5 * (*this)(1,a,b,c,d);
                  else
                     ward_1 = -0.5 * (*this)(1,a,b,c,d);

               }
               else
                  ward_1 = 0.0;

               output << alpha << "\t" << beta << "\t" << gamma << "\t" << delta << "\t" << ward_0 + ward_1 << std::endl;

            }

      }

}

/**
 * Input from file, with sp indices as row and column indices
 * @param filename Name of the inputfile
 */
void PHM::in_sp(const char *filename){

   ifstream input(filename);

   int a,b,c,d;
   int S;

   double value;

   int i,j;

   (*this) = 0;

   while(input >> S >> a >> b >> c >> d >> value){

      i = s2ph[a][b];
      j = s2ph[c][d];

      (*this)(S,i,j) = value;

   }

   this->symmetrize();

}
