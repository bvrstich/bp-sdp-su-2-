#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;

#include "include.h"

int PPHM::counter = 0;

int ***PPHM::pph2s;
int *****PPHM::s2pph;

double **PPHM::_6j;

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
PPHM::PPHM(int M,int N) : BlockMatrix(2) {

   this->N = N;
   this->M = M;

   //set the dimension and the degeneracies of the blocks
   this->setMatrixDim(0,M*M*M/8,2);//S=1/2 block
   this->setMatrixDim(1,M*M*(M - 2)/16,4);//S=3/2 block

   if(counter == 0)//make the lists
      construct_lists();

   ++counter;

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks, on for S=1/2 and one for S=3/2, and copies the content of the pphm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and pph basis.
 * @param pphm_c PPHM object to be copied into (*this)
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c) {

   this->N = pphm_c.gN();
   this->M = pphm_c.gM();

   if(counter == 0)
      construct_lists();

   ++counter;

}

/**
 * Destructor, if counter = 1 the lists will be deallocated.
 */
PPHM::~PPHM(){

   if(counter == 1){

      //first delete S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < M/2;++a){

            for(int b = 0;b < M/2;++b)
               delete [] s2pph[0][S_ab][a][b];

            delete [] s2pph[0][S_ab][a];

         }

         delete [] s2pph[0][S_ab];

      }

      //then the S = 3/2 part
      for(int a = 0;a < M/2;++a){

         for(int b = 0;b < M/2;++b)
            delete [] s2pph[1][0][a][b];

         delete [] s2pph[1][0][a];

      }

      delete [] s2pph[1][0];

      for(int S = 0;S < 2;++S)
         delete [] s2pph[S];

      delete [] s2pph;

      //now delete pph2s 
      for(int S = 0;S < 2;++S){

         for(int i = 0;i < this->gdim(S);++i)
            delete [] pph2s[S][i];

         delete [] pph2s[S];

      }

      delete [] pph2s;

      for(int S = 0;S < 2;++S)
         delete [] _6j[S];

      delete [] _6j;

   }

   --counter;

}


void PPHM::construct_lists(){

   //first allocation
   pph2s = new int ** [2];//two total spinblocks

   for(int S = 0;S < 2;++S){

      pph2s[S] = new int * [this->gdim(S)];//dimension of the blocks

      for(int i = 0;i < this->gdim(S);++i)
         pph2s[S][i] = new int [4];//amount of information stored for an index: (S_ab,a,b,c)

   }

   s2pph = new int **** [2];//two spinblocks

   s2pph[0] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

   for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

      s2pph[0][S_ab] = new int ** [M/2];

      for(int a = 0;a < M/2;++a){

         s2pph[0][S_ab][a] = new int * [M/2];

         for(int b = 0;b < M/2;++b)
            s2pph[0][S_ab][a][b] = new int [M/2];

      }

   }

   s2pph[1] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

   s2pph[1][0] = new int ** [M/2];

   for(int a = 0;a < M/2;++a){//loop and allocate

      s2pph[1][0][a] = new int * [M/2];

      for(int b = 0;b < M/2;++b)
         s2pph[1][0][a][b] = new int [M/2];

   }

   //initialize the lists
   int teller = 0;

   //S = 1/2 S_ab = 0: a <= b, c
   for(int a = 0;a < M/2;++a)
      for(int b = a;b < M/2;++b)
         for(int c = 0;c < M/2;++c){

            s2pph[0][0][a][b][c] = teller;

            pph2s[0][teller][0] = 0;//S_ab

            pph2s[0][teller][1] = a;
            pph2s[0][teller][2] = b;
            pph2s[0][teller][3] = c;

            ++teller;

         }

   //S = 1/2, S_ab = 1, a < b ,c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = 0;c < M/2;++c){

            s2pph[0][1][a][b][c] = teller;

            pph2s[0][teller][0] = 1;//S_ab

            pph2s[0][teller][1] = a;
            pph2s[0][teller][2] = b;
            pph2s[0][teller][3] = c;

            ++teller;

         }

   //re-init teller for block S = 3/2
   teller = 0;

   //S = 3/2, S_ab = 1, a < b, c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = 0;c < M/2;++c){

            s2pph[1][0][a][b][c] = teller;

            pph2s[1][teller][0] = 1;//S_ab

            pph2s[1][teller][1] = a;
            pph2s[1][teller][2] = b;
            pph2s[1][teller][3] = c;

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
 * @return nr of particles
 */
int PPHM::gN() const{

   return N;

}

/**
 * @return nr of sp orbitals
 */
int PPHM::gM() const{

   return M;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S The block index, when == 0 then access the block S = 1/2, for block == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the pph row index i together with b, c and S_ab in block S
 * @param b second sp index that forms the pph row index i together with a, c and S_ab in block S
 * @param c third sp index that forms the pph row index i together with a, b and S_ab in block S
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the pph column index j together with e, z and S_de in block S
 * @param e second sp index that forms the pph column index j together with d, z and S_de in block S
 * @param z third sp index that forms the pph column index j together with d, e and S_de in block S
 * @return the number on place PPHM(S,i,j) with the right phase.
 */
double PPHM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int i,j;

   int phase_i = get_inco(S,S_ab,a,b,c,i);

   if(phase_i == 0)
      return 0;

   int phase_j = get_inco(S,S_de,d,e,z,j);

   if(phase_j == 0)
      return 0;

   return phase_i*phase_j* (*this)(S,i,j);

}

/** 
 * Static member function that gets the pph-index and phase corresponding to the sp indices S,S_ab,a,b,c.
 * @param S block index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param S_ab intermediate spincoupling of a and b. = 0 or 1
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
int PPHM::get_inco(int S,int S_ab,int a,int b,int c,int &i){

   if(S == 0){//S = 1/2

      if(S_ab == 0){//symmetric in spatial sp's

         if(a <= b)
            i = s2pph[0][0][a][b][c];
         else
            i = s2pph[0][0][b][a][c];

         return 1;

      }
      else{//antisymmetric in spatial sp's

         if(a == b)
            return 0;

         if(a < b){

            i = s2pph[0][1][a][b][c];

            return 1;

         }
         else{

            i = s2pph[0][1][b][a][c];

            return -1;

         }

      }

   }
   else{//S = 3/2

      if(S_ab == 0)
         return 0;

      if(a == b)
         return 0;

      if(a < b){

         i = s2pph[1][0][a][b][c];

         return 1;

      }
      else{

         i = s2pph[1][0][b][a][c];

         return -1;

      }

   }

}

/**
 * The spincoupled T2 map, maps a TPM onto a PPHM object. See notes for more info
 * @param tpm Input TPM matrix
 */
void PPHM::T(const TPM &tpm){

   SPM spm(1.0/(N - 1.0),tpm);

   int a,b,c,d,e,z;
   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

   for(int i = 0;i < this->gdim(0);++i){

      S_ab = pph2s[0][i][0];

      a = pph2s[0][i][1];
      b = pph2s[0][i][2];
      c = pph2s[0][i][3];

      sign_ab = 1 - 2*S_ab;

      norm_ab = 1.0;

      if(a == b)
         norm_ab /= std::sqrt(2.0);

      for(int j = i;j < this->gdim(0);++j){

         S_de = pph2s[0][j][0];

         d = pph2s[0][j][1];
         e = pph2s[0][j][2];
         z = pph2s[0][j][3];

         sign_de = 1 - 2*S_de;

         norm_de = 1.0;

         if(d == e)
            norm_de /= std::sqrt(2.0);


         //start the map:
         (*this)(0,i,j) = 0.0;

         //tp(1)
         if(c == z)
            if(S_ab == S_de)
               (*this)(0,i,j) += tpm(S_ab,a,b,d,e);

         if(a == d){

            //sp(1) first term
            if(b == e)
               if(S_ab == S_de)
                  (*this)(0,i,j) += norm_ab * norm_de * spm(c,z);

            //tp(2)
            double ward = 0.0;

            for(int J = 0;J < 2;++J)
               for(int Z = 0;Z < 2;++Z)
                  ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,c,e,z,b);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= ward;

         }

         if(b == d){

            //sp(1) second term
            if(a == e)
               if(S_ab == S_de)
                  (*this)(0,i,j) += sign_ab * norm_ab * norm_de * spm(c,z);

            //tp(3)
            double ward = 0.0;

            for(int J = 0;J < 2;++J)
               for(int Z = 0;Z < 2;++Z)
                  ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,c,e,z,a);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * ward;

         }

         //tp(4)
         if(a == e){

            double ward = 0.0;

            for(int J = 0;J < 2;++J)
               for(int Z = 0;Z < 2;++Z)
                  ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,c,d,z,b);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_de * ward;

         }

         //tp(5)
         if(b == e){

            double ward = 0.0;

            for(int J = 0;J < 2;++J)
               for(int Z = 0;Z < 2;++Z)
                  ward += (2*J + 1.0) * (2*Z + 1.0) * _6j[J][S_ab] * _6j[J][S_de] * _6j[J][Z] * tpm(Z,c,d,z,a);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * sign_de * ward;

         }

      }

   }

   //the easier S = 3/2 part:
   for(int i = 0;i < this->gdim(1);++i){

      a = pph2s[1][i][1];
      b = pph2s[1][i][2];
      c = pph2s[1][i][3];

      for(int j = i;j < this->gdim(1);++j){

         d = pph2s[1][j][1];
         e = pph2s[1][j][2];
         z = pph2s[1][j][3];

         //init
         (*this)(1,i,j) = 0.0;

         //tp(1)
         if(c == z)
            (*this)(1,i,j) += tpm(1,a,b,d,e);

         if(a == d){

            //sp(1)
            if(b == e)
               (*this)(1,i,j) += spm(c,z);

            //tp(2)
            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,c,e,z,b);

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) -= ward;

         }

         //tp(3)
         if(b == d){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,c,e,z,a);

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) += ward;

         }

         //tp(5)
         if(b == e){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * _6j[1][Z] * tpm(Z,c,d,z,a);

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) -= ward;

         }

      }

   }

   this->symmetrize();

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int S = 0;S < pphm_p.gnr();++S){

      output << S << "\t" << pphm_p.gdim(S) << "\t" << pphm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(S);++i)
         for(int j = 0;j < pphm_p.gdim(S);++j){

            output << S << "\t" << i << "\t" << j << "\t|\t" << 
            
               pphm_p.pph2s[S][i][0] << "\t" << pphm_p.pph2s[S][i][1] << "\t" << pphm_p.pph2s[S][i][2] << "\t" << pphm_p.pph2s[S][i][3] << 

               "\t" << pphm_p.pph2s[S][j][0] << "\t" << pphm_p.pph2s[S][j][1] << "\t" << pphm_p.pph2s[S][j][2] << "\t" << pphm_p.pph2s[S][j][3] 
               
               << "\t" << pphm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * Input from file, with sp indices as row and column indices
 * @param filename Name of the inputfile
 */
void PPHM::in_sp(const char *filename){

   ifstream input(filename);

   int S_ab,S_de;
   int a,b,c,d,e,z;
   int S;

   double value;

   int i,j;

   (*this) = 0;

   while(input >> S >> S_ab >> a >> b >> c >> S_de >> d >> e >> z >> value){

      if(S == 0){

         i = s2pph[S][S_ab][a][b][(-c + M/2)%(M/2)];
         j = s2pph[S][S_de][d][e][(-z + M/2)%(M/2)];

      }
      else{

         i = s2pph[S][0][a][b][(-c + M/2)%(M/2)];
         j = s2pph[S][0][d][e][(-z + M/2)%(M/2)];

      }

      (*this)(S,i,j) = value;

   }

   this->symmetrize();

}
