#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;

#include "include.h"

int DPM::counter = 0;

int ***DPM::dp2s;
int *****DPM::s2dp;

double **DPM::_6j;

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
DPM::DPM(int M,int N) : BlockMatrix(2) {

   this->N = N;
   this->M = M;

   //set the dimension and the degeneracies of the blocks
   this->setMatrixDim(0,M/2*(M/2 - 1) + M/2*(M/2 - 1)*(M/2 - 2)/3,2);
   this->setMatrixDim(1,M/2*(M/2 - 1)*(M/2 - 2)/6,4);

   if(counter == 0)//make the lists
      construct_lists();

   ++counter;

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks, on for S=1/2 and one for S=3/2, and copies the content of the dpm_c blocks into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and dp basis.
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) {

   this->N = dpm_c.gN();
   this->M = dpm_c.gM();

   if(counter == 0)
      construct_lists();

   ++counter;

}

/**
 * destructor: if counter == 1 the memory for the static lists dp2s en s2dp will be deleted.
 */
DPM::~DPM(){

   if(counter == 1){

      //first delete S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < M/2;++a){

            for(int b = 0;b < M/2;++b)
               delete [] s2dp[0][S_ab][a][b];

            delete [] s2dp[0][S_ab][a];

         }

         delete [] s2dp[0][S_ab];

      }

      //then the S = 3/2 part
      for(int a = 0;a < M/2;++a){

         for(int b = 0;b < M/2;++b)
            delete [] s2dp[1][0][a][b];

         delete [] s2dp[1][0][a];

      }

      delete [] s2dp[1][0];

      for(int S = 0;S < 2;++S)
         delete [] s2dp[S];

      delete [] s2dp;

      //now delete dp2s 
      for(int S = 0;S < 2;++S){

         for(int i = 0;i < this->gdim(S);++i)
            delete [] dp2s[S][i];

               delete [] dp2s[S];

      }

      delete [] dp2s;

      for(int S = 0;S < 2;++S)
         delete [] _6j[S];

      delete [] _6j;

   }

   --counter;

}

/** 
 * Function that allocates and initializes the lists needed in the program, called when the first DPM object is constructed,
 */
void DPM::construct_lists(){

   //first allocation
   dp2s = new int ** [2];//two total spinblocks

   for(int S = 0;S < 2;++S){

      dp2s[S] = new int * [this->gdim(S)];//dimension of the blocks

      for(int i = 0;i < this->gdim(S);++i)
         dp2s[S][i] = new int [4];//amount of information stored for an index: (S_ab,a,b,c)

   }

   s2dp = new int **** [2];//two spinblocks

   s2dp[0] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

   for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

      s2dp[0][S_ab] = new int ** [M/2];

      for(int a = 0;a < M/2;++a){

         s2dp[0][S_ab][a] = new int * [M/2];

         for(int b = 0;b < M/2;++b)
            s2dp[0][S_ab][a][b] = new int [M/2];

      }

   }

   s2dp[1] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

   s2dp[1][0] = new int ** [M/2];

   for(int a = 0;a < M/2;++a){//loop and allocate

      s2dp[1][0][a] = new int * [M/2];

      for(int b = 0;b < M/2;++b)
         s2dp[1][0][a][b] = new int [M/2];

   }

   //initialize the lists
   int teller = 0;

   //first S == 1/2, S_ab == 0 and a == b != c
   for(int a = 0;a < M/2;++a){

      for(int c = 0;c < a;++c){

         s2dp[0][0][a][a][c] = teller;

         dp2s[0][teller][0] = 0;//S_ab

         dp2s[0][teller][1] = a;
         dp2s[0][teller][2] = a;
         dp2s[0][teller][3] = c;

         ++teller;

      }

      for(int c = a + 1;c < M/2;++c){

         s2dp[0][0][a][a][c] = teller;

         dp2s[0][teller][0] = 0;//S_ab

         dp2s[0][teller][1] = a;
         dp2s[0][teller][2] = a;
         dp2s[0][teller][3] = c;

         ++teller;

      }

   }

   //S and S_ab the same but a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[0][0][a][b][c] = teller;

            dp2s[0][teller][0] = 0;//S_ab

            dp2s[0][teller][1] = a;
            dp2s[0][teller][2] = b;
            dp2s[0][teller][3] = c;

            ++teller;

         }

   //S == 0, S_ab == 1, a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[0][1][a][b][c] = teller;

            dp2s[0][teller][0] = 1;//S_ab

            dp2s[0][teller][1] = a;
            dp2s[0][teller][2] = b;
            dp2s[0][teller][3] = c;

            ++teller;

         }

   //re-init teller
   teller = 0;

   //S == 1, S_ab == 1, a != b != c
   for(int a = 0;a < M/2;++a)
      for(int b = a + 1;b < M/2;++b)
         for(int c = b + 1;c < M/2;++c){

            s2dp[1][0][a][b][c] = teller;

            dp2s[1][teller][0] = 1;//S_ab

            dp2s[1][teller][1] = a;
            dp2s[1][teller][2] = b;
            dp2s[1][teller][3] = c;

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
 * @return number of particles
 */
int DPM::gN() const{

   return N;

}

/**
 * @return number of single particle oribals
 */
int DPM::gM() const{

   return M;

}

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymemtry relations are automatically accounted for:\n\n
 * DPM(S,S_ab,a,b,c,S_de,d,e,f) = sum_S_ac (some terms dependent on spin) DPM(S,S_ac,a,c,b,S_de,d,e,f) etc...
 * @param S The block index, when == 0 then access the block S = 1/2, for block == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the dp row index i together with b, c and S_ab in block S
 * @param b second sp index that forms the dp row index i together with a, c and S_ab in block S
 * @param c third sp index that forms the dp row index i together with a, b and S_ab in block S
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the dp column index j together with e, z and S_de in block S
 * @param e second sp index that forms the dp column index j together with d, z and S_de in block S
 * @param z third sp index that forms the dp column index j together with d, e and S_de in block S
 * @return the number on place DPM(S,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(S,S_ab,a,b,c,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(S,S_de,d,e,z,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += coef_i[I] * coef_j[J] * (*this)(S,i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

}

/** 
 * Static member function that gets the dp-indices and their coefficients of the sp indices S,S_ab,a,b,c.
 * @param S block index of the state
 * @param S_ab intermediate spincoupling of a and b.
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int S,int S_ab,int a,int b,int c,int *i,double *coef) const{

   //they cannot all be equal
   if(a == b && b == c)
      return 0;

   if(S == 0){//spin 1/2 block:

      //if normal basis:
      if(a == b){

         if(S_ab == 1)//spin has to be zero for a == b
            return 0;

         i[0] = s2dp[0][0][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else if (a < b && b < c){

         i[0] = s2dp[0][S_ab][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else{//anomal basis:

         int min,max,phase;

         //first order a and b for code saving reasons
         if(a < b){

            min = a;
            max = b;

            phase = 1;

         }
         else{

            min = b;
            max = a;

            phase = 1 - 2*S_ab;

            if(c > max){//we still have one simple dim = 1 term left: b < a < c

               i[0] = s2dp[0][S_ab][b][a][c];
               coef[0] = phase;

               return 1;

            }

         }

         //now we have four possibilities left:
         //don't forget to multiply every result by phase to get the right a and b for min and max!
         // 1) c < min < max
         // 2) c == min < max
         // 3) min < c < max
         // 4) min < max == c
         if(c < min){//c < min < max

            //the S_ca == 0 part:
            i[0] = s2dp[0][0][c][min][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            //the S_ca == 1 part:
            i[1] = s2dp[0][1][c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else if(c == min){//c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if a == c.

            i[0] = s2dp[0][0][c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * _6j[0][S_ab];

            return 1;

         }
         else if(c < max){//min < c < max

            //S_ac == 0 part:
            i[0] = s2dp[0][0][min][c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * _6j[0][S_ab];

            //S_ac == 1 part:
            i[1] = s2dp[0][1][min][c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * _6j[1][S_ab];

            return 2;

         }
         else{// min < c == max: also a 1 dim list, S_bc can only be 0 if b == c

            i[0] = s2dp[0][0][max][c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) *_6j[0][S_ab];

            return 1;

         }

      }

   }
   else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

      //only S_ab == 1 can couple to 3/2's.
      if(S_ab == 0)
         return 0;

      //if any of the sp orbs are equal, antisymmetry leads to zero:
      if(a == b || b == c || c == a)
         return 0;

      if(a < b){

         if(b < c){//a < b < c

            i[0] = s2dp[1][0][a][b][c];
            coef[0] = 1;

         }
         else if(c < a){//c < a < b

            i[0] = s2dp[1][0][c][a][b];
            coef[0] = 1;

         }
         else{//a < c < b

            i[0] = s2dp[1][0][a][c][b];
            coef[0] = -1;

         }

      }
      else{//b < a

         if(a < c){//b < a < c

            i[0] = s2dp[1][0][b][a][c];
            coef[0] = -1;

         }
         else if(c < b){//c < b < a

            i[0] = s2dp[1][0][c][b][a];
            coef[0] = -1;

         }
         else{//b < c < a

            i[0] = s2dp[1][0][b][c][a];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * The spincoupled T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,const TPM &tpm){

   //make sp matrix out of tpm
   SPM spm(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int a,b,c,d,e,z;
   int S_ab,S_de;

   int sign_ab,sign_de;

   double norm_ab,norm_de;

   double hard;

   //start with the S = 1/2 block, this is the most difficult one:
   for(int i = 0;i < this->gdim(0);++i){

      S_ab = dp2s[0][i][0];

      a = dp2s[0][i][1];
      b = dp2s[0][i][2];
      c = dp2s[0][i][3];

      sign_ab = 1 - 2*S_ab;

      norm_ab = 1.0;

      if(a == b)
         norm_ab /= std::sqrt(2.0);

      for(int j = i;j < this->gdim(0);++j){

         S_de = dp2s[0][j][0];

         d = dp2s[0][j][1];
         e = dp2s[0][j][2];
         z = dp2s[0][j][3];

         sign_de = 1 - 2*S_de;

         norm_de = 1.0;

         if(d == e)
            norm_de /= std::sqrt(2.0);

         hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * _6j[S_ab][S_de];

         //init
         (*this)(0,i,j) = 0.0;

         //the np part
         if(i == j)
            (*this)(0,i,j) = ward;

         //other parts are a bit more difficult.
         if(c == z){

            if(S_ab == S_de){

               //tp(1)
               (*this)(0,i,j) += A * tpm(S_ab,a,b,d,e);

               //sp(1) first term
               if(b == e)
                  (*this)(0,i,j) -= norm_ab * norm_de * spm(a,d);

               //sp(2) first term
               if(a == e)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_de * spm(b,d);

               //sp(4) first term
               if(b == d)
                  (*this)(0,i,j) -= sign_de * norm_ab * norm_de * spm(a,e);

               //sp(5) first term
               if(a == d)
                  (*this)(0,i,j) -= norm_ab * norm_de * spm(b,e);

            }

         }

         if(b == z){

            //tp(2)
            if(a == c)
               (*this)(0,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);
            else
               (*this)(0,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);

            //sp(1) second term
            if(c == e)
               (*this)(0,i,j) -= sign_ab * sign_de * norm_ab * norm_de * hard * spm(a,d);

            //sp(3)
            if(a == e)
               (*this)(0,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(c,d);

            //sp(4) second term
            if(c == d)
               (*this)(0,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(a,e);

            //sp(6)
            if(a == d)
               (*this)(0,i,j) -= sign_ab * sign_de * norm_ab * norm_de * hard * spm(c,e);

         }

         if(a == z){

            //tp(3)
            if(b == c)
               (*this)(0,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);
            else
               (*this)(0,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);

            //sp(2) second term
            if(c == e)
               (*this)(0,i,j) -= sign_de * norm_ab * norm_de * hard * spm(b,d);

            //sp(5) second term
            if(c == d)
               (*this)(0,i,j) -= norm_ab * norm_de * hard * spm(b,e);

         }

         if(c == e){

            //tp(4)
            if(d == z)
               (*this)(0,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);
            else
               (*this)(0,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);

            //sp(7) first term
            if(b == d)
               (*this)(0,i,j) -= norm_ab * norm_de * sign_de * hard * spm(a,z);

            //sp(8) first term
            if(a == d)
               (*this)(0,i,j) -= norm_ab * norm_de * sign_ab * sign_de * hard * spm(b,z);

         }

         if(b == e){

            //tp(5)
            double hulp = 0.0;

            //sum over intermediate spin
            for(int Z = 0;Z < 2;++Z)
               hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,a,c,d,z);

            //correct for norms of the tpm
            if(a == c)
               hulp *= std::sqrt(2.0);

            if(d == z)
               hulp *= std::sqrt(2.0);

            (*this)(0,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

            //sp(7) second term
            if(c == d)
               (*this)(0,i,j) -= norm_ab * norm_de * hard * spm(a,z);

            //sp(9) first term
            if(a == d)
               if(S_ab == S_de)
                  (*this)(0,i,j) -= norm_ab * norm_de * spm(c,z);

         }

         if(a == e){

            //tp(6)
            double hulp = 0.0;

            //sum over intermediate spin
            for(int Z = 0;Z < 2;++Z)
               hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,b,c,d,z);

            if(b == c)
               hulp *= std::sqrt(2.0);

            if(d == z)
               hulp *= std::sqrt(2.0);

            (*this)(0,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            //sp(8) second term
            if(c == d)
               (*this)(0,i,j) -= sign_ab * norm_ab * norm_de * hard * spm(b,z);

            //sp(9) second term
            if(b == d)
               if(S_ab == S_de)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_de * spm(c,z);

         }

         if(c == d){

            //tp(7)
            if(e == z)
               (*this)(0,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);
            else
               (*this)(0,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);

         }

         if(b == d){

            //tp(8)
            double hulp = 0.0;

            //sum over intermediate spin
            for(int Z = 0;Z < 2;++Z)
               hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,a,c,e,z);

            if(a == c)
               hulp *= std::sqrt(2.0);

            if(e == z)
               hulp *= std::sqrt(2.0);

            (*this)(0,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

         }

         if(a == d){

            //tp(8)
            double hulp = 0.0;

            //sum over intermediate spin
            for(int Z = 0;Z < 2;++Z)
               hulp += (2*Z + 1.0) * _6j[Z][S_ab] * _6j[Z][S_de] * tpm(Z,b,c,e,z);

            if(b == c)
               hulp *= std::sqrt(2.0);

            if(e == z)
               hulp *= std::sqrt(2.0);

            (*this)(0,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

         }

      }
   }


   //then the S = 3/2 block, this should be easy, totally antisymmetrical 
   for(int i = 0;i < this->gdim(1);++i){

      a = dp2s[1][i][1];
      b = dp2s[1][i][2];
      c = dp2s[1][i][3];

      for(int j = i;j < this->gdim(1);++j){

         d = dp2s[1][j][1];
         e = dp2s[1][j][2];
         z = dp2s[1][j][3];

         (*this)(1,i,j) = 0.0;

         if(i == j)
            (*this)(1,i,j) += ward;

         if(c == z){

            //tp(1)
            (*this)(1,i,j) += A * tpm(1,a,b,d,e);

            //sp(1) first part
            if(b == e)
               (*this)(1,i,j) -= spm(a,d);

            //sp(4) first part
            if(b == d)
               (*this)(1,i,j) += spm(a,e);

            //sp(5)
            if(a == d)
               (*this)(1,i,j) -= spm(b,e);

         }

         if(b == z){

            //tp(2)
            (*this)(1,i,j) -= A * tpm(1,a,c,d,e);

            //sp(1) second part
            if(c == e)
               (*this)(1,i,j) += spm(a,d);

            //sp(4) second part
            if(c == d)
               (*this)(1,i,j) -= spm(a,e);

            //sp(6)
            if(a == d)
               (*this)(1,i,j) += spm(c,e);

         }

         if(c == e){

            //tp(4)
            (*this)(1,i,j) -= A * tpm(1,a,b,d,z);

            //sp(7) first part
            if(b == d)
               (*this)(1,i,j) -= spm(a,z);

            //sp(8) first part
            if(a == d)
               (*this)(1,i,j) += spm(b,z);

         }

         if(b == e){

            //tp(5)
            (*this)(1,i,j) += A * tpm(1,a,c,d,z);

            //sp(7) second part
            if(c == d)
               (*this)(1,i,j) += spm(a,z);

            //sp(9) first part
            if(a == d)
               (*this)(1,i,j) -= spm(c,z);

         }

         //tp(7)
         if(c == d)
            (*this)(1,i,j) += A * tpm(1,a,b,e,z);

         //tp(8)
         if(b == d)
            (*this)(1,i,j) -= A * tpm(1,a,c,e,z);

         //tp(9)
         if(a == d)
            (*this)(1,i,j) += A * tpm(1,b,c,e,z);

      }
   }

   this->symmetrize();

}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(N*(N - 1.0));
   double c = 1.0/(N - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(M - 4.0);
   double b = 1.0/((M - 4.0)*(M - 3.0)*(M - 2.0));
   double c = 1.0/((M - 4.0)*(M - 3.0));

   this->T(a,b,c,tpm);

}

ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int S = 0;S < dpm_p.gnr();++S){

      output << S << "\t" << dpm_p.gdim(S) << "\t" << dpm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < dpm_p.gdim(S);++i)
         for(int j = 0;j < dpm_p.gdim(S);++j){

            output << S << "\t" << i << "\t" << j << "\t|\t" << 
            
               dpm_p.dp2s[S][i][0] << "\t" << dpm_p.dp2s[S][i][1] << "\t" << dpm_p.dp2s[S][i][2] << "\t" << dpm_p.dp2s[S][i][3] << 

               "\t" << dpm_p.dp2s[S][j][0] << "\t" << dpm_p.dp2s[S][j][1] << "\t" << dpm_p.dp2s[S][j][2] << "\t" << dpm_p.dp2s[S][j][3] << "\t" << dpm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * Input from file, with sp indices as row and column indices
 * @param filename Name of the inputfile
 */
void DPM::in_sp(const char *filename){

   ifstream input(filename);

   int S_ab,S_de;
   int a,b,c,d,e,z;
   int S;

   double value;

   int i,j;

   (*this) = 0;

   while(input >> S >> S_ab >> a >> b >> c >> S_de >> d >> e >> z >> value){

      if(S == 0){

         i = s2dp[S][S_ab][a][b][c];
         j = s2dp[S][S_de][d][e][z];

      }
      else{

         i = s2dp[S][0][a][b][c];
         j = s2dp[S][0][d][e][z];

      }

      (*this)(S,i,j) = value;

   }

   this->symmetrize();

}
