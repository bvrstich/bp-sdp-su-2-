#ifndef DPM_H
#define DPM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-02-2010\n\n
 * This class, DPM, is a class written for spinsymmetrical three-particle matrices (name comes from drie-particle matrix). It is written specially for the T_1 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give the relationship between the dp (three-particle) and the sp basis. This matrix falls apart in two blocks: S = 1/2 with degeneracy 2 and S = 3/2 with degeneracy 4. The basis is determined by the spatial orbital quantumnumbers a,b,c , an intermediate spincoupling quantumnumber S_ab = 0 or 1, and the total spin S.
 */
class DPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << dpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << dpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpm_p the DPM you want to print
    */
   friend ostream &operator<<(ostream &output,const DPM &dpm_p);

   public:
      
      //constructor
      DPM(int M,int N);

      //copy constructor
      DPM(const DPM &);

      //destructor
      virtual ~DPM();

      void construct_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int f) const;

      int get_inco(int S,int S_ab,int a,int b,int c,int *i,double *coef) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //generalized T1 map
      void T(double,double,double,const TPM &);

      //maak een DPM van een TPM via de T1 conditie
      void T(const TPM &);

      //maak een DPM van een TPM via de hat functie
      void hat(const TPM &);

      //input DPM from file
      void in_sp(const char *);

   private:

      //!static counter that counts the number of DPM objects running in the program
      static int counter;

      //!static list of dimension [2][dim[i]][4] that takes in a dp index i for block S and returns an intermediate spin: S_ab = dp2s[S][i][0] and three sp indices: a = dp2s[S][i][1], b = dp2s[S][i][1] and c = dp2s[S][i][2]
      static int ***dp2s;

      //!static list of dimension [2][2][M/2][M/2][M/2] that takes a block index S, an intermediate spin-index S_ab and three sp indices a,b and c and returns a dp index i: i = s2dp[S][S_ab][a][b][c]
      static int *****s2dp;

      //!list of 6j symbols needed.
      static double **_6j;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

};

#endif
