#ifndef PPHM_H
#define PPHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 03-05-2010\n\n
 * This class, PPHM, is a class written for spinsymmetrical two-particle-one-hole matrices. It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give the relationship between the pph (two-particle one hole) and the sp basis. This matrix has 2 blocks, one S = 1/2 block with degeneracy 2 and one S = 3/2 block with degeneracy 4.
 */
class PPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << pphm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << pphm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PPHM &pphm_p);

   public:
      
      //constructor
      PPHM(int M,int N);

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      void construct_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int f) const;

      static int get_inco(int S,int S_ab,int a,int b,int c,int &i);

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //maak een PPHM van een TPM via de T2 conditie
      void T(const TPM &);

      //input PPHM from file
      void in_sp(const char *);

   private:

      //!static counter that counts the number of PPHM objects running in the program
      static int counter;

      //!static list of dimension [2][dim[i]][4] that takes in a pph index i and a blockindex for spin, and returns three sp indices: a = pph2s[S][i][1], b = pph2s[S][i][2] and c = pph2s[S][i][3] and an intermediate spin S_ab = pph2s[S][i][0]
      static int ***pph2s;

      //!static list of dimension [2][2][M/2][M/2][M/2] that takes three sp indices a,b and c, a blockindex S for total spin, and an intermediate spinindex S_ab, and returns a pph index i: i = s2pph[S][S_ab][a][b][c]
      static int *****s2pph;

      //!list of 6j symbols needed.
      static double **_6j;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

};

#endif
