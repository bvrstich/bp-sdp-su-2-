#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>

using std::ostream;
using std::ofstream;

#include "BlockMatrix.h"

class SUP;
class PHM;
class DPM;
class PPHM;

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry included, it inherits alle the function from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:
      
      //constructor
      TPM(int M,int N);

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and with spin quantumnumer
      double operator()(int S,int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      void hubbard(double U);

      //Q afbeelding en zijn inverse
      void Q(int option,const TPM &);

      //Q like afbeelding Q(A,B,C,tpm_d)
      void Q(int option,double A,double B,double C,const TPM &);

      //overlapmatrix afbeelding en zijn inverse
      void S(int option,const TPM &);

      void init();

      void set_unit();

      void set_S_2();

      void proj_Tr();

      //de hessiaan afbeelding:
      void H(const TPM &b,const SUP &D);

      //los het stelsel op
      int solve(TPM &b,const SUP &D);

      void collaps(int option,const SUP &);

      void sp_pairing(double );

      void uncouple(const char *);

      //G down afbeelding
      void G(const PHM &);

      //trace one pair of indices of DPM
      void bar(const DPM &);

      //T1 down
      void T(const DPM &);

      //trace last pair of indices of PPHM
      void bar(const PPHM &);

      //T2 down
      void T(const PPHM &);

      //return the spin
      double spin() const;

      //input TPM from file
      void in_sp(const char *);

      void uncouple_ofstream(ofstream &output) const;

      static void constr_overlap(int,int);

      void S_L(int,const TPM &);

   private:

      //!static list of dimension [2][dim[i]][2] that takes in a tp index i and a spinquantumnumber S, and returns two sp indices: a = t2s[S][i][0] and b = t2s[S][i][1]
      static int ***t2s;

      //!static list of dimension [2][M/2][M/2] that takes two sp indices a,b and a spinquantumnumber S, and returns a tp index i: i = s2t[S][a][b]
      static int ***s2t;

      //!list of 6j symbols needed.
      static double **_6j;

      //!static variables needed for the overlapmatrix
      static double S_a,S_b,S_c;

      //!static counter that counts the number of TPM objects running in the program
      static int counter;

      //!nr of particles
      int N;

      //!dimension of sp hilbert space
      int M;

};

#endif
