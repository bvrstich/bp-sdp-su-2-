#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"
#include "PPHM.h"

/**
 * @author Brecht Verstichel
 * @date 23-04-2010\n\n
 * This class, PHM, is a class written for spinsymmetrical particle-hole matrices, it inherits all the functions from its mother class
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the ph basis.
 */
class PHM : public BlockMatrix {

    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << phm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << phm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the PHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHM &phm_p);

   public:
      
      //constructor
      PHM(int M,int N);

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //change the numbers in sp mode
      double &operator()(int S,int a,int b,int c,int d);

      //change the numbers in sp mode: read mode
      double operator()(int S,int a,int b,int c,int d) const;

      //geef N terug
      int gN() const;

      //geef N terug
      int gM() const;

      void G(const TPM &);

      void uncouple(const char *filename);

      //trace the first pair of indices of a PPHM object
      void bar(const PPHM &);

      //input PHM from file
      void in_sp(const char *);

   private:

      //!static counter that counts the number of PHM objects running in the program
      static int counter;

      //!static list of dimension [n_ph][2] that takes in a ph index i and returns two sp indices: a = ph2s[i][0] and b = ph2s[i][1]
      static int **ph2s;

      //!static list of dimension [M/2][M/2] that takes two sp indices a,b and returns a ph index i: i = s2ph[a][b]
      static int **s2ph;

      //!list of 6j symbols needed.
      static double **_6j;

      //!number of particles
      int N;

      //!dimension of sp hilbert space
      int M;

};

#endif
