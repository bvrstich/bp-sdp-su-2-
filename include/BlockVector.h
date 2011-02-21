#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;

#include "lapack.h"
#include "Vector.h"

template<class BlockMatrixType>
class BlockVector;

/**
 * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
 * ifstream object and type:\n\n
 * object << blockvector << endl;\n\n
 * For output onto the screen type: \n\n
 * cout << blockvector << endl;\n\n
 * @param output The stream to which you are writing (e.g. cout)
 * @param blockvector_p de BlockVector you want to print
 */
template<class BlockMatrixType>
ostream &operator<<(ostream &output,const BlockVector<BlockMatrixType> &blockvector_p);

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This is a class written for blockvectors. It will contain the eigenvalues of the TPM, etc. Matrices. It is a template class,
 * corresponding to the different BlockVectorType's that can be put in, it will automatically get the right dimension and degeneracy.
 * It is a double pointer to a Vector class and uses the functions of vector for its computations.
 */
template<class BlockMatrixType>
class BlockVector{

   public:

      //construct with as input a BlockMatrixType
      BlockVector(BlockMatrixType& );

      //copy constructor
      BlockVector(const BlockVector<BlockMatrixType> &);

      //destructor
      virtual ~BlockVector();

      //overload equality operator
      BlockVector &operator=(const BlockVector<BlockMatrixType> &);

      BlockVector &operator=(double );

      //overload += operator
      BlockVector &operator+=(const BlockVector<BlockMatrixType> &);

      //overload -= operator
      BlockVector &operator-=(const BlockVector<BlockMatrixType> &);

      BlockVector &daxpy(double alpha,const BlockVector<BlockMatrixType> &);

      BlockVector &operator/=(double );

      Vector &operator[](int block);

      //const version
      const Vector &operator[](int block) const;

      //easy to change the numbers
      double &operator()(int block,int index);

      //easy to access the numbers
      double operator()(int block,int index) const;

      void diagonalize(BlockMatrixType &);

      int gnr() const;

      int gdim(int) const;

      int gdeg(int) const;

      double sum() const;

      double log_product() const;

      double ddot(const BlockVector<BlockMatrixType> &) const;

      void dscal(double alpha);

      double min() const;
      
      double max() const;

      double centerpot(double ) const;

   private:

      //!double pointer to Vector, the blockvector
      Vector **blockvector;

      //!nr of blocks in the blockvector
      int nr;

      //!dimension of the blocks
      int *dim;

      //!degeneracy of the blocks
      int *degen;

};

/**
 * Construct and initialize the BlockVector object by diagonalizing a BlockMatrixType object
 * @param MT the matrix of BlockMatrixType that has to be diagonalized
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType>::BlockVector(BlockMatrixType &MT){

   //allocate
   this->nr = MT.gnr();

   blockvector = new Vector * [nr];

   dim = new int [nr];

   degen = new int [nr];

   for(int i = 0;i < nr;++i){

      blockvector[i] = new Vector(MT[i]);

      dim[i] = blockvector[i]->gn();

      degen[i] = MT.gdeg(i);

   }

}

/**
 * copy constructor 
 * @param vec_copy The blockvector you want to be copied into the object you are constructing, make sure that it is an allocated and filled blockvector!
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType>::BlockVector(const BlockVector<BlockMatrixType> &vec_copy){

   this->nr = vec_copy.gnr();

   blockvector = new Vector * [nr];

   dim = new int [nr];

   degen = new int [nr];

   for(int i = 0;i < nr;++i){

      blockvector[i] = new Vector(vec_copy[i]);

      dim[i] = blockvector[i]->gn();

      degen[i] = vec_copy.gdeg(i);

   }

}

/**
 * Destructor
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType>::~BlockVector(){

   for(int i = 0;i < nr;++i)
      delete blockvector[i];

   delete [] blockvector;

   delete [] dim;

   delete [] degen;

}

/**
 * overload the equality operator, make sure dimension and degeneracy of the blocks are the same.
 * @param blockvector_copy The blockvector you want to be copied into this
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::operator=(const BlockVector<BlockMatrixType> &blockvector_copy){

   for(int i = 0;i < nr;++i)
      *blockvector[i] = blockvector_copy[i];

   return *this;

}

/**
 * Make all the number in your blockvector equal to the number a, e.g. usefull for initialization (BlockVector M = 0)
 * @param a the number
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::operator=(double a){

   for(int i = 0;i < nr;++i)
      *blockvector[i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param blockvector_pl The blockvector you want to add to this
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::operator+=(const BlockVector<BlockMatrixType> &blockvector_pl){

   for(int i = 0;i < nr;++i)
      *blockvector[i] += blockvector_pl[i];

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param blockvector_pl The blockvector you want to deduct from this
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::operator-=(const BlockVector<BlockMatrixType> &blockvector_pl){

   for(int i = 0;i < nr;++i)
      *blockvector[i] -= blockvector_pl[i];

   return *this;

}

/**
 * add the blockvector blockvector_pl times the constant alpha to this
 * @param alpha the constant to multiply the blockvector_pl with
 * @param blockvector_pl the BlockVector to be multiplied by alpha and added to this
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::daxpy(double alpha,const BlockVector<BlockMatrixType> &blockvector_pl){

   for(int i = 0;i < nr;++i)
      blockvector[i]->daxpy(alpha,blockvector_pl[i]);

   return *this;

}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your blockvector through
 */
template<class BlockMatrixType>
BlockVector<BlockMatrixType> &BlockVector<BlockMatrixType>::operator/=(double c){

   for(int i = 0;i < nr;++i)
      *blockvector[i] /= c;

   return *this;

}

/**
 * Get the Vector located in block i
 * @param i block number
 * @return the Vector corresponding to this index
 */
template<class BlockMatrixType>
Vector &BlockVector<BlockMatrixType>::operator[](int i){

   return *blockvector[i];

}

/**
 * Get the Vector located in block i: the const version
 * @param i block number
 * @return the Vector corresponding to this index
 */
template<class BlockMatrixType>
const Vector &BlockVector<BlockMatrixType>::operator[](int i) const{

   return *blockvector[i];

}

/**
 * Write access to the number in block "block" and on index i.
 */
template<class BlockMatrixType>
double &BlockVector<BlockMatrixType>::operator()(int block,int i){

   return (*blockvector[block])[i];

}

/**
 * Read access to the number in block "block" and on index i.
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::operator()(int block,int i) const {

   return (*blockvector[block])[i];

}

/**
 * @return the nr of blocks in the blockVector
 */
template<class BlockMatrixType>
int BlockVector<BlockMatrixType>::gnr() const{

   return nr;

}

/**
 * @return the dimension of block "block" in the blockVector
 * @param block index of the block u want the dimension of.
 */
template<class BlockMatrixType>
int BlockVector<BlockMatrixType>::gdim(int block) const{

   return dim[block];

}

/**
 * @return the degeneracy of the block with index i in the BlockVector
 * @param i the index of the block
 */
template<class BlockMatrixType>
int BlockVector<BlockMatrixType>::gdeg(int i) const{

   return degen[i];

}

/**
 * @return the sum of all the elements in the blockvector, weighing the blocks with their degeneracy.
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::sum() const{

   double ward = 0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockvector[i]->sum();

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the blockvector (so the sum of all the logarithms), again weighed with their degeneracies
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::log_product() const{

   double ward = 0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockvector[i]->log_product();

   return ward;

}

/**
 * @return inproduct of (*this) blockvector with blockvector_i
 * @param blockvector_i input blockvector
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::ddot(const BlockVector &blockvector_i) const{

   double ward = 0.0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockvector[i]->ddot(blockvector_i[i]);

   return ward;

}

/**
 * Scale the blockvector (*this) with parameter alpha
 * @param alpha scalefactor
 */
template<class BlockMatrixType>
void BlockVector<BlockMatrixType>::dscal(double alpha){

   for(int i = 0;i < nr;++i)
      blockvector[i]->dscal(alpha);

}

template<class BlockMatrixType>
ostream &operator<<(ostream &output,const BlockVector<BlockMatrixType> &blockvector_p){

   for(int i = 0;i < blockvector_p.gnr();++i){

      std::cout << std::endl;
      std::cout << i << "\t" << blockvector_p.gdim(i) << "\t" << blockvector_p.gdeg(i) << std::endl;
      std::cout << std::endl;

      std::cout << blockvector_p[i] << std::endl;

   }

   return output;

}

/**
 * @return the minimal element present in this BlockVector object.
 * watch out, only works when BlockVector is filled with the eigenvalues of a diagonalized Matrix object
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::min() const{

   double ward = blockvector[0]->min();

   for(int i = 1;i < nr;++i)
      if(ward > blockvector[i]->min())
         ward = blockvector[i]->min();

   return ward;

}

/**
 * @return the maximal element present in this BlockVector object.
 * watch out, only works when BlockVector is filled with the eigenvalues of a diagonalized Matrix object
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::max() const{

   double ward = blockvector[0]->max();

   for(int i = 1;i < nr;++i)
      if(ward < blockvector[i]->max())
         ward = blockvector[i]->max();

   return ward;

}

/**
 * A function needed in the calcalation of the distance from the center, used the Vector::centerpot function 
 * but weighed with the degeneracies of the different blocks.
 * @return the result of the function
 */
template<class BlockMatrixType>
double BlockVector<BlockMatrixType>::centerpot(double alpha) const{

   double ward = 0.0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockvector[i]->centerpot(alpha);

   return ward;

}

/**
 * Diagonalize the BlockMatrix when the vector has allready been allocated
 * @param MT The BlockMatrixType you want to diagonalize
 */
template<class BlockMatrixType>
void BlockVector<BlockMatrixType>::diagonalize(BlockMatrixType &MT){

   for(int i = 0;i < nr;++i)
      blockvector[i]->diagonalize(MT[i]);

}

#endif
