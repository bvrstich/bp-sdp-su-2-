#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor: Watch out, the matrices themself haven't been allocated yet. 
 * Only the Blockmatrix itself is allocated and the array containing the dimensions, but not initialized.
 * @param nr nr of blocks in the blockmatrix
 */
BlockMatrix::BlockMatrix(int nr){

   this->nr = nr;

   blockmatrix = new Matrix * [nr];

   dim = new int [nr];

   flag = new int [nr];

   degen = new int [nr];

   //init flag:
   for(int i = 0;i < nr;++i)
      flag[i] = 0;

}

/**
 * copy constructor, make sure the input matrix and all the blocks have been allocated and filled before the copying
 * @param blockmat_copy The blockmatrix you want to be copied into the object you are constructing
 */
BlockMatrix::BlockMatrix(const BlockMatrix &blockmat_copy){

   this->nr = blockmat_copy.nr;

   blockmatrix = new Matrix * [nr];

   dim = new int [nr];

   flag = new int [nr];

   degen = new int [nr];

   for(int i = 0;i < nr;++i){

      flag[i] = 1;

      degen[i] = blockmat_copy.gdeg(i);

      dim[i] = blockmat_copy.gdim(i);

      blockmatrix[i] = new Matrix(blockmat_copy[i]);

   }

}

/**
 * Destructor
 */
BlockMatrix::~BlockMatrix(){

   //if the memory was allocated during its lifetime, deallocate the Matrix memory
   for(int i = 0;i < nr;++i)
      if(flag[i] == 1)
         delete blockmatrix[i];

   delete [] blockmatrix;

   delete [] flag;
   delete [] dim;
   delete [] degen;

}

/**
 * Function that allocates the memory of the blockmatrix block with dimension dim sets and the degeneracy of the matrix
 * @param block The index of the block that will be allocated
 * @param dim the dimension of the particular block to be allocated
 * @param degeneracy the degeneracy of block "block"
 */
void BlockMatrix::setMatrixDim(int block,int dim,int degeneracy){

   this->flag[block] = 1;

   this->dim[block] = dim;

   this->degen[block] = degeneracy;

   blockmatrix[block] = new Matrix(dim);

}

/**
 * [] overloaded, will return a reference to the block block in the blockmatrix.
 * @param block index of the block to be returned
 * @return A reference to the Matrix object located on blockmatrix[block].
 */
Matrix &BlockMatrix::operator[](int block){
   
   return *blockmatrix[block];

}

/**
 * [] overloaded, will return a reference to the block block in the blockmatrix: const version
 * @param block index of the block to be returned
 * @return A reference to the Matrix object located on blockmatrix[block].
 */
const Matrix &BlockMatrix::operator[](int block) const{
   
   return *blockmatrix[block];

}


/**
 * overload the equality operator: Make sure the blocks in both matrices have been allocated to the same dimensions and have the same degeneracy!
 * @param blockmat_copy The matrix you want to be copied into this
 */
BlockMatrix &BlockMatrix::operator=(const BlockMatrix &blockmat_copy){

   for(int i = 0;i < nr;++i)
      *blockmatrix[i] = blockmat_copy[i];

   return *this;

}

/**
 * Make all the numbers in your blockmatrix equal to the number a, e.g. usefull for initialization (BlockMatrix M = 0)
 * @param a the number
 */
BlockMatrix &BlockMatrix::operator=(double a){

   for(int i = 0;i < nr;++i)
      *blockmatrix[i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param blockmat_pl The matrix you want to add to this
 */
BlockMatrix &BlockMatrix::operator+=(const BlockMatrix &blockmat_pl){

   for(int i = 0;i < nr;++i)
      *blockmatrix[i] += blockmat_pl[i];

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param blockmat_pl The matrix you want to deduct from this
 */
BlockMatrix &BlockMatrix::operator-=(const BlockMatrix &blockmat_pl){

   for(int i = 0;i < nr;++i)
      *blockmatrix[i] -= blockmat_pl[i];

   return *this;

}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param blockmat_pl the BlockMatrix to be multiplied by alpha and added to this
 */
BlockMatrix &BlockMatrix::daxpy(double alpha,const BlockMatrix &blockmat_pl){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->daxpy(alpha,blockmat_pl[i]);

   return *this;

}
/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
BlockMatrix &BlockMatrix::operator/=(double c){

   for(int i = 0;i < nr;++i)
      *blockmatrix[i] /= c;

   return *this;

}

/**
 * write access to your blockmatrix, change the number in block "block", on row i and column j
 * @param block The index of the block you want to access
 * @param i row number
 * @param j column number
 * @return the entry on place block,i,j
 */
double &BlockMatrix::operator()(int block,int i,int j){

   return (*blockmatrix[block])(i,j);

}

/**
 * read access to your blockmatrix, read the number in block "block" on row i and column j
 * @param block The index of the block you want to access
 * @param i row number
 * @param j column number
 * @return the entry on place block,i,j
 */
double BlockMatrix::operator()(int block,int i,int j) const {

   return (*blockmatrix[block])(i,j);

}

/**
 * @return the nr of blocks
 */
int BlockMatrix::gnr() const{

   return nr;

}

/**
 * @return the dimension of the Matrix on the block with index i
 */
int BlockMatrix::gdim(int i) const{

   return dim[i];

}

/**
 * @return the degeneracy of the block with index i
 */
int BlockMatrix::gdeg(int i) const{

   return degen[i];

}

/**
 * @return the trace of the matrix, each block matrix is weighed with its degeneracy.
 */
double BlockMatrix::trace() const{

   double ward = 0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockmatrix[i]->trace();

   return ward;

}

/**
 * @return inproduct of (*this) blockmatrix with blockmatrix_in, defined as Tr (A B)
 * @param blockmatrix_in input matrix
 */
double BlockMatrix::ddot(const BlockMatrix &blockmatrix_in) const{

   double ward = 0.0;

   for(int i = 0;i < nr;++i)
      ward += degen[i]*blockmatrix[i]->ddot(blockmatrix_in[i]);

   return ward;

}

/**
 * Invert positive semidefinite symmetric blockmatrix which is stored in (*this), original matrix (*this) is destroyed
 */
void BlockMatrix::invert(){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->invert();

}

/**
 * Scale the blockmatrix (*this) with parameter alpha
 * @param alpha scalefactor
 */
void BlockMatrix::dscal(double alpha){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->dscal(alpha);

}

/**
 * Fill the matrix with random numbers.
 */
void BlockMatrix::fill_Random(){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->fill_Random();

}

/**
 * Take the square root out of the positive semidefinite blockmatrix, destroys original blocks, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void BlockMatrix::sqrt(int option){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->sqrt(option);

}

/**
 * Multiply symmetric blockmatrix object left en right with symmetric blockmatrix map to 
 * form another symmetric blockmatrix and put it in (*this): this = map*object*map
 * @param map BlockMatrix that will be multiplied to the left en to the right of matrix object
 * @param object central BlockMatrix
 */
void BlockMatrix::L_map(const BlockMatrix &map,const BlockMatrix &object){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->L_map(map[i],object[i]);
   
}

/**
 * BlockMatrix product of two general blockmatrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
BlockMatrix &BlockMatrix::mprod(const BlockMatrix &A, const BlockMatrix &B){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->mprod(A[i],B[i]);

   return *this;

}

/**
 * Copy upper triangle into lower triangle.
 */
void BlockMatrix::symmetrize(){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->symmetrize();

}

ostream &operator<<(ostream &output,const BlockMatrix &blockmatrix_p){

   for(int i = 0;i < blockmatrix_p.nr;++i){

      output << i << "\t" << blockmatrix_p.dim[i] << "\t" << blockmatrix_p.degen[i] << endl;
      output << endl;

      output << *blockmatrix_p.blockmatrix[i] << endl;

   }

   return output;

}

/**
 * Output the Blockmatrix to a file with name filename
 * @param filename string with filename
 */
void BlockMatrix::out(const char *filename) const{

   ofstream output(filename);

   output.precision(10);

   output << *this;

}

/**
 * Seperate SUP into two SUP's, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void BlockMatrix::sep_pm(BlockMatrix &m,BlockMatrix &n){

   for(int i = 0;i < nr;++i)
      blockmatrix[i]->sep_pm(m[i],n[i]);

}
