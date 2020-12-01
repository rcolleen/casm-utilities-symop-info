#ifndef WYCKOFF_HH
#define WYCKOFF_HH

#include <casmutils/definitions.hpp>
#include <vector>
#include <casmutils/sym/cartesian.hpp>


class Subspace
{
public:
    //put this constructor in private
    Subspace(const std::vector<Eigen::Vector3d>& basis_vectors, const Eigen::Vector3d& offset);
    //matrixXd and offset vector
    Subspace(const Eigen::Matrix3d& basis_col_matrix, const Eigen::Vector3d& offset);
    Subspace();//check if we need this, maybe get rid of it
    std::string formula() const; // Sesha suggested we overload the << operator
    Eigen::Vector3d offset()const {return this-> m_offset;};
    Eigen::Matrix3d basis_col_matrix() const {return this-> m_basis_col_matrix;}; 
                                         /* [100     [x
                                             010  dot y  + [offset]
                                             001]     z]*/
//    Subspace operator*(const casmutils::sym::CartOp& lhs);
    
private:
    std::vector<Eigen::Vector3d> orthonormalize_basis(const std::vector<Eigen::Vector3d>& input_vectors);
    Eigen::Vector3d m_offset;
    Eigen::Matrix3d m_basis_col_matrix;
};



#endif
