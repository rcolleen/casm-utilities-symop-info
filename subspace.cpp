//#include <casm/external/Eigen/Dense>
#include <casmutils/definitions.hpp>
//#include <casmutils/sym/cartesian.hpp>
//#include <casmutils/xtal/lattice.hpp>
//#include <casmutils/xtal/structure_tools.hpp>
#include <vector>
//#include <casmutil/xtal/structure.hpp>
//#include <filesystem>
#include "./subspace.hh"
#include "./vectorfunctions.hpp"
#include <cmath>


std::vector<Eigen::Vector3d> Subspace::orthonormalize_basis(const std::vector<Eigen::Vector3d>& input_vectors)
{
    int dim=input_vectors.size();
    std::vector<Eigen::Vector3d> output_vectors;
    for(int i=0;i<3; i++){
        if (i>dim-1){
            output_vectors.push_back(Eigen::Vector3d::Zero());}
        else
        {
             Eigen::Vector3d output_vec=input_vectors[i].normalized();
             if (i>0){
                 Eigen::Vector3d temp_proj=projector_operator(output_vec, output_vectors[0]);
                 if (i>1){
                     temp_proj=temp_proj+projector_operator(output_vec, output_vectors[1]);
                 }
                output_vec=output_vec-temp_proj;
             }
             
             output_vectors.push_back(output_vec.normalized());
        }
    }
             
    return output_vectors;
}

Subspace::Subspace(const std::vector<Eigen::Vector3d>& input_vectors, const Eigen::Vector3d& offset){
    std::vector<Eigen::Vector3d> basis_vectors=orthonormalize_basis( input_vectors);

    m_offset = offset -projector_operator(offset, basis_vectors[0]) - projector_operator(offset,basis_vectors[1]) - projector_operator(offset, basis_vectors[2]);
    for (int i = 0; i < 3; ++i)
    {
        m_basis_col_matrix(i, 0) = basis_vectors[0](i);
        m_basis_col_matrix(i, 1) = basis_vectors[1](i);
        m_basis_col_matrix(i, 2) = basis_vectors[2](i);
    }

    //This may still be a different subspace though order if given axis, or plane could be define in different orders
    //DONE: make sure planes are defined by orthoganol vectors?
    //      ake sure offset is not in subspace
}

Subspace::Subspace(const Eigen::Matrix3d& input_basis_vectors, const Eigen::Vector3d& input_offset)
    :Subspace({input_basis_vectors.col(0), input_basis_vectors.col(1), input_basis_vectors.col(2)}, input_offset) {}

Subspace::Subspace(): m_offset(Eigen::Vector3d::Zero()) , m_basis_col_matrix(Eigen::Matrix3d::Zero()){}

/*Subspace Subspace::operator*(const casmutils::sym::CartOp& lhs)
{
    Eigen::Matrix3d basis_vectors = lhs.matrix * m_basis_col_matrix;
    Eigen::Vector3d offset_vector = lhs * m_offset;
    Subspace product_wycoff_position(basis_vectors, offset_vector);
    return product_wycoff_position;
}
*/
std::string Subspace::formula() const
{
    double tol = 1e-5;
    std::string formula = "{";
    std::vector<std::string> xyz = {"x", "y", "z"};

    for (int i = 0; i < 3; ++i)
    {
        std::string temp_string;
        for (int j = 0; j < 3; ++j)
        {
            if (std::abs(this->m_basis_col_matrix(i,j)) > tol)
            {
                std::string sign = this->m_basis_col_matrix(i,j) > 0 ? "+" : "";
                temp_string += sign + std::to_string(this->m_basis_col_matrix(i,j)) + xyz[j];
            }
        }
        if (temp_string.size() == 0)
        {
            temp_string = "0";
        }
        formula += temp_string + ", ";
    }

    formula.pop_back();
    formula.pop_back();
    formula += "}";
    return formula;
}

