#ifndef CAT_SYMOP
#define CAT_SYMOP

#include <casmutils/definitions.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/symmetry.hpp>
#include <casmutils/sym/cartesian.hpp>
//#include "../avdv-factor-group/include.hpp"
#include "./subspace.hh"
#include <cmath>
#include <complex>
#include <iostream>
#include <optional>
#include <vector>

typedef casmutils::sym::CartOp SymOp;
typedef casmutils::xtal::Lattice Lattice;
bool compare_diff_to_prec(Eigen::Vector3d difference, double tol);
bool has_translation(const Eigen::Vector3d translation, const Eigen::Matrix3d lattice, double tol);
Eigen::Vector3d project_translation_onto_vectors(const std::vector<Eigen::Vector3d>& eigen_vectors,
                                                 Eigen::Vector3d translation);
bool almost_equal(double LHS, double RHS, double tol);
std::vector<Eigen::Vector3d> eigenvectors_with_positive_unit_eigenvalues(const Eigen::Matrix3d& cart_matrix,
                                                                         double tol);
enum class SYMOP_TYPE
{
    IDENTITY,
    INVERSION,
    GLIDE,
    SCREW,
    MIRROR,
    ROTATION,
    IMPROPER_ROTATION
};
Eigen::Matrix4d get_combined_matrix(casmutils::sym::CartOp symop);
SYMOP_TYPE check_op_type(const casmutils::sym::CartOp sym_op, const casmutils::xtal::Lattice lattice, double tol);
std::optional<Subspace> find_invariant_subspace(casmutils::sym::CartOp symop, casmutils::xtal::Lattice lattice, double tol);

void print_symop_type(SYMOP_TYPE type);
#endif
