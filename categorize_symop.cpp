#include "./categorize_symop.hh"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/sym/cartesian.hpp>

#include <functional>
#include <optional>
#include <stdexcept>

//#define tol 1e-6

bool compare_diff_to_prec(Eigen::Vector3d difference, double tol)
{
    int sig_diff = 0;
    for (int i = 0; i < 3; i++)
    {
        if (abs(difference(i)) > tol)
        {
            sig_diff++;
        }
    }
    if (sig_diff != 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool has_translation(const Eigen::Vector3d translation, const Lattice lattice_input, double tol)
{ // check if translation is 0 or integer multiple of lattice vectors-> false,
    // else truei
    Eigen::Matrix3d lattice = lattice_input.column_vector_matrix();
    Eigen::Vector3d lattice_displacement = lattice.inverse() * translation;
    Eigen::Vector3d modification(0.5, 0.5, 0.5);
    Eigen::Vector3d modif_displ = lattice_displacement + modification;
    Eigen::Vector3d trunc_displ;

    for (int i = 0; i < 3; i++)
    {
        trunc_displ(i) = (int)modif_displ(i);
    }
    Eigen::Vector3d difference = lattice_displacement - trunc_displ;
    bool check_trans = compare_diff_to_prec(difference, tol);

    if (check_trans == true)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Eigen::Vector3d project_translation_onto_vectors(const std::vector<Eigen::Vector3d>& eigen_vectors,
                                                 Eigen::Vector3d translation)
{
    // for each eigen vector this call the above project onto vector function and then added these vectors together
    Eigen::Vector3d projected_translation_vector = Eigen::Vector3d::Zero();
    for (const auto& vector : eigen_vectors)
    {
        double projected_component = vector.dot(translation);
        projected_translation_vector += (projected_component * vector);
    }
    return projected_translation_vector;
}

bool almost_equal(double LHS, double RHS, double tol) { return std::abs(LHS - RHS) < tol; }

/**
 * Returns a container of eigenvectors of the given matrix that have
 * positive unit eigenvalues
 */

std::vector<Eigen::Vector3d> eigenvectors_with_positive_unit_eigenvalues(const Eigen::Matrix3d& cart_matrix, double tol)
{
    // some function that evaluate eigne vlaue and counts how many ones
    // Now returns a matrix of the Eigen vectors with eigen values of 1, in column vectors
    // need to edit so that it returns only the real components of the eigenvector
    Eigen::EigenSolver<Eigen::Matrix3d> solver(cart_matrix, true);
    Eigen::Vector3<std::complex<double>> eigenvals = solver.eigenvalues();
    std::vector<Eigen::Vector3d> output_eigen_vectors;
    Eigen::Matrix3<std::complex<double>> eigenvectors = solver.eigenvectors();

    int ct = 0;
    for (int i = 0; i < 3; i++)
    {
        //      std::cout <<eigenvectors.col(i).real()<<std::endl;
        auto eigenval = eigenvals(i);
        //        std::cout<<eigenval.real()<<std::endl;
        if (almost_equal(eigenval.real(), 1, tol))
        {
            ct++;
            //            std::cout<<"EigenVal=1 Found!"<<std::endl;
            Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
            output_eigen_vectors.push_back(real_eigenvector);
        }
    }
    if (ct == 0)
    {
        for (int i = 0; i < 3; i++)
        {
            //        std::cout <<eigenvectors.col(i).real()<<std::endl;
            auto eigenval = eigenvals(i);
            //          std::cout<<eigenval.real()<<std::endl;
            if (almost_equal(eigenval.real(), -1, tol))
            {
                ct++;
                //              std::cout<<"EigenVal=1 Found!"<<std::endl;
                Eigen::Vector3d real_eigenvector = eigenvectors.col(i).real();
                output_eigen_vectors.push_back(real_eigenvector);
            }
        }
    }

    return output_eigen_vectors;
}

Eigen::Matrix4d get_combined_matrix(casmutils::sym::CartOp symop)
{
    
    Eigen::Matrix4d combined_matrix = Eigen::Matrix4d::Zero();

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            combined_matrix(i, j) = symop.matrix(i, j);
        }
        combined_matrix(i, 3) = symop.translation(i);
    }
    combined_matrix(3, 3) = 1;

    return combined_matrix;
}

SYMOP_TYPE check_op_type(const casmutils::sym::CartOp sym_op, const casmutils::xtal::Lattice lattice, double tol)
{ // take in sym_op returns string of op type
    Eigen::Matrix3d matrix = sym_op.matrix;
    int trace = sym_op.matrix.trace();
    std::string type;
    double det = sym_op.matrix.determinant();

    if ((3 - trace < tol))
    {
        return SYMOP_TYPE::IDENTITY;
    }
    if (trace + 3 < tol)
    {
        return SYMOP_TYPE::INVERSION;
    }
    std::vector<Eigen::Vector3d> eigen_vectors =
        eigenvectors_with_positive_unit_eigenvalues(sym_op.matrix, tol);
    //    std::cout<<eigen_vectors<<std::endl;
    //    std::cout << eigen_vectors.cols() << std::endl;
    if (eigen_vectors.size() == 2)
    {
        if (has_translation(project_translation_onto_vectors(eigen_vectors, sym_op.translation), lattice, tol))
        {
            return SYMOP_TYPE::GLIDE;
        }
        else
        {
            return SYMOP_TYPE::MIRROR;
        }
    }
    else if (eigen_vectors.size() == 1)
    {
        if (abs(det - 1) < tol)
        {
            if (has_translation(project_translation_onto_vectors(eigen_vectors, sym_op.translation), lattice,
                                tol))
            {
                return SYMOP_TYPE::SCREW;
            }
            else
            {
                return SYMOP_TYPE::ROTATION;
            }
        }
        else if (abs(det + 1) < tol)
        {
            return SYMOP_TYPE::IMPROPER_ROTATION;
        }
        else
        {
            throw std::runtime_error("Error type 1: Type not idenitified!!!");
        }
    }
    else
    {
        throw std::runtime_error("Error type 2: Type not idenitified!!!");
    }
}

std::optional<Subspace> find_invariant_subspace(casmutils::sym::CartOp symop, casmutils::xtal::Lattice lattice, double tol)
{
    auto type = check_op_type(symop, lattice, tol);
    if (type == SYMOP_TYPE::SCREW || type == SYMOP_TYPE::GLIDE)
    {
        return {};
    }
    auto symop_matrix = get_combined_matrix(symop);
    Eigen::Matrix3d basis_matrix = Eigen::Matrix3d::Zero();
    Eigen::Vector3d offset = Eigen::Vector3d::Zero();
    Eigen::EigenSolver<Eigen::Matrix4d> solver(symop_matrix, true);
    Eigen::Vector4<std::complex<double>> eigenvals = solver.eigenvalues();
    std::vector<Eigen::Vector4d> output_eigen_vectors;
    Eigen::Matrix4<std::complex<double>> eigenvectors = solver.eigenvectors();
    int dimension = 0;
    for (int i = 0; i < 4; i++)
    {
        if ((1 - abs(eigenvals(i).real())) < tol)
        {
            if (1 - abs(eigenvectors.col(i)(3)) < tol)
            {
                for (int k = 0; k < 3; k++)
                {
                    offset(k) = eigenvectors.col(i)(k).real();
                }
            }
            else
            {
                for (int k = 0; k < 3; k++)
                {
                    basis_matrix(k, dimension) = eigenvectors(k, i).real();
                }
                dimension++;
            }
        }
    }
    Subspace invariant_subspace = Subspace(basis_matrix, offset);

    return invariant_subspace;
}

void print_symop_type(SYMOP_TYPE type){
    if (type==SYMOP_TYPE::SCREW){std::cout<<"SCREW";}
    else if (type==SYMOP_TYPE::GLIDE){std::cout<<"GLIDE";}
    else if (type==SYMOP_TYPE::ROTATION){std::cout<<"ROTATION";}
    else if (type==SYMOP_TYPE::MIRROR){std::cout<<"MIRROR";}
    else if (type==SYMOP_TYPE::IMPROPER_ROTATION){std::cout<<"IMPROPER_ROTATION";}
    else if (type==SYMOP_TYPE::INVERSION){std::cout<<"INVERSION";}
    else if (type==SYMOP_TYPE::IDENTITY){std::cout<<"IDENTITY";}
    
    return;
}

