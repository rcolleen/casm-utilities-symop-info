#include <casmutils/xtal/lattice.hpp>
#include <iostream>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/sym/cartesian.hpp>
#include <CLI/CLI.hpp>
#include <filesystem>
#include <vector>
#include "./subspace.hh"
#include "vectorfunctions.hpp"
#include "categorize_symop.hh"

casmutils::sym::CartOp read_symop(const std::filesystem::path& input_path)
{
    std::ifstream mat_stream(input_path);

    Eigen::Matrix3d mat;
    Eigen::Vector3d translation;
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            mat_stream>>mat(i,j);
        }
    }

    for (int i=0; i<3; ++i)
    {
         mat_stream>>translation(i);
    }

    mat_stream.close();
    casmutils::sym::CartOp symop(mat, translation, false);
    
    return symop;
}

casmutils::xtal::Lattice read_lattice_file(const std::filesystem::path& input_path)
{
    std::ifstream mat_stream(input_path);

    Eigen::Matrix3d mat;
    Eigen::Vector3d translation;
    for(int i=0; i<3; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            mat_stream>>mat(i,j);
        }
    }

    mat_stream.close();
    return casmutils::xtal::Lattice::from_column_vector_matrix(mat);
    
}
int main(int argc, char* argv[])
{
    namespace cu=casmutils;

    // Parse CLI
    CLI::App app;

    // Path to input structure
    cu::fs::path input_symop_path;
    cu::fs::path input_lattice_path;
    cu::fs::path output_struc_path;

    app.add_option("-s,--symop_file", input_symop_path, "Path to input symop, file with column vector matrix followedby translation.")->required()->check(CLI::ExistingFile);
    app.add_option("-l,--lattice_file", input_lattice_path, "Path to input lattice, file with column vector matrix.")->required()->check(CLI::ExistingFile);
  //  app.add_option("-o,--output",output_struc_path,"Path to write transformed structure to.");

    CLI11_PARSE(app, argc, argv);

    // Read input symop from file
    cu::sym::CartOp input_symop = read_symop(input_symop_path);

    // Read in lattice
    cu::xtal::Lattice lattice = read_lattice_file(input_lattice_path);
    // print inputsi
    
    std::cout<<"Input SymOp:\n"<< input_symop.matrix<<"\n"<<input_symop.translation<<std::endl;
  //  std::cout<<"This is a :\n"<< check_op_type(input_symop, lattice, 1e-6)<<std::endl;

    std::cout<<"This Symmetry Operation is a ";
    SYMOP_TYPE type = check_op_type(input_symop, lattice, 1e-6);
    print_symop_type(type);
    std::cout<<"."<<std::endl;
    
    std::cout<<"The Invariant Subspace of this operations is: \n";
    if (type==SYMOP_TYPE::SCREW || type==SYMOP_TYPE::GLIDE){
        std::cout<<"Glide and Screw Symmetry Operations have no invariant subspace."<<std::endl;}
    else{
        auto subspace=find_invariant_subspace(input_symop, lattice, 1e-6);
        std::cout<<subspace.value().formula()<<std::endl;
    }


        

    return 0;
}
