#ifndef VECTOROPS_H
#define VECTOROPS_H

#include <casmutils/definitions.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <vector>

/*struct VectorCompare_f
{
	VectorCompare_f(const Eigen::Vector3d& vector, double prec);
	bool operator()(const Eigen::Vector3d& other) const;
	private:
	Eigen::Vector3d m_vector;
	double m_precision;
};

struct VectorPeriodicCompare_f
{
    VectorPeriodicCompare_f(const Eigen::Vector3d& vector, double prec, const Lattice& unit_cell);
    bool operator()(const Eigen::Vector3d& other) const;
private:
    Eigen::Vector3d m_vector;
    double m_precision;
    Lattice m_lattice;
};

Eigen::Vector3d operator*(const SymOp& transformation, const Eigen::Vector3d& vector);
*/
Eigen::Vector3d projector_operator(Eigen::Vector3d vect_a, Eigen::Vector3d vect_b);
#endif
