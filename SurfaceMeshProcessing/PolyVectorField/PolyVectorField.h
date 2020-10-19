#pragma once

#include "MeshDefinition.h"

class PolyVectorField
{
public:
	static void CalcBase(Mesh& mesh);
	static void Precompute(Mesh& mesh, int N, Eigen::SparseMatrix<std::complex<double>>& Afull);
	static void CalcField(Mesh& mesh, int N, Eigen::SparseMatrix<std::complex<double>>& Afull, Eigen::VectorXcd& polyVectorField);
	static void ShowField(Mesh& mesh, int N, Eigen::VectorXcd& polyVectorField);
};

