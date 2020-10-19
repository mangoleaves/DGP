#include "PolyVectorField.h"

// inputs:
//  mesh: data structure in OpenMesh
// outputs:
//  mesh property "base": complex base of a face
void PolyVectorField::CalcBase(Mesh& mesh)
{
	auto base = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, std::vector<OpenMesh::Vec3d>>(mesh, "base");
	for (auto fh : mesh.faces())
	{
		OpenMesh::SmartEdgeHandle eh = *mesh.fe_begin(fh);
		OpenMesh::Vec3d base1, base2;
		base1 = mesh.calc_edge_vector(eh).normalize();
		base2 = mesh.calc_face_normal(fh).cross(base1).normalize();
		base[fh].push_back(base1);
		base[fh].push_back(base2);

//		std::cout << "face " << fh.idx() << ": base1" << base1 << "\n base2" << base2 << std::endl;
	}
}

// inputs:
//  mesh: data structure in OpenMesh
//  N: the degree of the field
// output:
//  Afull: The resulting left-hand side matrices
void PolyVectorField::Precompute(Mesh& mesh, int N, Eigen::SparseMatrix<std::complex<double>>& Afull)
{
	// Build the sparse matrix, with an energy term for each edge and degree
	// Without constraints on the field, degree is N.
	int rowCounter = 0;
	std::vector<Eigen::Triplet<std::complex<double>>> AfullTriplets;
	auto base = OpenMesh::getProperty<OpenMesh::FaceHandle, std::vector<OpenMesh::Vec3d>>(mesh, "base");

	for (auto eh : mesh.edges())
	{
		if (eh.is_boundary())
		{
			continue;
		}

		OpenMesh::Vec3d e = mesh.calc_edge_vector(eh);
		OpenMesh::Vec2d vef = OpenMesh::Vec2d(e.dot(base[eh.halfedge(0).face()][0]), e.dot(base[eh.halfedge(0).face()][1])).normalize();
		std::complex<double> ef(vef[0], vef[1]);
		OpenMesh::Vec2d veg = OpenMesh::Vec2d(e.dot(base[eh.halfedge(1).face()][0]), e.dot(base[eh.halfedge(1).face()][1])).normalize();
		std::complex<double> eg(veg[0], veg[1]);

		AfullTriplets.push_back(Eigen::Triplet<std::complex<double>>(rowCounter, eh.halfedge(0).face().idx(), std::pow(std::conj(ef), N)));
		AfullTriplets.push_back(Eigen::Triplet<std::complex<double>>(rowCounter++, eh.halfedge(1).face().idx(), -1.0 * std::pow(std::conj(eg), N)));
	}

	Afull.conservativeResize(rowCounter, mesh.n_faces());
	Afull.setFromTriplets(AfullTriplets.begin(), AfullTriplets.end());
}

// inputs:
//  mesh: data structure in OpenMesh
//  N: the degree of the field
//  Afull: The resulting left-hand side matrices
// outputs:
//  polyvector
void PolyVectorField::CalcField(Mesh& mesh, int N, Eigen::SparseMatrix<std::complex<double>>& Afull, Eigen::VectorXcd& polyVectorField)
{
	Eigen::SparseMatrix<std::complex<double>> Lcomplex = Afull.adjoint() * Afull;
	Eigen::SparseMatrix<double> M;
	igl::speye(2 * Lcomplex.outerSize(), 2 * Lcomplex.outerSize(), M);
	Eigen::SparseMatrix<double> L(2 * Lcomplex.outerSize(), 2 * Lcomplex.outerSize());
	
	std::vector<Eigen::Triplet<double>> LTriplets;
	for (int k = 0; k < Lcomplex.outerSize(); k++)
	{
		for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(Lcomplex, k); it; ++it)
		{
			LTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value().real()));
			LTriplets.push_back(Eigen::Triplet<double>(it.row(), Lcomplex.cols() + it.col(), -it.value().imag()));
			LTriplets.push_back(Eigen::Triplet<double>(Lcomplex.rows() + it.row(), it.col(), it.value().imag()));
			LTriplets.push_back(Eigen::Triplet<double>(Lcomplex.rows() + it.row(), Lcomplex.cols() + it.col(), it.value().real()));
		}
	}
	L.setFromTriplets(LTriplets.begin(), LTriplets.end());
	Eigen::MatrixXd U;
	Eigen::VectorXd S;
	igl::eigs(L, M, 1, igl::EIGS_TYPE_SM, U, S);

	polyVectorField = U.block(0, 0, U.rows() / 2, 1).cast<std::complex<double> >().array() * std::complex<double>(1, 0) +
		U.block(U.rows() / 2, 0, U.rows() / 2, 1).cast<std::complex<double> >().array() * std::complex<double>(0, 1);

//	std::cout << "polyVectorField:\n" << polyVectorField << std::endl;
}

void PolyVectorField::ShowField(Mesh& mesh, int N, Eigen::VectorXcd& polyVectorField)
{
	auto base = OpenMesh::getProperty<OpenMesh::FaceHandle, std::vector<OpenMesh::Vec3d>>(mesh, "base");

	// calculate average edge length
	double edgeLenSum = 0;
	for (auto eh : mesh.edges())
	{
		edgeLenSum += mesh.calc_edge_length(eh);
	}
	double avgEdgeLen = edgeLenSum / mesh.n_edges() * 0.1;

	for (auto fh : mesh.faces())
	{
		Mesh::Point centroid = mesh.calc_centroid(fh);
		std::complex<double> coord = polyVectorField(fh.idx());

		std::vector<OpenMesh::Vec3d> vecs;
		vecs.push_back((coord.real() * base[fh][0] + coord.imag() * base[fh][1]).normalize() * avgEdgeLen);
		vecs.push_back((-coord.imag() * base[fh][0] + coord.real() * base[fh][1]).normalize() * avgEdgeLen);
		vecs.push_back((-coord.real() * base[fh][0] + -coord.imag() * base[fh][1]).normalize() * avgEdgeLen);
		vecs.push_back((coord.imag() * base[fh][0] - coord.real() * base[fh][1]).normalize() * avgEdgeLen);

		for (auto it : vecs)
		{
			auto cenh = mesh.add_vertex(centroid);
			auto vech = mesh.add_vertex(centroid + it);
			auto temph = mesh.add_vertex(centroid + it / 2);
			std::vector<OpenMesh::SmartVertexHandle> vhs;
			vhs.push_back(cenh);
			vhs.push_back(vech);
			vhs.push_back(temph);
			mesh.add_face(vhs);
		}
	}
}
