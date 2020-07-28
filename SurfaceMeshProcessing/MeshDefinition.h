#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#ifdef _DEBUG
#pragma comment(lib, "OpenMeshCored.lib")
#pragma comment(lib, "OpenMeshToolsd.lib")
#else
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
#endif
struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
};
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;

class MeshTools
{
private:
	typedef Eigen::Triplet<double> T;
	typedef Eigen::SparseMatrix<double> SpMat;

	static bool Parameterization(Mesh& mesh);
	static void CalcLocalCoordinate(Mesh& mesh);
	static void CalcSt(Mesh& mesh);
	static void CalcLt(Mesh& mesh);
	static void CalcCoeMat(Mesh& mesh, SpMat& coeMat);
	static void CalcBVec(Mesh& mesh, Eigen::VectorXd& bx, Eigen::VectorXd& by);
	static double CalcEnergy(Mesh& mesh);
public:
	static bool ReadMesh(Mesh & mesh, const std::string & filename);
	static bool ReadOBJ(Mesh & mesh, const std::string & filename);
	//static bool ReadOFF(Mesh & mesh, const std::string & filename);
	static bool WriteMesh(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static bool WriteOBJ(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static double Area(const Mesh & mesh);
	static double AverageEdgeLength(const Mesh & mesh);
	static bool HasBoundary(const Mesh & mesh);
	static bool HasOneComponent(const Mesh & mesh);
	static int Genus(const Mesh & mesh);
	static void BoundingBox(const Mesh & mesh, Mesh::Point & bmax, Mesh::Point & bmin);
	static void Reassign(const Mesh & mesh1, Mesh & mesh2);
	static bool ParameterizationARAP(Mesh& mesh, Mesh& paraMesh, int maxIter, double minEnergyVar);
};
