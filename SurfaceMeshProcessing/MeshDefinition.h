#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
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
	// pseudo inverse
	static void PseudoInverse(Eigen::Matrix3d& A, Eigen::Matrix3d& Apinv);
	// Initialize QEM
	static void InitQ(Mesh& mesh);
	// Calculate error
	static void CalcError(Eigen::Matrix4d& eQ, OpenMesh::Vec3d& vbar, double& error);
	// Min heap
	class minHeap
	{
	public:
		struct elem
		{
			double error;
			OpenMesh::SmartEdgeHandle eh;

			elem():error(0), eh(-1) {}
			elem(double e, OpenMesh::SmartEdgeHandle h): error(e), eh(h) {}
		};
		int heap_size;
		std::vector<elem> heap;
		OpenMesh::PropertyManager<OpenMesh::EPropHandleT<int>, int> pos;
	public:
		minHeap(Mesh& mesh) : heap_size(mesh.n_edges()), pos(mesh, "position") { heap.resize(mesh.n_edges() + 1); }

		inline bool less(elem e1, elem e2)
		{
			return e1.error < e2.error;
		}
		
		void make_heap(void);
		double min_value(void);
		elem extract_min(void);
		void set_value(OpenMesh::SmartEdgeHandle eh, double new_value);
		void delete_elem(OpenMesh::SmartEdgeHandle eh);
	private:
		inline int leftc(int i) { return 2 * i; }
		inline int rightc(int i) { return 2 * i + 1; }
		inline int parent(int i) { return i >> 1; }
		void up_move(int idx);
		void down_move(int idx);
	};
	// Simplify
	static void Simplify(Mesh& mesh, int Nv, double minCost);
};
