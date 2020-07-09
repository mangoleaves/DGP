#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
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
	struct ShortestPath {
		std::vector<OpenMesh::VertexHandle>* path;
		double distance;

		ShortestPath() :path(NULL), distance(0) {}
		ShortestPath(std::vector<OpenMesh::VertexHandle>* p, double dis) :path(p), distance(dis) {}
	};
	typedef std::vector<std::pair<OpenMesh::VertexHandle, OpenMesh::VertexHandle>> edges;

	static bool ReadMesh(Mesh& mesh, const std::string& filename);
	static bool ReadOBJ(Mesh& mesh, const std::string& filename);
	//static bool ReadOFF(Mesh & mesh, const std::string & filename);
	static bool WriteMesh(const Mesh& mesh, const std::string& filename, const std::streamsize& precision = 6);
	static bool WriteOBJ(const Mesh& mesh, const std::string& filename, const std::streamsize& precision = 6);
	static double Area(const Mesh& mesh);
	static double AverageEdgeLength(const Mesh& mesh);
	static bool HasBoundary(const Mesh& mesh);
	static bool HasOneComponent(const Mesh& mesh);
	static int Genus(const Mesh& mesh);
	static void BoundingBox(const Mesh& mesh, Mesh::Point& bmax, Mesh::Point& bmin);
	static void Reassign(const Mesh& mesh1, Mesh& mesh2);
	static ShortestPath FindShortestPath(Mesh& mesh, int beginIndex, int endIndex);
	static edges FindMST(Mesh& mesh, std::vector<int> vertexIdxs);
};
