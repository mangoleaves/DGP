#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshDefinition.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <cctype>

bool MeshTools::ReadMesh(Mesh & mesh, const std::string & filename)
{
	std::string path(".");
	auto slash = filename.find_last_of('/');
	auto backslash = filename.find_last_of('\\');
	if (slash != std::string::npos && backslash != std::string::npos)
	{
		slash = slash > backslash ? slash : backslash;
		path = filename.substr(0, slash);
	}
	else if (backslash != std::string::npos)
	{
		slash = backslash;
		path = filename.substr(0, slash);
	}
	else if (slash != std::string::npos)
	{
		path = filename.substr(0, slash);
	}

	std::string basename = filename.substr(slash + 1);
	auto point = basename.find_last_of('.');
	std::string ext("obj");
	if (point != std::string::npos)
	{
		ext = basename.substr(point + 1);
		basename = basename.substr(0, point);
	}
	std::transform(ext.begin(), ext.end(), ext.begin(), std::tolower);
	if (ext == "obj")
	{
		return ReadOBJ(mesh, path + "/" + basename + "." + ext);
	}
	else
	{
		//std::cout << "Error: the file extension " << ext << " is not supported." << std::endl;
		return OpenMesh::IO::read_mesh(mesh, filename);
	}
}

bool MeshTools::ReadOBJ(Mesh & mesh, const std::string & filename)
{
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return false;
	}
	std::string line;
	std::string keyWrd;
	Mesh::Point::value_type x, y, z;
	std::vector<Mesh::VertexHandle> vertexHandles;
	std::stringstream stream, lineData, tmp;

	// pass 1: read vertices
	while (ifs && !ifs.eof())
	{
		std::getline(ifs, line);
		if (ifs.bad())
		{
			std::cerr << "  Warning! Could not read file properly!\n";
			return false;
		}

		// Trim Both leading and trailing spaces
		auto start = line.find_first_not_of(" \t\r\n");
		auto end = line.find_last_not_of(" \t\r\n");
		if ((std::string::npos == start) || (std::string::npos == end))
			line = "";
		else
			line = line.substr(start, end - start + 1);

		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}
		stream.str(line);
		stream.clear();
		stream >> keyWrd;
		if (keyWrd == "v")
		{
			stream >> x >> y >> z;
			if (!stream.fail())
			{
				vertexHandles.push_back(mesh.add_vertex(Mesh::Point(x, y, z)));
			}
		}
	}

	// reset stream for second pass
	ifs.clear();
	ifs.seekg(0, std::ios::beg);

	int nCurrentPositions = 0;

	// pass 2: read faces
	while (ifs && !ifs.eof())
	{
		std::getline(ifs, line);
		if (ifs.bad())
		{
			std::cerr << "  Warning! Could not read file properly!\n";
			return false;
		}

		// Trim Both leading and trailing spaces
		auto start = line.find_first_not_of(" \t\r\n");
		auto end = line.find_last_not_of(" \t\r\n");
		line = ((std::string::npos == start) || (std::string::npos == end)) ? "" : line.substr(start, end - start + 1);

		// comment
		if (line.size() == 0 || line[0] == '#' || isspace(line[0]))
		{
			continue;
		}

		stream.str(line);
		stream.clear();
		stream >> keyWrd;

		// track current number of parsed vertex attributes,
		// to allow for OBJs negative indices
		if (keyWrd == "v")
		{
			++nCurrentPositions;
		}

		// faces
		else if (keyWrd == "f")
		{
			int value;

			// read full line after detecting a face
			std::string faceLine;
			std::getline(stream, faceLine);
			lineData.str(faceLine);
			lineData.clear();

			Mesh::FaceHandle fh;
			std::vector<Mesh::VertexHandle> faceVertices;

			// work on the line until nothing left to read
			while (!lineData.eof())
			{
				// read one block from the line ( vertex/texCoord/normal )
				std::string vertex;
				lineData >> vertex;

				//get the component (vertex/texCoord/normal)
				size_t found = vertex.find("/");

				// parts are seperated by '/' So if no '/' found its the last component
				if (found != std::string::npos)
				{
					// read the index value
					tmp.str(vertex.substr(0, found));
					tmp.clear();

					// If we get an empty string this property is undefined in the file
					if (vertex.substr(0, found).empty())
					{
						continue;
					}

					// Read current value
					tmp >> value;
				}
				else
				{
					// last component of the vertex, read it.
					tmp.str(vertex);
					tmp.clear();
					tmp >> value;

					// Nothing to read here ( garbage at end of line )
					if (tmp.fail())
					{
						continue;
					}
				}

				if (value < 0) {
					// Calculation of index :
					// -1 is the last vertex in the list
					// As obj counts from 1 and not zero add +1
					value = nCurrentPositions + value + 1;
				}
				// Obj counts from 1 and not zero .. array counts from zero therefore -1
				faceVertices.push_back(Mesh::VertexHandle(value - 1));
			}

			// note that add_face can possibly triangulate the faces, which is why we have to
			// store the current number of faces first
			size_t n_faces = mesh.n_faces();
			auto endIter = faceVertices.end();
			for (auto iter = faceVertices.begin(); iter != endIter; ++iter)
			{
				endIter = std::remove(iter + 1, endIter, *(iter));
			}
			faceVertices.erase(endIter, faceVertices.end());

			//A minimum of three vertices are required.
			if (faceVertices.size() > 2)
			{
				fh = mesh.add_face(faceVertices);
			}
			std::vector<Mesh::FaceHandle> newfaces;
			for (size_t i = 0; i < mesh.n_faces() - n_faces; ++i)
			{
				newfaces.push_back(Mesh::FaceHandle(int(n_faces + i)));
			}
		}
	}
	ifs.close();
	return true;
}

bool MeshTools::WriteMesh(const Mesh & mesh, const std::string & filename, const std::streamsize & precision)
{
	std::string path(".");
	auto slash = filename.find_last_of('/');
	auto backslash = filename.find_last_of('\\');
	if (slash != std::string::npos && backslash != std::string::npos)
	{
		slash = slash > backslash ? slash : backslash;
		path = filename.substr(0, slash);
	}
	else if (backslash != std::string::npos)
	{
		slash = backslash;
		path = filename.substr(0, slash);
	}
	else if (slash != std::string::npos)
	{
		path = filename.substr(0, slash);
	}

	std::string basename = filename.substr(slash + 1);
	auto point = basename.find_last_of('.');
	std::string ext("obj");
	if (point != std::string::npos)
	{
		ext = basename.substr(point + 1);
		basename = basename.substr(0, point);
	}
	std::transform(ext.begin(), ext.end(), ext.begin(), std::tolower);
	if (ext == "obj")
	{
		return WriteOBJ(mesh, path + "/" + basename + "." + ext, precision);
	}
	else
	{
		//std::cout << "Error: the file extension " << ext << " is not supported." << std::endl;
		return OpenMesh::IO::write_mesh(mesh, filename, OpenMesh::IO::Options::Default, precision);
	}
}

bool MeshTools::WriteOBJ(const Mesh & mesh, const std::string & filename, const std::streamsize & precision)
{
	std::ofstream ofs(filename);
	if (!ofs.is_open())
	{
		std::cerr << "Error: cannot open file " << filename << std::endl;
		return false;
	}
	ofs.precision(precision);
	for (const auto& vh : mesh.vertices())
	{
		ofs << "v " << mesh.point(vh) << std::endl;
	}
	for (const auto& fh : mesh.faces())
	{
		ofs << "f";
		for (const auto& fvh : mesh.fv_range(fh))
		{
			ofs << " " << fvh.idx() + 1;
		}
		ofs << std::endl;
	}
	ofs.close();
	return true;
}

double MeshTools::Area(const Mesh & mesh)
{
	double area = 0.0;
	for (const auto & f : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(f);
		const auto & p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto & p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto & p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		area += 0.5 * ((p1 - p0) % (p2 - p0)).norm();
	}
	return area;
}

double MeshTools::AverageEdgeLength(const Mesh & mesh)
{
	if (mesh.edges_empty()) return 0.0;
	double l = 0.0;
	for (const auto & eh : mesh.edges())
	{
		l += mesh.calc_edge_length(eh);
	}
	l /= mesh.n_edges();
	return l;
}

bool MeshTools::HasBoundary(const Mesh & mesh)
{
	if (mesh.halfedges_empty()) return false;
	for (const auto & heh : mesh.halfedges())
	{
		if (mesh.is_boundary(heh))
		{
			return true;
		}
	}
	return false;
}

bool MeshTools::HasOneComponent(const Mesh & mesh)
{
	if (mesh.faces_empty()) return false;
	std::vector<int> visit(mesh.n_faces(), 0);
	visit[0] = 1;
	std::queue<int> q;
	q.push(0);
	while (!q.empty())
	{
		int fid = q.front();
		q.pop();
		for (const auto & ffh : mesh.ff_range(mesh.face_handle(fid)))
		{
			if (!visit[ffh.idx()])
			{
				visit[ffh.idx()] = 1;
				q.push(ffh.idx());
			}
		}
	}
	for (const auto & v : visit)
	{
		if (!v)
		{
			return false;
		}
	}
	return true;
}

int MeshTools::Genus(const Mesh & mesh)
{
	if (HasBoundary(mesh) || !HasOneComponent(mesh)) return -1;
	return 1 - ((int)mesh.n_vertices() + (int)mesh.n_faces() - (int)mesh.n_edges()) / 2;
}

void MeshTools::BoundingBox(const Mesh & mesh, Mesh::Point & bmax, Mesh::Point & bmin)
{
	if (mesh.vertices_empty()) return;
	bmax = bmin = mesh.point(*mesh.vertices_begin());
	for (const auto vh : mesh.vertices())
	{
		bmax.maximize(mesh.point(vh));
		bmin.minimize(mesh.point(vh));
	}
}

// This function should be used after collapse or split, since the indices of
//   edges may be changed after a local modification. By using this function,
//   the indices of edges are the same as a mesh loaded from obj files.
void MeshTools::Reassign(const Mesh & mesh1, Mesh & mesh2)
{
	mesh2.clear();
	for (const auto & vh : mesh1.vertices())
	{
		mesh2.add_vertex(mesh1.point(vh));
	}
	for (const auto & fh : mesh1.faces())
	{
		std::vector<Mesh::VertexHandle> vhs;
		for (const auto & fvh : mesh1.fv_range(fh))
		{
			vhs.push_back(mesh2.vertex_handle(fvh.idx()));
		}
		mesh2.add_face(vhs);
	}
}

void MeshTools::PseudoInverse(Eigen::Matrix3d& A, Eigen::Matrix3d& Apinv)
{
	Eigen::BDCSVD<Eigen::Matrix3d> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Index rank = svd.rank();
	Apinv = svd.matrixV().leftCols(rank) *
		svd.singularValues().head(rank).asDiagonal().inverse() *
		svd.matrixU().leftCols(rank).adjoint();
}

void MeshTools::InitQ(Mesh& mesh)
{
	auto faceQ = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Eigen::Matrix4d>(mesh, "faceQ");
	auto vertexQ = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, Eigen::Matrix4d>(mesh, "vertexQ");
	auto edgeQ = OpenMesh::getOrMakeProperty<OpenMesh::EdgeHandle, Eigen::Matrix4d>(mesh, "edgeQ");

	// Init faceQ
	for (auto fh : mesh.faces())
	{
		OpenMesh::Vec3d n = mesh.calc_face_normal(fh);
		double di = n.dot(mesh.point(*mesh.fv_begin(fh)));
		Eigen::Vector4d nbar(n[0], n[1], n[2], -di);
		faceQ[fh] = nbar * nbar.transpose();
	}
	// Init vertexQ
	for (auto vh : mesh.vertices())
	{
		vertexQ[vh].setZero();
		for (auto vfIter = mesh.vf_begin(vh); vfIter.is_valid(); vfIter++)
		{
			vertexQ[vh] += faceQ[*vfIter];
		}
	}
	// Init edgeQ
	for (auto eh : mesh.edges())
	{
		edgeQ[eh] = vertexQ[eh.v0()] + vertexQ[eh.v1()];
	}
}

void MeshTools::CalcError(Eigen::Matrix4d& eQ, OpenMesh::Vec3d& vbar, double& error)
{
	Eigen::Matrix3d A, Apinv;
	Eigen::Vector3d b;
	Eigen::Vector4d v;

	A = eQ.block(0, 0, 3, 3);
	b = -eQ.block(0, 3, 3, 1);
	PseudoInverse(A, Apinv);
	v << Apinv * b, 1;

	error = v.transpose() * eQ * v;
	vbar[0] = v[0];
	vbar[1] = v[1];
	vbar[2] = v[2];
}

void MeshTools::Simplify(Mesh& mesh, int Nv, double minCost)
{
	// 初始化所有面、顶点、边的QEM矩阵
	InitQ(mesh);

	auto faceQ = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Eigen::Matrix4d>(mesh, "faceQ");
	auto vertexQ = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, Eigen::Matrix4d>(mesh, "vertexQ");
	auto edgeQ = OpenMesh::getOrMakeProperty<OpenMesh::EdgeHandle, Eigen::Matrix4d>(mesh, "edgeQ");
	auto vbar = OpenMesh::getOrMakeProperty<OpenMesh::EdgeHandle, OpenMesh::Vec3d>(mesh, "vbar");
	// 计算边上的Error，并建堆
	minHeap mh(mesh);
	for (auto eh : mesh.edges())
	{
		double error = 0;
		CalcError(edgeQ[eh], vbar[eh], error);
		mh.heap[eh.idx() + 1] = minHeap::elem(error, eh);
	}
	mh.make_heap();
	// collapse循环
	int n = mesh.n_vertices();
	while (n > Nv && mh.min_value() < minCost)
	{
		// 取得最小元素
		auto min_elem = mh.extract_min();
		auto heh = min_elem.eh.h0();
		auto v1 = heh.to();
		auto vb = vbar[min_elem.eh];

		// 从堆中删去将要被删去的三条边
		if (!heh.is_boundary())
		{
			mh.delete_elem(heh.next().next().edge());
		}
		if (!heh.opp().is_boundary())
		{
			mh.delete_elem(heh.opp().next().edge());
		}
		// collapse
		mesh.collapse(heh);
		n--;
		// 更改新顶点位置
		mesh.set_point(v1, vb);
		// 更新新顶点一环邻面的Q
		for (auto vfIter = mesh.vf_begin(v1); vfIter.is_valid(); vfIter++)
		{
			OpenMesh::Vec3d n = mesh.calc_face_normal(*vfIter);
			double di = n.dot(mesh.point(*mesh.fv_begin(*vfIter)));
			Eigen::Vector4d nbar(n[0], n[1], n[2], -di);
			faceQ[*vfIter] = nbar * nbar.transpose();
		}
		// 更新相关顶点的Q
		vertexQ[v1].setZero();
		for (auto vfIter = mesh.vf_begin(v1); vfIter.is_valid(); vfIter++)
		{
			vertexQ[v1] += faceQ[*vfIter];
		}
		for (auto vvIter = mesh.vv_begin(v1); vvIter.is_valid(); vvIter++)
		{
			vertexQ[*vvIter].setZero();
			for (auto vfIter = mesh.vf_begin(*vvIter); vfIter.is_valid(); vfIter++)
			{
				vertexQ[*vvIter] += faceQ[*vfIter];
			}
		}
		// 更新相关边的Q，并计算Error，更新其在堆中的值和位置
		for (auto vvIter = mesh.vv_begin(v1); vvIter.is_valid(); vvIter++)
		{
			for (auto veIter = mesh.ve_begin(*vvIter); veIter.is_valid(); veIter++)
			{
				edgeQ[*veIter] = vertexQ[(*veIter).v0()] + vertexQ[(*veIter).v1()];
				double error = 0;
				CalcError(edgeQ[*veIter], vbar[*veIter], error);
				mh.set_value(*veIter, error);
			}
		}
	}
}

void MeshTools::minHeap::make_heap(void)
{
	// 已有n个元素位于 heap[1,...,n]中，初始化为堆
	for (int i = 1; i <= heap_size; i++)
	{
		pos[heap[i].eh] = i;
	}
	for (int i = heap_size >> 1; i > 0; i--)
	{
		down_move(i);
	}
}

double MeshTools::minHeap::min_value(void)
{
	return heap[1].error;
}

MeshTools::minHeap::elem MeshTools::minHeap::extract_min(void)
{
	elem res = heap[1];
	pos[heap[1].eh] = -1;
	heap[1] = heap[heap_size];
	heap_size -= 1;
	down_move(1);
	return res;
}

void MeshTools::minHeap::set_value(OpenMesh::SmartEdgeHandle eh, double new_value)
{
	double old_value = heap[pos[eh]].error;
	heap[pos[eh]].error = new_value;
	if (new_value < old_value)
	{
		up_move(pos[eh]);
	}
	else
	{
		down_move(pos[eh]);
	}
}

void MeshTools::minHeap::delete_elem(OpenMesh::SmartEdgeHandle eh)
{
	set_value(eh, -INFINITY);
	extract_min();
}

void MeshTools::minHeap::up_move(int idx)
{
	// 将指定元素向根节点方向调整，使其满足二叉堆性质
	heap[0] = heap[idx];
	int p = parent(idx);
	while (p > 0 && less(heap[0], heap[p]))
	{
		heap[idx] = heap[p];
		pos[heap[idx].eh] = idx;
		idx = p;
		p = parent(idx);
	}
	heap[idx] = heap[0];
	pos[heap[idx].eh] = idx;
}

void MeshTools::minHeap::down_move(int idx)
{
	// 将指定元素向叶节点方向调整，使其满足二叉堆性质
	heap[0] = heap[idx];
	while (idx <= (heap_size >> 1))
	{
		int minIdx = idx;
		int lc = leftc(idx);
		int rc = rightc(idx);
		if (less(heap[lc], heap[0]))
		{
			if (rc <= heap_size && less(heap[rc], heap[lc]))
			{
				minIdx = rc;
			}
			else
			{
				minIdx = lc;
			}
		}
		else if (rc <= heap_size && less(heap[rc], heap[0]))
		{
			minIdx = rc;
		}

		if (minIdx != idx)
		{
			heap[idx] = heap[minIdx];
			pos[heap[idx].eh] = idx;
			idx = minIdx;
		}
		else
		{
			break;
		}
	}
	heap[idx] = heap[0];
	pos[heap[idx].eh] = idx;
}
