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

Mesh::Point MeshTools::CircumCenter(Mesh& mesh, OpenMesh::FaceHandle faceHandle)
{
	// 获取面的三个顶点
	std::vector<Mesh::Point> X;
	for (auto fvIter = mesh.fv_begin(faceHandle); fvIter.is_valid(); fvIter++)
	{
		X.push_back(mesh.point(*fvIter));
	}
	// 计算面的法向量
	auto n = mesh.calc_face_normal(faceHandle);
	// 构造线性方程组，计算外心
	Eigen::Matrix3d A;
	Eigen::Vector3d b;
	A << 2 * (X[1].data()[0] - X[0].data()[0]), 2 * (X[1].data()[1] - X[0].data()[1]), 2 * (X[1].data()[2] - X[0].data()[2]),
		2 * (X[2].data()[0] - X[0].data()[0]), 2 * (X[2].data()[1] - X[0].data()[1]), 2 * (X[2].data()[2] - X[0].data()[2]),
		n.data()[0], n.data()[1], n.data()[2];
	b << pow(X[1].data()[0], 2) + pow(X[1].data()[1], 2) + pow(X[1].data()[2], 2) - pow(X[0].data()[0], 2) - pow(X[0].data()[1], 2) - pow(X[0].data()[2], 2),
		pow(X[2].data()[0], 2) + pow(X[2].data()[1], 2) + pow(X[2].data()[2], 2) - pow(X[0].data()[0], 2) - pow(X[0].data()[1], 2) - pow(X[0].data()[2], 2),
		n.data()[0] * X[0].data()[0] + n.data()[1] * X[0].data()[1] + n.data()[2] * X[0].data()[2];

	Eigen::Vector3d cc = A.colPivHouseholderQr().solve(b);
	return Mesh::Point(cc.x(), cc.y(), cc.z());
	
}

void MeshTools::LocalAveragingRegion(Mesh& mesh)
{
	try
	{
		auto vertexLAR = OpenMesh::getProperty<OpenMesh::VertexHandle, double>(mesh, "vertexLAR");
		for (auto vHandle : mesh.vertices())
		{
			vertexLAR[vHandle] = 0.0;
		}

		for (auto faceHandle : mesh.faces())
		{
			// 判断是否是钝角三角形
			bool isObtuseAngle = false;
			OpenMesh::VertexHandle obtuseVertexHandle;
			for (auto fhIter = mesh.fh_begin(faceHandle); fhIter.is_valid(); fhIter++)
			{
				auto angle = mesh.calc_sector_angle(*fhIter);
				if (angle > M_PI / 2)
				{
					isObtuseAngle = true;
					obtuseVertexHandle = mesh.to_vertex_handle(*fhIter);
					break;
				}
			}
			// 计算面积
			if (isObtuseAngle)
			{
				double faceArea = mesh.calc_face_area(faceHandle);
				for (auto fvIter = mesh.fv_begin(faceHandle); fvIter.is_valid(); fvIter++)
				{
					if (*fvIter == obtuseVertexHandle)
					{
						vertexLAR[*fvIter] += faceArea / 2;
					}
					else
					{
						vertexLAR[*fvIter] += faceArea / 4;
					}
				}
			}
			else
			{
				auto cc = CircumCenter(mesh, faceHandle);
				for (auto fhIter = mesh.fh_begin(faceHandle); fhIter.is_valid(); fhIter++)
				{
					auto edgeMidpoint = mesh.calc_edge_midpoint(*fhIter);
					auto edgeLength = mesh.calc_edge_length(*fhIter);
					double partArea = 0.5 * edgeLength * (edgeMidpoint - cc).norm();
					vertexLAR[mesh.to_vertex_handle(*fhIter)] += 0.5 * partArea;
					vertexLAR[mesh.from_vertex_handle(*fhIter)] += 0.5 * partArea;
				}
			}
		}
	}
	catch (const std::exception& x)
	{
		std::cout << x.what() << std::endl;
		return;
	}
}

void MeshTools::MeanCurvature(Mesh& mesh)
{
	try
	{
		auto meanCurvature = OpenMesh::getProperty<OpenMesh::VertexHandle, OpenMesh::Vec3d>(mesh, "meanCurvature");
		auto vertexLAR = OpenMesh::getProperty<OpenMesh::VertexHandle, double>(mesh, "vertexLAR");
		for (auto vHandle : mesh.vertices())
		{
			meanCurvature[vHandle] = OpenMesh::Vec3d(0, 0, 0);
			for (auto vvIter = mesh.vv_begin(vHandle); vvIter.is_valid(); vvIter++)
			{
				auto firstHalfedgeHandle = mesh.find_halfedge(vHandle, *vvIter);
				auto nextHalfedgeHandle = firstHalfedgeHandle.next();

				double firstAngle = mesh.calc_sector_angle(firstHalfedgeHandle);
				double nextAngle = mesh.calc_sector_angle(nextHalfedgeHandle);

				auto firstVertexHandle = mesh.to_vertex_handle(firstHalfedgeHandle);
				auto nextVertexHandle = mesh.to_vertex_handle(nextHalfedgeHandle);

				meanCurvature[vHandle] += OpenMesh::Vec3d(mesh.point(nextVertexHandle) - mesh.point(vHandle)) / tan(firstAngle);
				meanCurvature[vHandle] += OpenMesh::Vec3d(mesh.point(firstVertexHandle) - mesh.point(vHandle)) / tan(nextAngle);
			}
			meanCurvature[vHandle] /= 4 * vertexLAR[vHandle];
		}
	}
	catch (const std::exception& x)
	{
		std::cout << x.what() << std::endl;
		return;
	}
}

void MeshTools::AbsoluteMeanCurvature(Mesh& mesh)
{
	try
	{
		auto meanCurvature = OpenMesh::getProperty<OpenMesh::VertexHandle, OpenMesh::Vec3d>(mesh, "meanCurvature");
		auto absoluteMeanCurvature = OpenMesh::getProperty<OpenMesh::VertexHandle, double>(mesh, "absoluteMeanCurvature");
		for (auto vHandle : mesh.vertices())
		{
			absoluteMeanCurvature[vHandle] = meanCurvature[vHandle].norm();
		}
	}
	catch (const std::exception& x)
	{
		std::cout << x.what() << std::endl;
		return;
	}
}

void MeshTools::GaussianCurvature(Mesh& mesh)
{
	try
	{
		auto vertexLAR = OpenMesh::getProperty<OpenMesh::VertexHandle, double>(mesh, "vertexLAR");
		auto gaussianCurvature = OpenMesh::getProperty<OpenMesh::VertexHandle, double>(mesh, "gaussianCurvature");

		for (auto vHandle : mesh.vertices())
		{
			gaussianCurvature[vHandle] = 2 * M_PI;
			for (auto vvIter = mesh.vv_iter(vHandle); vvIter.is_valid(); vvIter++)
			{
				auto outHalfedgeHandle = mesh.find_halfedge(vHandle, *vvIter);
				auto inHalfedgeHandle = outHalfedgeHandle.next().next();
				double theta = mesh.calc_sector_angle(inHalfedgeHandle);
				gaussianCurvature[vHandle] -= theta;
			}
			gaussianCurvature[vHandle] /= vertexLAR[vHandle];
		}
	}
	catch (const std::exception& x)
	{
		std::cout << x.what() << std::endl;
	}
}

double MeshTools::BoxMuller(double mu, double sigma)
{
	double u1 = (double)rand() / RAND_MAX;
	double u2 = (double)rand() / RAND_MAX;
	double z = sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
	return mu + z * sigma;
}

void MeshTools::GaussianNoise(Mesh& mesh, double proportion)
{
	if (mesh.vertices_empty())
	{
		return;
	}

	srand(time(NULL));

	double sumEdgeLength = 0.0;
	double avgEdgeLength = 0.0;

	for (auto eh : mesh.edges())
	{
		sumEdgeLength += mesh.calc_edge_length(eh);
	}
	avgEdgeLength = sumEdgeLength / mesh.n_edges();

	double mu = 0.0;
	double sigma = avgEdgeLength * proportion;

	for (auto vh : mesh.vertices())
	{
		auto vp = mesh.point(vh);
		double dx = BoxMuller(mu, sigma);
		double dy = BoxMuller(mu, sigma);
		double dz = BoxMuller(mu, sigma);
		vp += Mesh::Point(dx, dy, dz);
		mesh.set_point(vh, vp);
	}
}

void MeshTools::NormalUpdating(Mesh& mesh, double sigmaS, int maxIter)
{
	try
	{
		// 初始化
		auto updatedNormal = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Mesh::Normal>(mesh, "updatedNormal");
		auto faceCentroid = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh);
		auto faceArea = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, double>(mesh);
		for (auto fh : mesh.faces())
		{
			updatedNormal[fh] = mesh.calc_face_normal(fh);
			faceCentroid[fh] = mesh.calc_face_centroid(fh);
			faceArea[fh] = mesh.calc_face_area(fh);
		}
		// 计算sigmaC，使用相邻面的中心间的平均距离
		double sigmaC = 0.0;
		for (auto fh : mesh.faces())
		{
			auto centroid = mesh.calc_face_centroid(fh);
			for (auto ffIter = mesh.ff_begin(fh); ffIter.is_valid(); ffIter++)
			{
				sigmaC += (centroid - mesh.calc_face_centroid(*ffIter)).norm();
			}
		}
		sigmaC /= mesh.n_faces() * 3;
		// 为简便表示和计算，作处理 sigma = -2 * sigma^2
		sigmaS = -2 * pow(sigmaS, 2);
		sigmaC = -2 * pow(sigmaC, 2);
		// 迭代更新面法向量
		for (int i = 0; i < maxIter; i++)
		{
			for (auto fh : mesh.faces())
			{
				// 先求和式部分，累加到un，并将weight * Ws * Wc累加到K
				// 最后更新法向量为 un / K
				Mesh::Normal un(0.0, 0.0, 0.0);
				double K = 0.0;
				auto centroid = faceCentroid[fh];
				for (auto ffIter = mesh.ff_begin(fh); ffIter.is_valid(); ffIter++)
				{
					double weight = faceArea[*ffIter];
					double Ws = exp(pow((updatedNormal[fh] - updatedNormal[*ffIter]).norm(), 2) / sigmaS);
					double Wc = exp(pow((centroid - faceCentroid[*ffIter]).norm(), 2) / sigmaC);
					K += weight * Ws * Wc;
					un += weight * Ws * Wc * mesh.calc_face_normal(*ffIter);
				}
				updatedNormal[fh] = un / K;
				updatedNormal[fh] /= updatedNormal[fh].norm();
			}
		}
	}
	catch (const std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}

void MeshTools::VertexUpdating(Mesh& mesh, int maxIter)
{
	try
	{
		// 获得updatedNormal属性
		auto updatedNormal = OpenMesh::getProperty<OpenMesh::FaceHandle, Mesh::Normal>(mesh, "updatedNormal");
		// 计算原体积
		double volume = 0.0;
		for (auto fh : mesh.faces())
		{
			auto faceNormal = mesh.calc_face_normal(fh);
			auto faceArea = mesh.calc_face_area(fh);
			auto fvIter = mesh.fv_begin(fh);
			volume += mesh.point(*fvIter).dot(faceNormal) * faceArea;
		}
		// 迭代更新顶点位置
		for (int i = 0; i < maxIter; i++)
		{
			for (auto vh : mesh.vertices())
			{
				auto oldPoint = mesh.point(vh);
				auto increment = Mesh::Point(0.0, 0.0, 0.0);
				double N = 0.0;
				for (auto vfIter = mesh.vf_begin(vh); vfIter.is_valid(); vfIter++)
				{
					auto centroid = mesh.calc_face_centroid(*vfIter);
					increment += updatedNormal[*vfIter] * (updatedNormal[*vfIter].dot(centroid - oldPoint));
					N += 1.0;
				}
				increment /= N;
				mesh.set_point(vh, oldPoint + increment);
			}
			// 保持体积
			double newVolume = 0.0;
			for (auto fh : mesh.faces())
			{
				auto faceNormal = mesh.calc_face_normal(fh);
				auto faceArea = mesh.calc_face_area(fh);
				auto fvIter = mesh.fv_begin(fh);
				newVolume += mesh.point(*fvIter).dot(faceNormal) * faceArea;
			}
			double proportion = pow(volume / newVolume, 1.0 / 3.0);

			for (auto vh : mesh.vertices())
			{
				mesh.set_point(vh, mesh.point(vh) * proportion);
			}
		}
	}
	catch (const std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}

void MeshTools::Denoise(Mesh& mesh, double sigmaS, int normalMaxIter, int vertexMaxIter)
{
	if (mesh.vertices_empty())
	{
		return;
	}
	NormalUpdating(mesh, sigmaS, normalMaxIter);
	VertexUpdating(mesh, vertexMaxIter);
}

double MeshTools::Error(Mesh& originalMesh, Mesh& processedMesh)
{
	double error = 0.0;
	double sumArea = 0.0;
	for (auto vh : originalMesh.vertices())
	{
		auto nvp = processedMesh.point(processedMesh.vertex_handle(vh.idx()));
		for (auto vfIter = originalMesh.vf_begin(vh); vfIter.is_valid(); vfIter++)
		{
			auto fv = originalMesh.point(*originalMesh.fv_begin(*vfIter));
			double area = originalMesh.calc_face_area(*vfIter);
			double dist = pow((nvp-fv).dot(originalMesh.calc_face_normal(*vfIter)), 2);
			sumArea += area;
			error += area * dist;
		}
	}
	return sqrt(error / sumArea);
}
