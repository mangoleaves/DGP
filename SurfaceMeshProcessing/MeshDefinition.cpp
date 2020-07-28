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

// Tutte+Barycentric Mapping，参数化坐标储存在mesh的UV属性中
bool MeshTools::Parameterization(Mesh& mesh)
{
	// 初始化
	auto isBoundary = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, bool>(mesh);
	std::vector<int> boundaryIndices;
	OpenMesh::VertexHandle firstBoundaryVh;

	for (auto vh : mesh.vertices())
	{
		isBoundary[vh] = false;
	}
	// 寻找第一个边界点
	for (auto vh : mesh.vertices())
	{
		if (mesh.is_boundary(vh))
		{
			firstBoundaryVh = vh;
			boundaryIndices.push_back(vh.idx());
			isBoundary[vh] = true;
			break;
		}
	}
	if (boundaryIndices.empty())
	{
		std::cout << "网格无边界。" << std::endl;
		return false;
	}
	// 依次寻找其余边界点
	bool isContinue;
	OpenMesh::VertexHandle nextBoundaryVh = firstBoundaryVh;
	do
	{
		isContinue = false;
		for (auto vvIter = mesh.vv_begin(nextBoundaryVh); vvIter.is_valid(); vvIter++)
		{
			if (mesh.is_boundary(*vvIter) && !isBoundary[*vvIter])
			{
				isBoundary[*vvIter] = true;
				boundaryIndices.push_back((*vvIter).idx());
				nextBoundaryVh = *vvIter;
				isContinue = true;
				break;
			}
		}
	} while (isContinue);
	// 寻找是否还有其他边界点，以确定该网格和圆盘同胚
	for (auto vh : mesh.vertices())
	{
		if (mesh.is_boundary(vh) && !isBoundary[vh])
		{
			std::cout << "网格有多个边界。" << std::endl;
			return false;
		}
	}
	// 将边界顶点均匀映射到圆上
	auto paraCoordinate = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, Eigen::Vector2d>(mesh, "paraCoordinate");

	double r = 1000.0;
	double theta = 0.0;
	double thetaIncrement = 2 * M_PI / (double)boundaryIndices.size();

	for (auto idx : boundaryIndices)
	{
		auto vh = mesh.vertex_handle(idx);
		paraCoordinate[vh][0] = r * cos(theta);
		paraCoordinate[vh][1] = r * sin(theta);
		theta += thetaIncrement;
	}
	// 为内部顶点确定连续的行索引，以便构造线性方程组
	auto vertexRow = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, int>(mesh);
	int row = 0;

	for (auto vh : mesh.vertices())
	{
		if (!isBoundary[vh])
		{
			vertexRow[vh] = row;
			row++;
		}
	}
	// 构造稀疏线性方程组
	int nInnerVertex = mesh.n_vertices() - boundaryIndices.size();
	SpMat A(nInnerVertex, nInnerVertex);
	std::vector<T> coefficients;
	Eigen::VectorXd bu(nInnerVertex);
	Eigen::VectorXd bv(nInnerVertex);

	for (int i = 0; i < nInnerVertex; i++)
	{
		bu[i] = 0;
		bv[i] = 0;
	}
	for (auto vh : mesh.vertices())
	{
		if (!isBoundary[vh])
		{
			double N = 0;
			for (auto vvIter = mesh.vv_begin(vh); vvIter.is_valid(); vvIter++)
			{
				N += 1.0;
				if (isBoundary[*vvIter])
				{
					bu[vertexRow[vh]] += paraCoordinate[*vvIter][0];
					bv[vertexRow[vh]] += paraCoordinate[*vvIter][1];
				}
				else
				{
					coefficients.push_back(T(vertexRow[vh], vertexRow[*vvIter], -1));
				}
			}
			coefficients.push_back(T(vertexRow[vh], vertexRow[vh], N));
		}
	}
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	A.makeCompressed();
	Eigen::SimplicialCholesky<SpMat> chol(A);
	Eigen::VectorXd u = chol.solve(bu);
	Eigen::VectorXd v = chol.solve(bv);
	// 置内部顶点的参数坐标
	for (auto vh : mesh.vertices())
	{
		if (!isBoundary[vh])
		{
			paraCoordinate[vh][0] = u[vertexRow[vh]];
			paraCoordinate[vh][1] = v[vertexRow[vh]];
		}
	}
	return true;
}

void MeshTools::CalcLocalCoordinate(Mesh& mesh)
{
	auto localCoordinate = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, std::map<int, Eigen::Vector2d>>(mesh, "localCoordinate");

	for (auto fh : mesh.faces())
	{
		auto fvIter = mesh.fv_begin(fh);
		auto Vi = *fvIter;
		auto Vj = *(++fvIter);
		auto Vk = *(++fvIter);
		// 计算X轴、Y轴方向的单位向量
		OpenMesh::Vec3d vecX = (mesh.point(Vj) - mesh.point(Vi)).normalize();
		OpenMesh::Vec3d vecY = mesh.calc_face_normal(fh).cross(mesh.point(Vj) - mesh.point(Vi)).normalize();
		// 将3个点投影到X、Y轴上
		localCoordinate[fh][Vi.idx()] = Eigen::Vector2d(0.0, 0.0);
		localCoordinate[fh][Vj.idx()] = Eigen::Vector2d(vecX.dot(mesh.point(Vj) - mesh.point(Vi)), 0.0);
		localCoordinate[fh][Vk.idx()] = Eigen::Vector2d(vecX.dot(mesh.point(Vk) - mesh.point(Vi)), vecY.dot(mesh.point(Vk) - mesh.point(Vi)));
	}
}

// S_t = Sum_{i=0}^{2} cot(theta^i)(u^i - u^{i+1})(x^i - x^{i+1})^T
void MeshTools::CalcSt(Mesh& mesh)
{
	auto faceSt = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(mesh, "faceSt");
	auto localCoordinate = OpenMesh::getProperty<OpenMesh::FaceHandle, std::map<int, Eigen::Vector2d>>(mesh, "localCoordinate");
	auto paraCoordinate = OpenMesh::getProperty<OpenMesh::VertexHandle, Eigen::Vector2d>(mesh, "paraCoordinate");

	for (auto fh : mesh.faces())
	{
		faceSt[fh] << 0, 0, 0, 0;
		for (auto fhIter = mesh.fh_begin(fh); fhIter.is_valid(); fhIter++)
		{
			auto Vi = (*fhIter).from();
			auto Vj = (*fhIter).to();
			double theta = mesh.calc_sector_angle((*fhIter).next());
			faceSt[fh] += (paraCoordinate[Vi] - paraCoordinate[Vj]) *
				(localCoordinate[fh][Vi.idx()] - localCoordinate[fh][Vj.idx()]).transpose() / tan(theta);
		}
	}
}

// 对S_t signed SVD分解，S_t = U * S * V^T，则L_t = U * V^T
void MeshTools::CalcLt(Mesh& mesh)
{
	auto faceSt = OpenMesh::getProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(mesh, "faceSt");
	auto faceLt = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(mesh, "faceLt");

	for (auto fh : mesh.faces())
	{
		auto SVD = faceSt[fh].jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
		auto U = SVD.matrixU();
		auto V = SVD.matrixV();
		Eigen::Matrix2d Lt = U * V.transpose();
		if (Lt.determinant() < 0)
		{
			V(0, 1) = -V(0, 1);
			V(1, 1) = -V(1, 1);
			Lt = U * V.transpose();
		}
		faceLt[fh] = Lt;
	}
}

void MeshTools::CalcCoeMat(Mesh& mesh, SpMat& coeMat)
{
	std::vector<T> coefficients;
	for (auto vh : mesh.vertices())
	{
		int idx_i = vh.idx();
		double coe_ii = 0.0;
		for (auto vvIter = mesh.vv_begin(vh); vvIter.is_valid(); vvIter++)
		{
			int idx_j = (*vvIter).idx();
			auto he_ij = mesh.find_halfedge(vh, *vvIter);
			auto he_ji = he_ij.opp();
			double coe = 0.0;
			if (!he_ij.is_boundary())
			{
				double theta_ij = mesh.calc_sector_angle(he_ij.next());
				coe += 1.0 / tan(theta_ij);
			}
			if (!he_ji.is_boundary())
			{
				double theta_ji = mesh.calc_sector_angle(he_ji.next());
				coe += 1.0 / tan(theta_ji);
			}
			coe_ii += coe;
			coefficients.push_back(T(idx_i, idx_j, -coe));
		}
		coefficients.push_back(T(idx_i, idx_i, coe_ii));
	}
	coeMat.resize(mesh.n_vertices(), mesh.n_vertices());
	coeMat.setFromTriplets(coefficients.begin(), coefficients.end());
	coeMat.makeCompressed();
}

void MeshTools::CalcBVec(Mesh& mesh, Eigen::VectorXd& bx, Eigen::VectorXd& by)
{
	auto localCoordinate = OpenMesh::getProperty<OpenMesh::FaceHandle, std::map<int, Eigen::Vector2d>>(mesh, "localCoordinate");
	auto faceLt = OpenMesh::getProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(mesh, "faceLt");

	for (auto vh : mesh.vertices())
	{
		int idx_i = vh.idx();
		Eigen::Vector2d b_i(0.0, 0.0);
		for (auto vvIter = mesh.vv_begin(vh); vvIter.is_valid(); vvIter++)
		{
			int idx_j = (*vvIter).idx();
			auto he_ij = mesh.find_halfedge(vh, *vvIter);
			auto he_ji = he_ij.opp();
			if (!he_ij.is_boundary())
			{
				auto x_i = localCoordinate[he_ij.face()][idx_i];
				auto x_j = localCoordinate[he_ij.face()][idx_j];
				double theta_ij = mesh.calc_sector_angle(he_ij.next());
				auto Lt_ij = faceLt[he_ij.face()];
				b_i += Lt_ij * (x_i - x_j) / tan(theta_ij);
			}
			if (!he_ji.is_boundary())
			{
				auto x_i = localCoordinate[he_ji.face()][idx_i];
				auto x_j = localCoordinate[he_ji.face()][idx_j];
				double theta_ji = mesh.calc_sector_angle(he_ji.next());
				auto Lt_ji = faceLt[he_ji.face()];
				b_i += Lt_ji * (x_i - x_j) / tan(theta_ji);
			}
		}
		bx[idx_i] = b_i[0];
		by[idx_i] = b_i[1];
	}
}

double MeshTools::CalcEnergy(Mesh& mesh)
{
	auto localCoordinate = OpenMesh::getProperty<OpenMesh::FaceHandle, std::map<int, Eigen::Vector2d>>(mesh, "localCoordinate");
	auto paraCoordinate = OpenMesh::getProperty<OpenMesh::VertexHandle, Eigen::Vector2d>(mesh, "paraCoordinate");
	auto faceLt = OpenMesh::getProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(mesh, "faceLt");

	double E = 0.0;
	for (auto heh : mesh.halfedges())
	{
		if (!heh.is_boundary())
		{
			auto vi = heh.from();
			auto vj = heh.to();
			auto fh = heh.face();
			double theta = mesh.calc_sector_angle(heh.next());
			E += pow(((paraCoordinate[vi] - paraCoordinate[vj]) - faceLt[fh]
				* (localCoordinate[fh][vi.idx()] - localCoordinate[fh][vj.idx()])).norm(), 2.0) / tan(theta);
		}
	}
	return E * 0.5;
}

bool MeshTools::ParameterizationARAP(Mesh& mesh, Mesh& paraMesh, int maxIter, double minEnergyVar)
{
	if (maxIter < 0 || minEnergyVar < 0.0)
	{
		return false;
	}
	// 初始参数化
	if (!Parameterization(mesh))
	{
		return false;
	}
	// 局部坐标
	CalcLocalCoordinate(mesh);
	// 稀疏系数矩阵
	SpMat coeMat;
	CalcCoeMat(mesh, coeMat);
	Eigen::SimplicialCholesky<SpMat> chol(coeMat);
	// 右侧向量
	Eigen::VectorXd bx, by;
	bx.resize(mesh.n_vertices());
	by.resize(mesh.n_vertices());
	// 解
	Eigen::VectorXd u, v;

	auto paraCoordinate = OpenMesh::getProperty<OpenMesh::VertexHandle, Eigen::Vector2d>(mesh, "paraCoordinate");
	
	double prevE, E;

	// 迭代一次计算初始能量

	// Local phase
	CalcSt(mesh);
	CalcLt(mesh);
	// Global phase
	CalcBVec(mesh, bx, by);
	u = chol.solve(bx);
	v = chol.solve(by);
	// 置内部顶点的参数坐标
	for (auto vh : mesh.vertices())
	{
		paraCoordinate[vh][0] = u[vh.idx()];
		paraCoordinate[vh][1] = v[vh.idx()];
	}
	prevE = CalcEnergy(mesh);

	// 开始迭代
	for (int i = 1; i < maxIter; i++)
	{
		// Local phase
		CalcSt(mesh);
		CalcLt(mesh);
		// Global phase
		CalcBVec(mesh, bx, by);
		u = chol.solve(bx);
		v = chol.solve(by);
		// 置内部顶点的参数坐标
		for (auto vh : mesh.vertices())
		{
			paraCoordinate[vh][0] = u[vh.idx()];
			paraCoordinate[vh][1] = v[vh.idx()];
		}
		// 判断能量变化
		E = CalcEnergy(mesh);
		if (prevE - E < minEnergyVar)
		{
			break;
		}
		prevE = E;
	}
	// 输出
	paraMesh.assign(mesh);
	for (auto vh : paraMesh.vertices())
	{
		auto originVh = mesh.vertex_handle(vh.idx());
		paraMesh.set_point(vh, Mesh::Point(paraCoordinate[originVh][0], paraCoordinate[originVh][1], 0));
	}
	return true;
}


