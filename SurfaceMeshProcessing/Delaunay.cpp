#include "Delaunay.h"

void Delaunay::CalcCircumcenter(Mesh& mesh, OpenMesh::SmartFaceHandle fh, Mesh::Point& cc)
{
	Eigen::Matrix2d A;
	Eigen::Vector2d b;
	// 构造第一条边的垂直平分线方程
	auto heh = fh.halfedge();
	auto p1 = mesh.point(heh.from());
	auto p2 = mesh.point(heh.to());
	A(0, 0) = 2.0 * (p2[0] - p1[0]);
	A(0, 1) = 2.0 * (p2[1] - p1[1]);
	b(0) = pow(p2[0], 2) + pow(p2[1], 2) - pow(p1[0], 2) - pow(p1[1], 2);
	// 构造第二条边的垂直平分线方程
	auto p3 = mesh.point(heh.next().to());
	A(1, 0) = 2.0 * (p3[0] - p2[0]);
	A(1, 1) = 2.0 * (p3[1] - p2[1]);
	b(1) = pow(p3[0], 2) + pow(p3[1], 2) - pow(p2[0], 2) - pow(p2[1], 2);
	// 解方程得外心
	auto solve = A.colPivHouseholderQr().solve(b);
	cc[0] = solve(0);
	cc[1] = solve(1);
}

void Delaunay::Lowson(Mesh& mesh)
{

}
