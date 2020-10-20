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

void Delaunay::CalcCircumcenter(Mesh& mesh, OpenMesh::SmartFaceHandle fh, Mesh::Point& cc, double& radius)
{
	CalcCircumcenter(mesh, fh, cc);
	radius = (cc - mesh.point(fh.halfedge().from())).norm();
}

void Delaunay::Lowson(Mesh& mesh)
{
	// 标记一条边是否需要检查其是否满足空圆性质
	auto needCheck = OpenMesh::makeTemporaryProperty<OpenMesh::EdgeHandle, bool>(mesh);
	for (auto eh : mesh.edges())
	{
		needCheck[eh] = true;
	}
	// 检查每条边相邻的两个圆是否满足空圆性质
	bool hasCheck;
	do
	{
		hasCheck = false;
		for (auto eh : mesh.edges())
		{
			if (needCheck[eh])
			{
				hasCheck = true;
				// 计算相邻的一面的外心cc和半径radius
				Mesh::Point cc;
				double radius;
				CalcCircumcenter(mesh, eh.h0().face(), cc, radius);
				// 若第4个点和cc的距离小于radius，则需要flip
				if ((mesh.point(eh.h1().next().to()) - cc).norm() < radius)
				{
					mesh.flip(eh);
					// 受影响的边为外围的4条，标记为needCheck
					needCheck[eh.h0().next().edge()] = true;
					needCheck[eh.h0().next().next().edge()] = true;
					needCheck[eh.h1().next().edge()] = true;
					needCheck[eh.h1().next().next().edge()] = true;
				}
				needCheck[eh] = false;
			}
		}
	} while (hasCheck);
}

void Delaunay::DelaunayTriangulation(Mesh& mesh)
{
	Lowson(mesh);
}
