#include "Remeshing.h"


void Remeshing::DoRemeshing(Mesh& mesh, double targetLength)
{
	// 长度上下界
	double low = 4.0 / 5.0 * targetLength;
	double high = 4.0 / 3.0 * targetLength;

	// AABB tree
	Tree tree;
	std::list<Triangle> triangles;
	ConstructAABBTree(mesh, tree, triangles);

	// 主循环
	for (int i = 0; i < 10; i++)
	{
		SplitLongEdges(mesh, high);
		CollapseShortEdges(mesh, low, high);
		EqualizeValences(mesh);
		TangentialRelaxation(mesh);
		ProjectToSurface(mesh, tree);
	}
}

void Remeshing::SplitLongEdges(Mesh& mesh, double high)
{
	bool hasNewEdge = false;	// 是否有新边生成
	do {
		hasNewEdge = false;
		for (auto eh : mesh.edges())
		{
			if (mesh.calc_edge_length(eh) > high)
			{
				hasNewEdge = true;
				Mesh::Point midPoint = (mesh.point(eh.v0()) + mesh.point(eh.v1())) / 2.0;
				auto vh = mesh.add_vertex(midPoint);
				mesh.split_edge(eh, vh);
			}
		}
	} while (hasNewEdge);
}

void Remeshing::CollapseShortEdges(Mesh& mesh, double low, double high)
{
	bool hasNewEdge = false;
	do {
		hasNewEdge = false;
		for (auto eh : mesh.edges())
		{
			if (mesh.calc_edge_length(eh) < low)
			{
				bool collapseOK = true;
				auto va = eh.v0();
				auto vb = eh.v1();
				auto heh = mesh.find_halfedge(va, vb);
				for (auto vvIter = mesh.vv_begin(va); vvIter.is_valid(); vvIter++)
				{
					if ((mesh.point(*vvIter) - mesh.point(vb)).length() > high)
					{
						collapseOK = false;
					}
				}
				if (collapseOK && mesh.is_collapse_ok(heh))
				{
					hasNewEdge = true;
					mesh.collapse(heh);
				}
			}
		}
	} while (hasNewEdge);
}

int Remeshing::TargetValence(OpenMesh::SmartVertexHandle vh)
{
	return vh.is_boundary() ? 4 : 6;
}

void Remeshing::EqualizeValences(Mesh& mesh)
{
	for (auto eh : mesh.edges())
	{
		auto va = eh.v0();
		auto vb = eh.v1();
		auto vc = eh.h0().next().to();
		auto vd = eh.h1().next().to();

		int vala = mesh.valence(va);
		int valb = mesh.valence(vb);
		int valc = mesh.valence(vc);
		int vald = mesh.valence(vd);

		int tara = TargetValence(va);
		int tarb = TargetValence(vb);
		int tarc = TargetValence(vc);
		int tard = TargetValence(vd);

		int deviationPre = abs(vala - tara) + abs(valb - tarb) + abs(valc - tarc) + abs(vald - tard);
		int deviationPost = abs(vala - 1 - tara) + abs(valb - 1 - tarb) + abs(valc + 1 - tarc) + abs(vald + 1 - tard);

		if (deviationPost < deviationPre && mesh.is_flip_ok(eh))
		{
			mesh.flip(eh);
		}
	}
}

void Remeshing::TangentialRelaxation(Mesh& mesh)
{
	auto q = OpenMesh::makeTemporaryProperty<OpenMesh::VertexHandle, Mesh::Point>(mesh);
	assert(mesh.has_vertex_normals());

	for (auto vh : mesh.vertices())
	{
		q[vh] = Mesh::Point(0, 0, 0);
		double N = 0;
		for (auto vvIter = mesh.vv_begin(vh); vvIter.is_valid(); vvIter++)
		{
			q[vh] += mesh.point(*vvIter);
			N += 1;
		}
		q[vh] /= N;
	}

	for (auto vh : mesh.vertices())
	{
		Mesh::Point pv = mesh.point(vh);
		Mesh::Normal nv = mesh.calc_vertex_normal(vh);
		mesh.set_point(vh, q[vh] + nv.dot(pv - q[vh]) * nv);
	}
}

void Remeshing::ConstructAABBTree(Mesh& mesh, Tree& tree, std::list<Triangle>& triangles)
{
	// 遍历原网格，将所有三角面存放在triangles中，用于初始化AABB tree
	for (auto fh : mesh.faces())
	{
		std::vector<Point> ps;
		for (auto vh : fh.vertices())
		{
			Mesh::Point p = mesh.point(vh);
			ps.push_back(Point(p[0], p[1], p[2]));
		}
		triangles.push_back(Triangle(ps[0], ps[1], ps[2]));
	}
	
	tree.insert(triangles.begin(), triangles.end());
	tree.build();
	tree.accelerate_distance_queries();
}

void Remeshing::ProjectToSurface(Mesh& mesh, Tree& tree)
{
	for (auto vh : mesh.vertices())
	{
		Mesh::Point vp = mesh.point(vh);
		Point pointQuery(vp[0], vp[1], vp[2]);
		Point closetPoint = tree.closest_point(pointQuery);
		mesh.set_point(vh, Mesh::Point(closetPoint[0], closetPoint[1], closetPoint[2]));
	}
}


