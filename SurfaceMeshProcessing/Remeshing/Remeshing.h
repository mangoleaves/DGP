#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include "MeshDefinition.h"

class Remeshing
{
public:
	// AABB tree types
	typedef CGAL::Simple_cartesian<double> K;
	typedef K::FT FT;
	typedef K::Ray_3 Ray;
	typedef K::Line_3 Line;
	typedef K::Point_3 Point;
	typedef K::Triangle_3 Triangle;
	typedef std::list<Triangle>::iterator Iterator;
	typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
	typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
	typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

	static void DoRemeshing(Mesh& mesh, double targetLength);
	static void SplitLongEdges(Mesh& mesh, double high);
	static void CollapseShortEdges(Mesh& mesh, double low, double high);
	static inline int TargetValence(OpenMesh::SmartVertexHandle vh);
	static void EqualizeValences(Mesh& mesh);
	static void TangentialRelaxation(Mesh& mesh);
	static void ConstructAABBTree(Mesh& mesh, Tree& tree, std::list<Triangle>& triangles);
	static void ProjectToSurface(Mesh& mesh, Tree& tree);
};

