#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include <Eigen/Dense>
#include "MeshDefinition.h"

class Delaunay
{
public:
	static void CalcCircumcenter(Mesh& mesh, OpenMesh::SmartFaceHandle fh, Mesh::Point& cc);
	static void Lowson(Mesh& mesh);
};

