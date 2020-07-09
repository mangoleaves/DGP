#pragma once

#include <vector>
#include <algorithm>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshDefinition.h>



/* 使用Dijkstra算法寻找两点间的最短路径 */
std::vector<OpenMesh::VertexHandle>* FindShortestPath(Mesh& mesh, int beginIndex, int endIndex)
{

}
