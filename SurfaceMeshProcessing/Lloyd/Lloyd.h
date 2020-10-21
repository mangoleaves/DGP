#pragma once
#include <iostream>
#include <cmath>
#include <queue>
#include "MeshDefinition.h"

class Lloyd
{
public:
	Lloyd() = default;
	Lloyd(Mesh& m, int k) :mesh(m), K(k) {}

	void DoLloyd();
private:
	// types
	typedef struct pxy {
		int idx;
		Mesh::Point X;
		Mesh::Normal N;

		pxy() = default;
		pxy(int i, Mesh::Point x, Mesh::Normal n) :idx(i), X(x), N(n) {}
	}Proxy;

	typedef std::vector<Proxy> Proxys;
	typedef std::vector<OpenMesh::FaceHandle> SeedTriangles;

	typedef struct qe {
		OpenMesh::FaceHandle fh;
		int label;
		double priority;

		qe() = default;
		qe(OpenMesh::FaceHandle f, int l, double p) :fh(f), label(l), priority(p) {}
	}QueueElem;

	class compare {
	public:
		bool operator()(QueueElem e1, QueueElem e2)
		{
			return e1.priority > e2.priority;
		}
	};

	typedef std::priority_queue<QueueElem, compare> TriangleQueue;

	// private data
	Mesh mesh;
	OpenMesh::PropertyManager<OpenMesh::FPropHandleT<MeshTraits::Normal>, int>* centroid;
	OpenMesh::PropertyManager<OpenMesh::FPropHandleT<double>, int>* area;
	int K;
	Proxys proxys;
	SeedTriangles seedtris;
	OpenMesh::PropertyManager<OpenMesh::FPropHandleT<int>, int>* partition;

	// private functions
	double CalcError(OpenMesh::FaceHandle fh, Proxy proxy);
	void InitialSeeding();
	void GeometryPartition();
	bool ProxyFitting();
	void GenerateSeeds();
};



