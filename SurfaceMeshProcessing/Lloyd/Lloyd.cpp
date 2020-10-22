#include "Lloyd.h"

void Lloyd::DoLloyd(Mesh& outputMesh, int nIter)
{
	InitialSeeding();
	GeometryPartitioning();
	ProxyFitting();
	for (int i = 0; i < nIter; i++)
	{
		GenerateSeeds();
		GeometryPartitioning();
		ProxyFitting();
	}
	/*
	do
	{
		GeometryPartition();
		if (ProxyFitting())
			GenerateSeeds();
		else
			break;
	} while (true);
	*/
	auto par = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	auto outpar = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(outputMesh, "partition");
	for (auto fh : mesh.faces())
	{
		outpar[outputMesh.face_handle(fh.idx())] = par[fh];
	}
}

void Lloyd::InitialSeeding()
{
	// 预先计算面的法向量，重心，面积
	mesh.request_face_normals();
	auto centroid = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh, "centroid");
	auto area = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	for (auto fh : mesh.faces())
	{
		centroid[fh] = mesh.calc_face_centroid(fh);
		area[fh] = mesh.calc_face_area(fh);
	}

	// 每个面的分类，值为Proxy的idx，总是大于等于0，为-1表示未分类。
	auto partition = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	for (auto fh : mesh.faces())
	{
		partition[fh] = -1;
	}

	// 按索引等间距选取K个面作为seed triangle，生成proxy并分配
	int nface = mesh.n_faces();
	int proxyIdx, i;
	for (i = 0, proxyIdx = 0; i < nface; i += (nface - 1) / (K - 1), proxyIdx++)
	{
		auto fh = mesh.face_handle(i);
		seedtris.push_back(fh);
		proxys.push_back(Proxy(proxyIdx, centroid[fh], mesh.normal(fh)));
		partition[fh] = proxyIdx;
	}
}

double Lloyd::CalcError(OpenMesh::SmartFaceHandle& fh, Proxy& proxy)
{
	auto area = OpenMesh::getProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	return area[fh] * pow((mesh.normal(fh) - proxy.N).norm(), 2);
}

void Lloyd::GeometryPartitioning()
{
	auto area = OpenMesh::getProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	auto partition = OpenMesh::getProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	TriangleQueue triq;
	// 将每个seed triangle的三个邻面加入优先队列，label为seed triangle的partition idx，priority是三角面和proxy之间的error
	for (auto fh : seedtris)
	{
		for (auto ffIter = mesh.ff_begin(fh); ffIter.is_valid(); ffIter++)
		{
			triq.push(QueueElem(*ffIter, partition[fh], CalcError(*ffIter, proxys[partition[fh]])));
		}
	}
	// 优先队列非空时，取最小权值的面。
	// 若其已被分配，则跳过；
	// 若未被分配，则分配到label对应的proxy中，并将相邻的三角面按相同的方法加入优先队列。
	while (!triq.empty())
	{
		QueueElem minElem = triq.top();
		triq.pop();

		if (partition[minElem.fh] >= 0)
		{
			continue;
		}
		else
		{
			partition[minElem.fh] = minElem.label;
			for (auto ffIter = mesh.ff_begin(minElem.fh); ffIter.is_valid(); ffIter++)
			{
				if (partition[*ffIter] == -1)
				{
					triq.push(QueueElem(*ffIter, partition[minElem.fh], CalcError(*ffIter, proxys[partition[minElem.fh]])));
				}
			}
		}
	}
}

bool Lloyd::ProxyFitting()
{
	auto centroid = OpenMesh::getProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh, "centroid");
	auto area = OpenMesh::getProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	auto partition = OpenMesh::getProperty<OpenMesh::FaceHandle, int>(mesh, "partition");

	std::vector<double> regionArea(K, 0.0);
	std::vector<Mesh::Point> regionX(K, Mesh::Point(0, 0, 0));
	std::vector<Mesh::Normal> regionN(K, Mesh::Normal(0, 0, 0));

	// 计算每个region的重心、法向量、面积和
	for (auto fh : mesh.faces())
	{
		int proxyIdx = partition[fh];
		regionArea[proxyIdx] += area[fh];
		regionX[proxyIdx] += area[fh] * centroid[fh];
		regionN[proxyIdx] += area[fh] * mesh.normal(fh);
	}
	// 更新proxy
	bool isUpdated = false;
	for (int i = 0; i < K; i++)
	{
		if (proxys[i].X != regionX[i] / regionArea[i])
		{
			proxys[i].X = regionX[i] / regionArea[i];
			isUpdated = true;
		}
		if (proxys[i].N != regionN[i] / regionArea[i])
		{
			proxys[i].N = regionN[i] / regionArea[i];
			isUpdated = true;
		}
	}
	return isUpdated;
}

void Lloyd::GenerateSeeds()
{
	auto partition = OpenMesh::getProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	// 为每个proxy选取最接近它的triangle为seed
	std::vector<double> minError(K, INFINITY);

	// 遍历每个面，找到与每个proxy之间error最小的面，作为该proxy的seed triangle
	// error最小的面一定是在Geometry Partition后属于proxy的面
	for (auto fh : mesh.faces())
	{
		int proxyIdx = partition[fh];
		double error = CalcError(fh, proxys[proxyIdx]);
		if (error < minError[proxyIdx])
		{
			minError[proxyIdx] = error;
			seedtris[proxyIdx] = fh;
		}
	}
	// 重置partition
	for (auto fh : mesh.faces())
	{
		if (fh.idx() != seedtris[partition[fh]].idx())
		{
			partition[fh] = -1;
		}
	}
}

