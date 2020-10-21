#include "Lloyd.h"

void Lloyd::DoLloyd()
{
	InitialSeeding();
	bool isUpdated;
	do
	{
		isUpdated = false;
		GeometryPartition();
		isUpdated = ProxyFitting();
		GenerateSeeds();
	} while (isUpdated);
}

void Lloyd::InitialSeeding()
{
	// 预先计算面的法向量，重心，面积
	mesh.request_face_normals();
	centroid = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh, "centroid");
	area = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	for (auto fh : mesh.faces())
	{
		(*centroid)[fh] = mesh.calc_face_centroid(fh);
		(*area)[fh] = mesh.calc_face_area(fh);
	}

	// 每个面的分类，值为Proxy的idx，总是大于等于0，为-1表示未分类。
	partition = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	for (auto fh : mesh.faces())
	{
		(*partition)[fh] = -1;
	}

	// 按索引等间距选取K个面作为seed triangles
	int nface = mesh.n_faces();
	int proxyIdx = 0;
	for (int i = 0; i < nface; i += nface / K)
	{
		auto fh = mesh.face_handle(i);
		seedtris.push_back(fh);
		proxys.push_back(Proxy(proxyIdx, (*centroid)[fh], mesh.normal(fh)));
		(*partition)[fh] = proxyIdx;
		proxyIdx++;
	}
}

double Lloyd::CalcError(OpenMesh::FaceHandle fh, Proxy proxy)
{
	return (*area)[fh] * (mesh.normal(fh) - proxy.N).norm();
}

void Lloyd::GeometryPartition()
{
	TriangleQueue triq;
	// 将每个seed triangle的三个邻面加入优先队列，label为seed triangle的partition idx，priority是三角面和proxy之间的error
	for (auto fh : seedtris)
	{
		for (auto ffIter = mesh.ff_begin(fh); ffIter.is_valid(); ffIter++)
		{
			triq.push(QueueElem(*ffIter, (*partition)[fh], CalcError(*ffIter, proxys[(*partition)[fh]])));
		}
	}
	
	while (!triq.empty())
	{
		QueueElem elem = triq.top();
		triq.pop();

		if ((*partition)[elem.fh] >= 0)
		{
			continue;
		}
		else
		{
			(*partition)[elem.fh] = elem.label;
			for (auto ffIter = mesh.ff_begin(elem.fh); ffIter.is_valid(); ffIter++)
			{
				if ((*partition)[*ffIter] == -1)
				{
					triq.push(QueueElem(*ffIter, (*partition)[elem.fh], CalcError(*ffIter, proxys[(*partition)[elem.fh]])));
				}
			}
		}
	}
}

bool Lloyd::ProxyFitting()
{
	std::vector<double> regionArea(K, 0.0);
	std::vector<Mesh::Point> regionX(K, Mesh::Point(0, 0, 0));
	std::vector<Mesh::Normal> regionN(K, Mesh::Normal(0, 0, 0));

	// 计算每个region的重心、法向量、面积和
	for (auto fh : mesh.faces())
	{
		int proxyIdx = (*partition)[fh];
		regionArea[proxyIdx] += (*area)[fh];
		regionX[proxyIdx] += (*centroid)[fh];
		regionN[proxyIdx] += mesh.normal(fh);
	}
	// 更新proxy
	bool isUpdated = false;
	for (int i = 0; i < proxys.size(); i++)
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
	// 为每个proxy选取最接近它的triangle为seed
	seedtris.clear();
	std::vector<double> minError(K, INFINITY);

	// 遍历每个面，找到与每个proxy之间error最小的面，作为该proxy的seed triangle
	// error最小的面一定是在Geometry Partition后属于proxy的面
	for (auto fh : mesh.faces())
	{
		int proxyIdx = (*partition)[fh];
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
		if (fh.idx() != seedtris[(*partition)[fh]].idx())
		{
			(*partition)[fh] = -1;
		}
	}
}

