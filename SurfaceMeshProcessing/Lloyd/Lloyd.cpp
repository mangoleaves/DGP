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
	// Ԥ�ȼ�����ķ����������ģ����
	mesh.request_face_normals();
	centroid = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh, "centroid");
	area = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	for (auto fh : mesh.faces())
	{
		(*centroid)[fh] = mesh.calc_face_centroid(fh);
		(*area)[fh] = mesh.calc_face_area(fh);
	}

	// ÿ����ķ��ֵ࣬ΪProxy��idx�����Ǵ��ڵ���0��Ϊ-1��ʾδ���ࡣ
	partition = &OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	for (auto fh : mesh.faces())
	{
		(*partition)[fh] = -1;
	}

	// �������ȼ��ѡȡK������Ϊseed triangles
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
	// ��ÿ��seed triangle����������������ȶ��У�labelΪseed triangle��partition idx��priority���������proxy֮���error
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

	// ����ÿ��region�����ġ��������������
	for (auto fh : mesh.faces())
	{
		int proxyIdx = (*partition)[fh];
		regionArea[proxyIdx] += (*area)[fh];
		regionX[proxyIdx] += (*centroid)[fh];
		regionN[proxyIdx] += mesh.normal(fh);
	}
	// ����proxy
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
	// Ϊÿ��proxyѡȡ��ӽ�����triangleΪseed
	seedtris.clear();
	std::vector<double> minError(K, INFINITY);

	// ����ÿ���棬�ҵ���ÿ��proxy֮��error��С���棬��Ϊ��proxy��seed triangle
	// error��С����һ������Geometry Partition������proxy����
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
	// ����partition
	for (auto fh : mesh.faces())
	{
		if (fh.idx() != seedtris[(*partition)[fh]].idx())
		{
			(*partition)[fh] = -1;
		}
	}
}

