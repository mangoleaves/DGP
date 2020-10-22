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
	// Ԥ�ȼ�����ķ����������ģ����
	mesh.request_face_normals();
	auto centroid = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, Mesh::Point>(mesh, "centroid");
	auto area = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, double>(mesh, "area");
	for (auto fh : mesh.faces())
	{
		centroid[fh] = mesh.calc_face_centroid(fh);
		area[fh] = mesh.calc_face_area(fh);
	}

	// ÿ����ķ��ֵ࣬ΪProxy��idx�����Ǵ��ڵ���0��Ϊ-1��ʾδ���ࡣ
	auto partition = OpenMesh::getOrMakeProperty<OpenMesh::FaceHandle, int>(mesh, "partition");
	for (auto fh : mesh.faces())
	{
		partition[fh] = -1;
	}

	// �������ȼ��ѡȡK������Ϊseed triangle������proxy������
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
	// ��ÿ��seed triangle����������������ȶ��У�labelΪseed triangle��partition idx��priority���������proxy֮���error
	for (auto fh : seedtris)
	{
		for (auto ffIter = mesh.ff_begin(fh); ffIter.is_valid(); ffIter++)
		{
			triq.push(QueueElem(*ffIter, partition[fh], CalcError(*ffIter, proxys[partition[fh]])));
		}
	}
	// ���ȶ��зǿ�ʱ��ȡ��СȨֵ���档
	// �����ѱ����䣬��������
	// ��δ�����䣬����䵽label��Ӧ��proxy�У��������ڵ������水��ͬ�ķ����������ȶ��С�
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

	// ����ÿ��region�����ġ��������������
	for (auto fh : mesh.faces())
	{
		int proxyIdx = partition[fh];
		regionArea[proxyIdx] += area[fh];
		regionX[proxyIdx] += area[fh] * centroid[fh];
		regionN[proxyIdx] += area[fh] * mesh.normal(fh);
	}
	// ����proxy
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
	// Ϊÿ��proxyѡȡ��ӽ�����triangleΪseed
	std::vector<double> minError(K, INFINITY);

	// ����ÿ���棬�ҵ���ÿ��proxy֮��error��С���棬��Ϊ��proxy��seed triangle
	// error��С����һ������Geometry Partition������proxy����
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
	// ����partition
	for (auto fh : mesh.faces())
	{
		if (fh.idx() != seedtris[partition[fh]].idx())
		{
			partition[fh] = -1;
		}
	}
}

