#include <QtCore>
#include <QMessageBox>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewerWidget.h"

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false),
	isLoadSource(false),
	isLoadTarget(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();
	bool read_OK = MeshTools::ReadMesh(mesh, filename);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

bool MeshViewerWidget::LoadSource(const std::string& filename)
{
	if (LoadMesh(filename))
	{
		sourceMesh.assign(mesh);
		isLoadSource = true;
		return true;
	}
	else
	{
		return false;
	}
}

void MeshViewerWidget::ShowSource()
{
	if (isLoadSource)
	{
		mesh.assign(sourceMesh);
		UpdateMesh();
		update();
	}
	else
	{
		QMessageBox::critical(NULL, "Warning", "No Source Mesh.");
	}
}

bool MeshViewerWidget::LoadTarget(const std::string& filename)
{
	if (LoadMesh(filename))
	{
		targetMesh.assign(mesh);
		isLoadTarget = true;
		return true;
	}
	else
	{
		return false;
	}
}

void MeshViewerWidget::ShowTarget()
{
	if (isLoadTarget)
	{
		mesh.assign(targetMesh);
		UpdateMesh();
		update();
	}
	else
	{
		QMessageBox::critical(NULL, "Warning", "No Target Mesh.");
	}
}

void MeshViewerWidget::DoMorphing()
{
	// 求source到target每个三角面的变形矩阵A
	// 对A分解，得到R·S，记录插值参数
	MeshTools::CalcFactorMatA(sourceMesh, targetMesh);

	// 选取预定义点，构造H
	int pdVidx;
	Eigen::MatrixXd H;
	MeshTools::CalcMatH(sourceMesh, pdVidx, H);
	//std::cout << "H:" << std::endl << H << std::endl;
	// 随t变化，构造G矩阵，解方程组，更新坐标
	auto theta = OpenMesh::getProperty<OpenMesh::FaceHandle, double>(sourceMesh, "rotationAngle");
	auto matS = OpenMesh::getProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(sourceMesh, "matS");
	double deltat = 0.01;
	for (double t = 0; t <= 1; t += deltat)
	{
		// 插值A(t)
		auto matAt = OpenMesh::makeTemporaryProperty<OpenMesh::FaceHandle, Eigen::Matrix2d>(sourceMesh);
		for (auto fh : sourceMesh.faces())
		{
			double itpTheta = theta[fh] * t;
			Eigen::Matrix2d R;
			R << cos(itpTheta), sin(itpTheta),
				-sin(itpTheta), cos(itpTheta);
			Eigen::Matrix2d S;
			S = (1 - t) * Eigen::Matrix2d::Identity() + t * matS[fh];
			matAt[fh] = R * S;
	//		std::cout << "interpolation A:" << std::endl << matAt[fh] << std::endl;
		}
		// 插值预定义点
		Mesh::Point pdP = sourceMesh.point(sourceMesh.vertex_handle(pdVidx));
		Mesh::Point targetP = targetMesh.point(targetMesh.vertex_handle(pdVidx));
		pdP += (targetP - pdP) * t;
	//	std::cout << "interpolation point:" << std::endl << pdP << std::endl;
		// 构造G
		Eigen::VectorXd Gx = Eigen::VectorXd::Zero(sourceMesh.n_vertices() - 1);
		Eigen::VectorXd Gy = Eigen::VectorXd::Zero(sourceMesh.n_vertices() - 1);
		for (auto Vi : sourceMesh.vertices())
		{
			if (Vi.idx() == pdVidx)
			{
				continue;
			}
			for (auto vvIter = sourceMesh.vv_begin(Vi); vvIter.is_valid(); vvIter++)
			{
				auto heH = sourceMesh.find_halfedge(Vi, *vvIter);
				if (sourceMesh.is_boundary(heH))
				{
					continue;
				}
				auto Vj = *vvIter;
				auto Vk = heH.next().to();
				auto pi = sourceMesh.point(Vi);
				auto pj = sourceMesh.point(Vj);
				auto pk = sourceMesh.point(Vk);

				double denominator = (pi[0] - pj[0]) * (pj[1] - pk[1]) - (pi[1] - pj[1]) * (pj[0] - pk[0]);

				double coeix = (pk[0] - pj[0]) / denominator;
				double coeiy = (pj[1] - pk[1]) / denominator;

				auto At = matAt[heH.face()];

				Gx[Vi.idx()] += At(0, 0) * coeiy + At(0, 1) * coeix;
				Gy[Vi.idx()] += At(1, 0) * coeiy + At(1, 1) * coeix;
				if (Vj.idx() == pdVidx)
				{
					double coejx = (pi[0] - pk[0]) / denominator;
					double coejy = (pk[1] - pi[1]) / denominator;
					Gx[Vi.idx()] -= pdP[0] * (coejx * coeix + coejy * coeiy);
					Gy[Vi.idx()] -= pdP[1] * (coejx * coeix + coejy * coeiy);
				}
				else if (Vk.idx() == pdVidx)
				{
					double coekx = (pj[0] - pi[0]) / denominator;
					double coeky = (pi[1] - pj[1]) / denominator;
					Gx[Vi.idx()] -= pdP[0] * (coekx * coeix + coeky * coeiy);
					Gy[Vi.idx()] -= pdP[1] * (coekx * coeix + coeky * coeiy);
				}
			}
		}
	//	std::cout << "Gx:" << std::endl << Gx << std::endl << "Gy:" << std::endl << Gy << std::endl;
		// 解方程组
		Eigen::VectorXd solx = H.colPivHouseholderQr().solve(Gx);
		Eigen::VectorXd soly = H.colPivHouseholderQr().solve(Gy);
	//	std::cout << "solution x:" << std::endl << solx << std::endl << "solution y:" << std::endl << soly << std::endl;
		{
			/*
			auto Vi = sourceMesh.vertex_handle(0);
			auto Vj = sourceMesh.vertex_handle(1);
			auto Vk = sourceMesh.vertex_handle(2);
			auto pi = sourceMesh.point(Vi);
			auto pj = sourceMesh.point(Vj);
			auto pk = sourceMesh.point(Vk);

			double denominator = (pi[0] - pj[0]) * (pj[1] - pk[1]) - (pi[1] - pj[1]) * (pj[0] - pk[0]);

			double coeix = (pk[0] - pj[0]) / denominator;
			double coejx = (pi[0] - pk[0]) / denominator;
			double coekx = (pj[0] - pi[0]) / denominator;

			double coeiy = (pj[1] - pk[1]) / denominator;
			double coejy = (pk[1] - pi[1]) / denominator;
			double coeky = (pi[1] - pj[1]) / denominator;

			double b1 = solx[0] * coeiy + solx[1] * coejy + pdP[0] * coeky;
			double b2 = solx[0] * coeix + solx[1] * coejx + pdP[0] * coekx;
			double b3 = soly[0] * coeiy + soly[1] * coejy + pdP[1] * coeky;
			double b4 = soly[0] * coeix + soly[1] * coejx + pdP[1] * coekx;
			*/
		//	std::cout << "B:" << std::endl << b1 << "," << b2 << std::endl << b3 << "," << b4 << std::endl;
		}
		// 更新顶点位置
		mesh.set_point(mesh.vertex_handle(pdVidx), pdP);
		for (int idx = 0; idx < pdVidx; idx++)
		{
			mesh.set_point(mesh.vertex_handle(idx), Mesh::Point(solx[idx], soly[idx], 0));
		}
		UpdateMesh();
		update();
		Sleep(100);
	}
}

void MeshViewerWidget::Clear(void)
{
	mesh.clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	mesh.update_normals();
	if (mesh.vertices_empty())
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : mesh.edges())
	{
		double len = mesh.calc_edge_length(eh);
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return MeshTools::WriteMesh(mesh, filename, DBL_DECIMAL_DIG);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (!mesh.vertices_empty())
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (mesh.n_vertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::DrawPoints(void) const
{
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		glNormal3dv(mesh.normal(vh).data());
		glVertex3dv(mesh.point(vh).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : mesh.faces())
	{
		glNormal3dv(mesh.normal(fh).data());
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glVertex3dv(mesh.point(fvh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawSmooth(void) const
{
	glColor3d(0.8, 0.8, 0.8);
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	for (const auto& fh : mesh.faces())
	{
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			auto vh0 = mesh.from_vertex_handle(heh);
			auto vh1 = mesh.to_vertex_handle(heh);
			glNormal3dv(mesh.normal(vh0).data());
			glVertex3dv(mesh.point(vh0).data());
			glNormal3dv(mesh.normal(vh1).data());
			glVertex3dv(mesh.point(vh1).data());
		}
	}
	glEnd();
	glLineWidth(linewidth);
}
