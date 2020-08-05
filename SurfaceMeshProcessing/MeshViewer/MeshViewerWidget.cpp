#include <QtCore>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewerWidget.h"

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
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
		auto vertexState = OpenMesh::getOrMakeProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
		for (auto vh : mesh.vertices())
		{
			vertexState[vh] = NotSelected;
		}
		selectMode = NoSelect;
		UpdateMesh();
		update();
		return true;
	}
	return false;
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
	avgEdgeLength = avelen / mesh.n_vertices();

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

void MeshViewerWidget::SetSMFixed(void)
{
	selectMode = SelectFixed;
}

void MeshViewerWidget::SetSMCustom(void)
{
	selectMode = SelectCustom;
}

void MeshViewerWidget::SetSMMove(void)
{
	selectMode = Move;
}

void MeshViewerWidget::SetSMNoSelect(void)
{
	selectMode = NoSelect;
}

void MeshViewerWidget::ClearSelected(void)
{
	auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
	for (auto vh : mesh.vertices())
	{
		vertexState[vh] = NotSelected;
	}
	update();
}

bool MeshViewerWidget::event(QEvent* _event)
{
	if (_event->type() == QEvent::MouseButtonPress)
	{
		if (selectMode == Move)
		{
			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			double depth;
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), objCor, depth);

			OpenMesh::VertexHandle minVh;
			if (NearestVertex(objCor, minVh))
			{
				auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
				if (vertexState[minVh] == Custom)
				{
					isMovable = true;
					moveDepth = depth;
					lastObjCor = objCor;
					return true;
				}
			}
		}
	}
	else if (_event->type() == QEvent::MouseMove)
	{
		if (selectMode == Move && isMovable)
		{
			auto e = static_cast<QMouseEvent*>(_event);
			QPoint winCor = e->pos();
			OpenMesh::Vec3d objCor;
			WinCor2ObjCor(winCor.x(), winCor.y(), moveDepth, objCor);

			auto moveVec = objCor - lastObjCor;
			lastObjCor = objCor;
			Mesh deformedMesh;
			deformedMesh.assign(mesh);
			auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
			for (auto vh_ : deformedMesh.vertices())
			{
				if (vertexState[mesh.vertex_handle(vh_.idx())] == Custom)
				{
					deformedMesh.set_point(vh_, deformedMesh.point(vh_) + moveVec);
				}
			}
			MeshTools::Deform(mesh, deformedMesh);
			MeshTools::AssignPoints(mesh, deformedMesh);
			update();
			return true;
		}
	}
	else if (_event->type() == QEvent::MouseButtonRelease)
	{
		if (selectMode == Move && isMovable)
		{
			isMovable = false;
			return true;
		}
	}
	return QGLViewerWidget::event(_event);
}

void MeshViewerWidget::mouseDoubleClickEvent(QMouseEvent* _event)
{
	switch (selectMode)
	{
	case NoSelect:
		break;
	case SelectFixed:
	case SelectCustom: 
	{
		QPoint winCor = _event->pos();
		double depth;
		OpenMesh::Vec3d objCor;
		WinCor2ObjCor(winCor.x(), winCor.y(), objCor, depth);
		OpenMesh::VertexHandle minVh;

		if (NearestVertex(objCor,minVh))
		{
			auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
			if (selectMode == SelectFixed)
			{
				if (vertexState[minVh] == NotSelected)
				{
					vertexState[minVh] = Fixed;
				}
				else if (vertexState[minVh] == Fixed)
				{
					vertexState[minVh] = NotSelected;
				}
			}
			else if (selectMode == SelectCustom)
			{
				if (vertexState[minVh] == NotSelected)
				{
					vertexState[minVh] = Custom;
				}
				else if (vertexState[minVh] == Custom)
				{
					vertexState[minVh] = NotSelected;
				}
			}
			update();
		}
		break;
	}
	case Move:
		break;
	default:
		break;
	}
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

bool MeshViewerWidget::NearestVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh)
{
	double maxAllowedDis = avgEdgeLength * 0.5;
	double minDis = INFINITY;
	OpenMesh::VertexHandle mv;

	for (auto vh : mesh.vertices())
	{
		auto vp = mesh.point(vh);
		double dis = (objCor - vp).norm();
		if (dis < minDis)
		{
			minDis = dis;
			mv = vh;
		}
	}

	if (minDis <= maxAllowedDis)
	{
		minVh = mv;
		return true;
	}
	else
	{
		return false;
	}
}

void MeshViewerWidget::DrawPoints(void)
{
	auto vertexState = OpenMesh::getProperty<OpenMesh::VertexHandle, VertexState>(mesh, "vertexState");
	glColor3d(0.2, 0.2, 0.2);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		if (vertexState[vh] == NotSelected)
		{
			glNormal3dv(mesh.normal(vh).data());
			glVertex3dv(mesh.point(vh).data());
		}
	}
	glEnd();
	glColor3d(0.0, 0.0, 0.7);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		if (vertexState[vh] == Fixed)
		{
			glNormal3dv(mesh.normal(vh).data());
			glVertex3dv(mesh.point(vh).data());
		}
	}
	glColor3d(0.7, 0.7, 0.0);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		if (vertexState[vh] == Custom)
		{
			glNormal3dv(mesh.normal(vh).data());
			glVertex3dv(mesh.point(vh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void)
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
	DrawPoints();
}

void MeshViewerWidget::DrawHiddenLines()
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

void MeshViewerWidget::DrawFlatLines(void)
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

void MeshViewerWidget::DrawFlat(void)
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
	DrawPoints();
}

void MeshViewerWidget::DrawSmooth(void)
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
	DrawPoints();
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
