#pragma once
#include <QString>
#include <QEvent>
#include <QMouseEvent>
#include "QGLViewerWidget.h"
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include "MeshDefinition.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	enum SelectMode {
		NoSelect,
		SelectFixed,
		SelectCustom,
		Move
	};
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
	void SetSMFixed(void);
	void SetSMCustom(void);
	void SetSMMove(void);
	void SetSMNoSelect(void);
	void ClearSelected(void);
protected:
	virtual bool event(QEvent* _event) override;
	virtual void mouseDoubleClickEvent(QMouseEvent* _event) override;
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);
	bool NearestVertex(OpenMesh::Vec3d objCor, OpenMesh::VertexHandle& minVh);

private:
	void DrawPoints(void);
	void DrawWireframe(void);
	void DrawHiddenLines(void);
	void DrawFlatLines(void);
	void DrawFlat(void);
	void DrawSmooth(void);
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:
	Mesh mesh;
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	double avgEdgeLength;
	SelectMode selectMode;
	bool isMovable;
	double moveDepth;
	OpenMesh::Vec3d lastObjCor;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
};
