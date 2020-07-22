#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "MeshDefinition.h"
#include "ColorMap.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
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
	// Denoise
	void ShowOrigin(void);
	void AddNoise(double proportion);
	void Denoise(double sS, int nMI, int vMI);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawSmooth(void) const;
	void DrawMeanCurvature(void);
	void DrawAbsoluteMeanCurvature(void);
	void DrawGaussianCurvature(void);
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:
	Mesh mesh;
	// Denoise
	Mesh originalMesh;
	Mesh noisyMesh;
	Mesh denoisedMesh;
	double noiseProportion;
	double sigmaS;
	int normalMaxIter, vertexMaxIter;
	// Curvature
	bool hasCalcLAR;
	bool hasCalcMeanCurvature;
	bool hasCalcAbsoluteMeanCurvature;
	bool hasCalcGaussianCurvature;
	bool hasCalcMCColor;
	bool hasCalcAMCColor;
	bool hasCalcGCColor;

	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
};
