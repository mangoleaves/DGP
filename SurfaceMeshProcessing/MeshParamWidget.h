#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>

class MeshParamWidget : public QWidget
{
	Q_OBJECT

public:
	MeshParamWidget(QWidget *parent = 0);
	~MeshParamWidget(void);
private:
	void CreateTabWidget(void);
	void CreateCurvatureWidget(void);
	void CreateDenoiseWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void ShowMC();
	void ShowAMC();
	void ShowGC();
	void ShowOrigin();
	void AddNoise();
	void Denoise();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	// Curvature
	QPushButton* MCButton;
	QPushButton* AMCButton;
	QPushButton* GCButton;
	QWidget* curvatureWidget;
	// Denoise
	QPushButton* showOriginButton;
	QPushButton* addNoiseButton;
	QPushButton* denoiseButton;
	QWidget* denoiseWidget;
public:		// 偷一下懒，就直接把LineEdit暴露了
	QLineEdit* noiseProportionLE;
	QLineEdit* sigmaSLE;
	QLineEdit* nMILE;
	QLineEdit* vMILE;
};
