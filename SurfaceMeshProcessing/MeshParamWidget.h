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
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void ShowMC();
	void ShowAMC();
	void ShowGC();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QPushButton* MCButton;
	QPushButton* AMCButton;
	QPushButton* GCButton;
	QWidget* curvatureWidget;
};
