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
	void CreateMorphWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void LoadSourceSignal();
	void ShowSourceSignal();
	void LoadTargetSignal();
	void ShowTargetSignal();
	void InitMorphingSignal(double interval);
	void DoMorphingSignal();
private slots:
	void EmitInitSignal(void);
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QWidget* morphWidget;
	QPushButton* loadSource, * showSource, * loadTarget, * showTarget, * initMorphing, * doMorphing;
	QLineEdit* intervalLE;
};
