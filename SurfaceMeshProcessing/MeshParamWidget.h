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
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void DoLloydSignal(int K, int nIter);
private slots:
	void EmitDoLloydSignal();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QLineEdit* kLE;
	QLineEdit* nIterLE;
	QPushButton* DoLloydBtn;
};
