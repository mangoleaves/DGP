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
	void CreateSimpWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void SimpSignal(int Nv, double minCost);
private slots:
	void EmitSimpSignal(void);
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QWidget* simpWidget;
	QLineEdit* NvLE,* minCostLE;
	QPushButton* simpBtn;
};
