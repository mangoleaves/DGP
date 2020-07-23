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
	void CreateParaWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void ShowOriginSignal();
	void ParaSignal();
private:
	QPushButton* paraButton;
	QPushButton* showOriginButton;
	QWidget* paraWidget;

	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;

};
