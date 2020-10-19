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
	void RemeshingSignal(double targetLength);
	void ShowOriginSignal();
private slots:
	void EmitRemeshingSignal(void);
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QLineEdit* lengthLE;
	QPushButton* remeshingBtn;
	QPushButton* showOriginBtn;
};
