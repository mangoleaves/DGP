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
	void CreateDeformWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void NoSelectSignal();
	void SelectFixedSignal();
	void SelectCustomSignal();
	void MoveSignal();
	void ClearSignal();
private slots:
	void DeformBtnClicked(int id);
	void ResetCheck();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QPushButton* noSelectBtn;
	QPushButton* selectFixedBtn;
	QPushButton* selectCustomBtn;
	QPushButton* moveBtn;
	QPushButton* clearBtn;
	QButtonGroup* deformBtnGroup;
	QWidget* deformWidget;
};
