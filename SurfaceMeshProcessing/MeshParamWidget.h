#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>
#include <QSizePolicy>

class MeshParamWidget : public QWidget
{
	Q_OBJECT

public:
	MeshParamWidget(QWidget *parent = 0);
	~MeshParamWidget(void);
	QString getBeginIdx(void);
	QString getEndIdx(void);
	QString getVertexIdxs(void);
private:
	void CreateTabWidget(void);
	void CreateFindShortestPathWidget(void);
	void CreateMSTWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void FindShortestPathSignal();
	void FindMSTSignal();
private:
	// Shortest path
	QWidget* fspWidget;		// Find shortest path widget
	QLabel* fspTitleLabel;
	QLabel* beginLabel;
	QLabel* endLabel;
	QLineEdit* bpiLineEdit;	// Begin point index line edit
	QLineEdit* epiLineEdit;	// End point index line edit
	QPushButton* fspButton;	// Find shortest path button
	// MST
	QWidget* mstWidget;
	QLabel* mstTitleLabel;
	QLabel* vertexIdxsLabel;
	QLineEdit* vertexIdxsLineEdit;
	QPushButton* mstButton;

	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;

};
