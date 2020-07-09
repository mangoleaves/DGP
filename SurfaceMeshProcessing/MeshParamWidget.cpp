#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateFindShortestPathWidget();
	CreateMSTWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

QString MeshParamWidget::getBeginIdx(void)
{
	return bpiLineEdit->text();
}

QString MeshParamWidget::getEndIdx(void)
{
	return epiLineEdit->text();
}

QString MeshParamWidget::getVertexIdxs(void)
{
	return vertexIdxsLineEdit->text();
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateFindShortestPathWidget(void)
{
	fspWidget = new QWidget();
	fspTitleLabel = new QLabel(tr("Find Shortest Path"));
	fspTitleLabel->setFont(QFont("Microsoft YaHei", 10, 75));
	beginLabel = new QLabel(tr("Begin Point Index:"));
	endLabel = new QLabel(tr("End Point Index:"));
	bpiLineEdit = new QLineEdit();
	epiLineEdit = new QLineEdit();
	fspButton = new QPushButton(tr("Find"));
	connect(fspButton, SIGNAL(clicked()), SIGNAL(FindShortestPathSignal()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(fspTitleLabel);
	layout->addWidget(beginLabel);
	layout->addWidget(bpiLineEdit);
	layout->addWidget(endLabel);
	layout->addWidget(epiLineEdit);
	layout->addWidget(fspButton);
	fspWidget->setLayout(layout);
}

void MeshParamWidget::CreateMSTWidget(void)
{
	mstWidget = new QWidget();
	mstTitleLabel = new QLabel(tr("Minimal Spanning Tree"));
	mstTitleLabel->setFont(QFont("Microsoft YaHei", 10, 75));
	vertexIdxsLabel = new QLabel(tr("Vertex indexs:"));
	vertexIdxsLineEdit = new QLineEdit();
	mstButton = new QPushButton(tr("Find"));
	connect(mstButton, SIGNAL(clicked()), SIGNAL(FindMSTSignal()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(mstTitleLabel);
	layout->addWidget(vertexIdxsLabel);
	layout->addWidget(vertexIdxsLineEdit);
	layout->addWidget(mstButton);
	mstWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(fspWidget, 1, 0, 1, 1);
	layout->addWidget(mstWidget, 2, 0, 1, 1);
	this->setLayout(layout);
}
