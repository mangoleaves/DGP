#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));
	triBtn = new QPushButton(tr("Delaunay Triangulation"));
	connect(triBtn, SIGNAL(clicked()), SIGNAL(TriSignal()));

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(triBtn);
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}
