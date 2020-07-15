#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateCurvatureWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
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

void MeshParamWidget::CreateCurvatureWidget(void)
{
	MCButton = new QPushButton(tr("Mean Curvature"));
	AMCButton = new QPushButton(tr("Absolute Mean Curvature"));
	GCButton = new QPushButton(tr("Gaussian Curvature"));
	connect(MCButton, SIGNAL(clicked()), SIGNAL(ShowMC()));
	connect(AMCButton, SIGNAL(clicked()), SIGNAL(ShowAMC()));
	connect(GCButton, SIGNAL(clicked()), SIGNAL(ShowGC()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(MCButton);
	layout->addWidget(AMCButton);
	layout->addWidget(GCButton);
	
	curvatureWidget = new QWidget();
	curvatureWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(curvatureWidget, 1, 0, 1, 1);
	this->setLayout(layout);
}
