#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateMorphWidget();
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

void MeshParamWidget::CreateMorphWidget(void)
{
	loadSource = new QPushButton(tr("Load Source Mesh"));
	showSource = new QPushButton(tr("Show Source Mesh"));
	loadTarget = new QPushButton(tr("Load Target Mesh"));
	showTarget = new QPushButton(tr("Show Target Mesh"));
	doMorphing = new QPushButton(tr("Morphing"));
	connect(loadSource, SIGNAL(clicked()), SIGNAL(LoadSourceSignal()));
	connect(showSource, SIGNAL(clicked()), SIGNAL(ShowSourceSignal()));
	connect(loadTarget, SIGNAL(clicked()), SIGNAL(LoadTargetSignal()));
	connect(showTarget, SIGNAL(clicked()), SIGNAL(ShowTargetSignal()));
	connect(doMorphing, SIGNAL(clicked()), SIGNAL(DoMorphingSignal()));
	
	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(loadSource);
	layout->addWidget(showSource);
	layout->addWidget(loadTarget);
	layout->addWidget(showTarget);
	layout->addWidget(doMorphing);

	morphWidget = new QWidget();
	morphWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(morphWidget, 1, 0, 1, 1);
	this->setLayout(layout);
}
