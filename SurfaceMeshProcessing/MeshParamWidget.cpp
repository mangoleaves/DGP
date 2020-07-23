#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateParaWidget();
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

void MeshParamWidget::CreateParaWidget(void)
{
	paraButton = new QPushButton(tr("Parameterization"));
	connect(paraButton, SIGNAL(clicked()), SIGNAL(ParaSignal()));
	showOriginButton = new QPushButton(tr("Show Original Mesh"));
	connect(showOriginButton, SIGNAL(clicked()), SIGNAL(ShowOriginSignal()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(showOriginButton);
	layout->addWidget(paraButton);
	paraWidget = new QWidget();
	paraWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(paraWidget, 0, 0, 1, 1);
	layout->addWidget(twParam, 1, 0, 1, 1);
	this->setLayout(layout);
}
