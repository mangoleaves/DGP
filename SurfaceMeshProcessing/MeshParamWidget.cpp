#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateSimpWidget();
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

void MeshParamWidget::CreateSimpWidget(void)
{
	simpBtn = new QPushButton(tr("Simplify"));
	auto NvLbl = new QLabel("minimal vertex number:");
	NvLE = new QLineEdit();
	NvLE->setPlaceholderText("100");
	auto minCostLbl = new QLabel("minimal cost:");
	minCostLE = new QLineEdit();
	minCostLE->setPlaceholderText("100.000");
	auto NvVld = new QIntValidator(0, INT_MAX);
	auto minCostVld = new QDoubleValidator(0.0, INFINITY, 3);
	NvLE->setValidator(NvVld);
	minCostLE->setValidator(minCostVld);
	connect(simpBtn, SIGNAL(clicked()), this, SLOT(EmitSimpSignal()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(NvLbl);
	layout->addWidget(NvLE);
	layout->addWidget(minCostLbl);
	layout->addWidget(minCostLE);
	layout->addWidget(simpBtn);
	simpWidget = new QWidget();
	simpWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(simpWidget, 1, 0, 1, 1);
	this->setLayout(layout);
}

void MeshParamWidget::EmitSimpSignal(void)
{
	int Nv;
	double minCost;
	if (NvLE->text().isEmpty())
	{
		Nv = 100;
	}
	else
	{
		Nv = NvLE->text().toInt();
	}
	if (minCostLE->text().isEmpty())
	{
		minCost = 100.0;
	}
	else
	{
		minCost = minCostLE->text().toDouble();
	}
	emit SimpSignal(Nv, minCost);
}
