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
	showOriginButton = new QPushButton(tr("Show Original Mesh"));
	connect(showOriginButton, SIGNAL(clicked()), SIGNAL(ShowOriginSignal()));
	auto iterLabel = new QLabel(tr("Maximum iterations:"));
	iterLE = new QLineEdit();
	iterLE->setPlaceholderText(tr("1000"));
	auto energyLabel = new QLabel(tr("Minimum energy vatiation:"));
	energyLE = new QLineEdit();
	energyLE->setPlaceholderText(tr("0.000001"));
	paraButton = new QPushButton(tr("Parameterization"));
	connect(paraButton, SIGNAL(clicked()), this, SLOT(ParaBtnClicked()));
	

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(showOriginButton);
	layout->addWidget(iterLabel);
	layout->addWidget(iterLE);
	layout->addWidget(energyLabel);
	layout->addWidget(energyLE);
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

void MeshParamWidget::ParaBtnClicked(void)
{
	try
	{
		QString maxIterStr = iterLE->text();
		QString energyStr = energyLE->text();
		int maxIter = 1000;
		double minEnergyVar = 0.000001;
		if (!maxIterStr.isEmpty())
		{
			maxIter = maxIterStr.toInt();
		}
		if (!energyStr.isEmpty())
		{
			minEnergyVar = energyStr.toDouble();
		}
		emit ParaSignal(maxIter, minEnergyVar);
	}
	catch (const std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}
