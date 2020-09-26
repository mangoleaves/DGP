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
	initMorphing = new QPushButton(tr("Init"));
	doMorphing = new QPushButton(tr("Morphing"));

	connect(loadSource, SIGNAL(clicked()), SIGNAL(LoadSourceSignal()));
	connect(showSource, SIGNAL(clicked()), SIGNAL(ShowSourceSignal()));
	connect(loadTarget, SIGNAL(clicked()), SIGNAL(LoadTargetSignal()));
	connect(showTarget, SIGNAL(clicked()), SIGNAL(ShowTargetSignal()));
	connect(initMorphing, SIGNAL(clicked()), this, SLOT(EmitInitSignal()));
	connect(doMorphing, SIGNAL(clicked()), SIGNAL(DoMorphingSignal()));

	QLabel* intervalLbl = new QLabel("Interval:");
	intervalLE = new QLineEdit();
	QDoubleValidator* dv = new QDoubleValidator(0.0, 1.0, 3, this);
	intervalLE->setPlaceholderText("0.0~1.0, default 0.1");
	intervalLE->setValidator(dv);
	
	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(loadSource);
	layout->addWidget(showSource);
	layout->addWidget(loadTarget);
	layout->addWidget(showTarget);
	layout->addWidget(intervalLbl);
	layout->addWidget(intervalLE);
	layout->addWidget(initMorphing);
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

void MeshParamWidget::EmitInitSignal(void)
{
	QString intervalStr = intervalLE->text();
	if (intervalStr.isEmpty())
	{
		emit InitMorphingSignal(0.1);
	}
	else
	{
		double interval = intervalStr.toDouble();
		emit InitMorphingSignal(interval);
	}
}
