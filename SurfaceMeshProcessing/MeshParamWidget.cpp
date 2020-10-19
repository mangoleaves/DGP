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
	QLabel* lengthLbl = new QLabel(tr("Target Length"));
	lengthLE = new QLineEdit();
	QDoubleValidator* lengthValidator = new QDoubleValidator(0.0, INFINITY, 6);
	lengthLE->setValidator(lengthValidator);
	lengthLE->setPlaceholderText(tr("> 0.0"));
	lengthLbl->setBuddy(lengthLE);
	remeshingBtn = new QPushButton(tr("Remeshing"));
	connect(remeshingBtn, SIGNAL(clicked()), this, SLOT(EmitRemeshingSignal()));
	showOriginBtn = new QPushButton(tr("Show Origin"));
	connect(showOriginBtn, SIGNAL(clicked()), SIGNAL(ShowOriginSignal()));

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(lengthLbl);
	layout->addWidget(lengthLE);
	layout->addWidget(remeshingBtn);
	layout->addWidget(showOriginBtn);
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

void MeshParamWidget::EmitRemeshingSignal(void)
{
	if (lengthLE->hasAcceptableInput())
	{
		double targetLength = lengthLE->text().toDouble();
		emit RemeshingSignal(targetLength);
	}
}
