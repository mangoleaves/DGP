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
	kLE = new QLineEdit();
	QIntValidator *kiv = new QIntValidator(1, 100);
	kLE->setValidator(kiv);
	kLE->setPlaceholderText(tr("K"));
	nIterLE = new QLineEdit();
	QIntValidator* nIteriv = new QIntValidator(0, 1000);
	nIterLE->setValidator(nIteriv);
	nIterLE->setPlaceholderText(tr("Iteration Number"));
	DoLloydBtn = new QPushButton(tr("Do Lloyd"));
	connect(DoLloydBtn, SIGNAL(clicked()), this, SLOT(EmitDoLloydSignal()));

	QVBoxLayout *layout = new QVBoxLayout();
	layout->addWidget(pbPrintInfo);
	layout->addWidget(kLE);
	layout->addWidget(nIterLE);
	layout->addWidget(DoLloydBtn);
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

void MeshParamWidget::EmitDoLloydSignal(void)
{
	if (kLE->hasAcceptableInput() && nIterLE->hasAcceptableInput())
	{
		int K = kLE->text().toInt();
		int nIter = nIterLE->text().toInt();
		emit DoLloydSignal(K, nIter);
	}
}
