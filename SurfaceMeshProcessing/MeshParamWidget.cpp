#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateDeformWidget();
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

void MeshParamWidget::CreateDeformWidget(void)
{
	deformBtnGroup = new QButtonGroup();
	noSelectBtn = new QPushButton(tr("No Select"));
	selectFixedBtn = new QPushButton(tr("Select Fixed"));
	selectCustomBtn = new QPushButton(tr("Select Custom"));
	moveBtn = new QPushButton(tr("Move"));
	clearBtn = new QPushButton(tr("Clear"));

	noSelectBtn->setCheckable(true);
	noSelectBtn->setChecked(true);
	selectFixedBtn->setCheckable(true);
	selectCustomBtn->setCheckable(true);
	moveBtn->setCheckable(true);
	deformBtnGroup->addButton(noSelectBtn, 0);
	deformBtnGroup->addButton(selectFixedBtn, 1);
	deformBtnGroup->addButton(selectCustomBtn, 2);
	deformBtnGroup->addButton(moveBtn, 3);
	deformBtnGroup->setExclusive(true);

	connect(deformBtnGroup, SIGNAL(buttonClicked(int)), this, SLOT(DeformBtnClicked(int)));
	connect(clearBtn, SIGNAL(clicked()), SIGNAL(ClearSignal()));
	connect(clearBtn, SIGNAL(clicked()), SLOT(ResetCheck()));

	auto layout = new QVBoxLayout();
	layout->addWidget(noSelectBtn);
	layout->addWidget(selectFixedBtn);
	layout->addWidget(selectCustomBtn);
	layout->addWidget(moveBtn);
	layout->addWidget(clearBtn);

	deformWidget = new QWidget();
	deformWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(deformWidget, 1, 0, 1, 1);
	this->setLayout(layout);
}

void MeshParamWidget::DeformBtnClicked(int id)
{
	switch (id)
	{
	case 0:
		emit NoSelectSignal();
	case 1:
		emit SelectFixedSignal();
		break;
	case 2:
		emit SelectCustomSignal();
		break;
	case 3:
		emit MoveSignal();
		break;
	default:
		break;
	}
}

void MeshParamWidget::ResetCheck()
{
	noSelectBtn->setChecked(true);
	selectFixedBtn->setChecked(false);
	selectCustomBtn->setChecked(false);
	moveBtn->setChecked(false);
}
