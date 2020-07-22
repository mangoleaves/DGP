#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateCurvatureWidget();
	CreateDenoiseWidget();
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
	auto titleLbl = new QLabel(tr("Curvature"));
	titleLbl->setFont(QFont("Microsoft YaHei", 10, 75));
	MCButton = new QPushButton(tr("Mean Curvature"));
	AMCButton = new QPushButton(tr("Absolute Mean Curvature"));
	GCButton = new QPushButton(tr("Gaussian Curvature"));
	connect(MCButton, SIGNAL(clicked()), SIGNAL(ShowMC()));
	connect(AMCButton, SIGNAL(clicked()), SIGNAL(ShowAMC()));
	connect(GCButton, SIGNAL(clicked()), SIGNAL(ShowGC()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(titleLbl);
	layout->addWidget(MCButton);
	layout->addWidget(AMCButton);
	layout->addWidget(GCButton);
	
	curvatureWidget = new QWidget();
	curvatureWidget->setLayout(layout);
}

void MeshParamWidget::CreateDenoiseWidget(void)
{
	QLabel* titleLbl = new QLabel(tr("Denoise"));
	titleLbl->setFont(QFont("Microsoft YaHei", 10, 75));
	QLabel* noiseProportionLbl = new QLabel(tr("Noise Proportion"));
	QLabel* sigmaSLbl = new QLabel(tr("sigma S"));
	QLabel* nMILbl = new QLabel(tr("Normal Updating Iterations"));
	QLabel* vMILbl = new QLabel(tr("Vertex Updating Iterations"));
	noiseProportionLE = new QLineEdit();
	sigmaSLE = new QLineEdit();
	nMILE = new QLineEdit();
	vMILE = new QLineEdit();
	showOriginButton = new QPushButton(tr("Show Original Mesh"));
	addNoiseButton = new QPushButton(tr("Add Noise / Show Noisy Mesh"));
	denoiseButton = new QPushButton(tr("Denoise / Show Denoised Mesh"));
	connect(showOriginButton, SIGNAL(clicked()), SIGNAL(ShowOrigin()));
	connect(addNoiseButton, SIGNAL(clicked()), SIGNAL(AddNoise()));
	connect(denoiseButton, SIGNAL(clicked()), SIGNAL(Denoise()));

	QVBoxLayout* layout = new QVBoxLayout();
	layout->addWidget(titleLbl);
	layout->addWidget(showOriginButton);
	layout->addWidget(noiseProportionLbl);
	layout->addWidget(noiseProportionLE);
	layout->addWidget(addNoiseButton);
	layout->addWidget(sigmaSLbl);
	layout->addWidget(sigmaSLE);
	layout->addWidget(nMILbl);
	layout->addWidget(nMILE);
	layout->addWidget(vMILbl);
	layout->addWidget(vMILE);
	layout->addWidget(denoiseButton);

	denoiseWidget = new QWidget();
	denoiseWidget->setLayout(layout);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	layout->addWidget(curvatureWidget, 1, 0, 1, 1);
	layout->addWidget(denoiseWidget, 2, 0, 1, 1);
	this->setLayout(layout);
}
