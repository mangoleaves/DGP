#include "surfacemeshprocessing.h"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	QSurfaceFormat format;
	format.setSamples(0);
	format.setDepthBufferSize(24);
	QSurfaceFormat::setDefaultFormat(format);
	SurfaceMeshProcessing mainWin;
	mainWin.show();
	return app.exec();
}
