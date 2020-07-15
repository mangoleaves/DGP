#pragma once
#include "MeshDefinition.h"

class ColorMap
{
private:
	double maxValue, minValue;
public:
	ColorMap();
	ColorMap(double maxV, double minV);

	double GetMaxValue();
	double GetMinValue();
	void SetMaxMinValue(double maxV, double minV);
	
	OpenMesh::Vec3d MapToColor(double value);
};

