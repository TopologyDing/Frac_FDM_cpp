#pragma once
#include<iostream>
#include <iomanip>
#include<vector>

#include "Material.h"
#include "Point.h"
#include "Line.h"

#include<math.h>
#include "mkl.h"

#include  <Eigen/Sparse>
class Geometry {
public:
	Geometry();
	Geometry(double lx, double ly, int Nx, int Ny);
	~Geometry();
	void ConfigureGeometry();
	void GenerateMesh();
	void ProcessInterface();
	void ShowGeometry();
	void ConfigureGhostPoint();
	void ShowGhostPoint();
	void ShowGhostPoint1();
	void ComputeFirstOrdDeriCentral();
	void ComputeFirstOrdDeriForward();
	void ComputeFirstOrdDeriBackward();
	void ConfigureInterfacePoint(int x, int y, std::vector<std::vector<int>>& ninePointScheme);

	void NinePointInterpolation(int x, int y, double posi_x, double posi_y, 
		std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, bool dir);
	void ThreePointInterpolation(std::vector<double>& dispGridPointPortion, std::vector<double>& dispGhostPointPortion,
		std::vector<double>& deriGridPointPortion, std::vector<double>& deriGhostPointPortion, bool dir);

	void GetNeighborPointInfo(int x, int y, int parent_x, int parent_y, std::vector<int>& pointInfo);
	 
	void ComputeFracDeri(double posi_x, double posi_y, short material);
	void ComputeFirstOrdPrecisionFracDeri(int x, int y, std::vector<double>&, std::vector<double>&, bool dir);
	void ComputeSecOrdPrecisionFracDeri(int x, int y, std::vector<double>&, std::vector<double>&, bool dir);
	void GetNlcInfo(std::vector<double>&, std::vector<std::vector<int>>&, std::vector<double>&,
		int x, int y, double lengthScale, double pposi_x = -1, double pposi_y = -1);

private:
	double R;
	double lx;
	double ly;
	int Nx;
	int Ny;
	int GridSize;
	double dx;
	double dy;
	std::vector<std::vector<Point>> pointMatrix;
	// Parent point of ghost point: No repetitive elements
	std::vector<std::vector<int>> parentGhostPoint;
	//std::vector<std::vector<int>> parentGhostPointIndex;
	std::vector<GhostPoint> ghostPoints;
	std::vector<GhostPoint> ghostPoints1;
	// Ghost point index: Include diagonal ghost point
	std::vector<std::vector<int>> ghostPointIndex;
	// Ghost point coordinates: Stored according to ghost point index sequence
	std::vector<std::vector<int>> ghostPointCoord;
	// Ghost point index: Exclude diagonal ghost point
	std::vector<std::vector<int>> ghostPointIndex1;
	// Ghost point coordinates: Stored according to ghost point index sequence
	std::vector<std::vector<int>> ghostPointCoord1;
	std::vector<std::vector<double>> ghostPointGridPointPortion;
	std::vector<std::vector<double>> ghostPointGhostPointPortion;
	// First order derivative: Central (x direction and y direction)
	std::vector<std::vector<double>> firstOrdDeriCentralGridPointPortion_x;
	std::vector<std::vector<double>> firstOrdDeriCentralGhostPointPortion_x;
	std::vector<std::vector<double>> firstOrdDeriCentralGridPointPortion_y;
	std::vector<std::vector<double>> firstOrdDeriCentralGhostPointPortion_y;
	// First order derivative: Forward
	std::vector<std::vector<double>> firstOrdDeriForwardGridPointPortion;
	std::vector<std::vector<double>> firstOrdDeriForwardGhostPointPortion;
	// First order derivative: Backward
	std::vector<std::vector<double>> firstOrdDeriBackwardGridPointPortion;
	std::vector<std::vector<double>> firstOrdDeriBackwardGhostPointPortion;
	// Interface point displacement
	std::vector<std::vector<double>> interfacePointDispGridPointPortion;
	std::vector<std::vector<double>> interfacePointDispGhostPointPortion;
	// Interface point derivative
	std::vector<std::vector<double>> interfacePointDeriGridPointPortion_x;
	std::vector<std::vector<double>> interfacePointDeriGhostPointPortion_x;
	std::vector<std::vector<double>> interfacePointDeriGridPointPortion_y;
	std::vector<std::vector<double>> interfacePointDeriGhostPointPortion_y;
	// Fractional derivative
	std::vector<std::vector<double>> fracOrdDeriGridPointPortion_x;
	std::vector<std::vector<double>> fracOrdDeriGhostPointPortion_x;
	std::vector<std::vector<double>> fracOrdDeriGridPointPortion_y;
	std::vector<std::vector<double>> fracOrdDeriGhostPointPortion_y;

	std::vector<Material> materials;
	std::vector<Line> curves;
};