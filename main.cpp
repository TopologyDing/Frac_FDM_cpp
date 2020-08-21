#include "Geometry.h"

int main() {
	Geometry myGeometry(1.5, 1.5, 15, 15);
	myGeometry.GenerateMesh();
	myGeometry.ProcessInterface();
	//myGeometry.ShowGeometry();
	//myGeometry.ShowGhostPoint();
	//myGeometry.ShowGhostPoint1();
	myGeometry.ConfigureGhostPoint();
	myGeometry.ComputeFirstOrdDeriCentral();
	return 0;
}


