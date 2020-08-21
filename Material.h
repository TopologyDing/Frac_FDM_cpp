#pragma once

class Material {
public:
	Material();
	Material(double pmu, double plambda, double pv, short pmaterial_ID, double pnlc_range, double palpha);
	~Material();

	double mu; // Lame parameter mu
	double lambda; // Lame parameter lambda
	double v; // Poisson ratio
	double alpha;
	double nlc_range;
	short material_ID; //ID of the material
};