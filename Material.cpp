#include "Material.h"

Material::Material() {
	this->mu = 0;
	this->lambda = 0;
	this->v = 0;
	this->material_ID = -1;
	this->nlc_range = 0;
	this->alpha = 1;
}

Material::Material(double pmu, double plambda, double pv, short pmaterial_ID, double pnlc_range, double palpha) {
	this->mu = pmu;
	this->lambda = plambda;
	this->v = pv;
	this->material_ID = pmaterial_ID;
	this->nlc_range = pnlc_range;
	this->alpha = palpha;
}

Material::~Material() {

}



