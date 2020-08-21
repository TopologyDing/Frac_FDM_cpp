#pragma once
#include "Line.h"

class Domain {
public:
	Domain();
	~Domain();

	std::vector<int> lineIndex;
	int DomainIndex;
	int materialIndex;
};
