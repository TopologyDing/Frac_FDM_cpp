#pragma once
#include "Point.h"

enum class LineType {
	InterfaceStraightLine,
	InterfaceCircle,
	BoundaryStraightLine,
	BoundaryCircle
};


class Line {
public:
	Line();
	~Line();
	// vector stores vertex points which connect with other lines
	std::vector<Point> VertexPoints;
	// vector stores the other corresponding line index of vertex points
	std::vector<int> VertexLineIndex;
	// vector stores domainIndex which current line belongs to
	std::vector<int> DomainIndex;
	LineType lineType;
	int lineIndex;
};
