#pragma once
#include<cstring>
#include<vector>
class Point {
public:
	Point();
	Point(double px, double py);
	Point(const Point& p);
	~Point();

	double x;
	double y;
	short offset_x;
	short offset_y;
	short point_material;
	/*
	Point type: 0(normal point)
	Point type: -1(1 diagonal fictitious point)
	Point type: 1(1 fictitious point)
	Point type: 2(2 fictitious point, 1 diagonal fictitious point)
	*/
    std::vector<short> point_type;
	/*
	Interface information index:
	Right:0; Top:1; Left:2; Bottom:3;
	*/
	std::vector<double>point_interface_posi_x;
	std::vector<double>point_interface_posi_y;
};

class GhostPoint {
public:
	GhostPoint();
	GhostPoint(double px, double py, int pindex_x, int pindex_y, short pori);
	GhostPoint(const GhostPoint& p);
	~GhostPoint();
	
	double x;
	double y;
	int index_x;
	int index_y;
	/*
	Orientation of ghost point to parent point:
	Right:0; Top:1; Left:2; Bottom:3;
	*/
	short orientation;
};