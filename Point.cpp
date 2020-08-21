#include "Point.h"

Point::Point() {
	this->x = -1;
	this->y = -1;
	this->point_type = { -1,-1 };
	this->point_interface_posi_x = { -1,-1,-1,-1 };
	this->point_interface_posi_y = { -1,-1,-1,-1 };
	this->point_material = -1;
}

Point::Point(double px, double py) {
	this->x = px;
	this->y = py;
	this->point_type = { -1,-1 };
	this->point_interface_posi_x = { -1,-1,-1,-1 };
	this->point_interface_posi_y = { -1,-1,-1,-1 };
	this->point_material = -1;
}

Point::Point(const Point& p) {
	this->x = p.x;
	this->y = p.y;
	this->offset_x = p.offset_x;
	this->offset_y = p.offset_y;
	this->point_type = p.point_type;
	this->point_interface_posi_x = p.point_interface_posi_x;
	this->point_interface_posi_y = p.point_interface_posi_y;
	this->point_material = p.point_material;
}

Point::~Point() {
}

GhostPoint::GhostPoint() {
	this->x = -1;
	this->y = -1;
	this->index_x = -1;
	this->index_y = -1;
	this->orientation = -1;
}

GhostPoint::GhostPoint(double px, double py, int pindex_x, int pindex_y, short pori) {
	this->x = px;
	this->y = py;
	this->index_x = pindex_x;
	this->index_y = pindex_y;
	this->orientation = pori;
}

GhostPoint::GhostPoint(const GhostPoint& p) {
	this->x = p.x;
	this->y = p.y;
	this->index_x = p.index_x;
	this->index_y = p.index_y;
	this->orientation = p.orientation;
}

GhostPoint::~GhostPoint() {

}
