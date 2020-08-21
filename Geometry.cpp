#include "Geometry.h"
Geometry::Geometry(double pLx, double pLy, int pNx, int pNy) {
	this->lx = pLx;
	this->ly = pLy;
	this->Nx = pNx;
	this->Ny = pNy;
	this->GridSize = (Nx + 1) * (Ny + 1);
	this->dx = lx / Nx;
	this->dy = ly / Ny;
	this->R = 0.875;
	const MKL_INT row = Nx + 1;
	const MKL_INT col = Ny + 1;
}

void Geometry::ConfigureGeometry() {
	double R = 0.875;
	Line line;
	line.lineType = LineType::BoundaryStraightLine;
	line.VertexPoints.push_back(Point(0,0));
	line.VertexPoints.push_back(Point(0,R));
	line.lineIndex = 1;
	curves.push_back(line);
}



void Geometry::GenerateMesh() {
	std::vector<Point> rowPoints;
	Point curPoint;
	for (int i = 0; i < Nx + 1; ++i) {
		for (int j = 0; j < Ny + 1; ++j) {
			curPoint = Point(i * dx, j * dy);
			double r = sqrt(pow(i * dx, 2) + pow(j * dy, 2));
			// Set material type
			if (r < this->R)
				curPoint.point_material = 1;
			else
				curPoint.point_material = 2;

			// Set offset (nine point scheme)
			if (i == 0) curPoint.offset_x = 1;
			else if (i == Nx) curPoint.offset_x = -1;
			else curPoint.offset_x = 0;

			if (j == 0) curPoint.offset_y = 1;
			else if (j == Ny) curPoint.offset_y = -1;
			else curPoint.offset_y = 0;

			rowPoints.push_back(curPoint);
		}
		this->pointMatrix.push_back(rowPoints);
		rowPoints.clear();
	}
}

void Geometry::ProcessInterface() {
	int a[2] = { -1, 1 };
	GhostPoint curGhostPoint;
	std::vector<int> rowGhostPointIndex;
	std::vector<int> rowGhostPointIndex1;
	std::vector<int> rowGhostPointCoord;
	std::vector<int> rowGhostPointCoord1;
	// Preallocate capacity for row vectors
	rowGhostPointIndex.reserve(1);
	rowGhostPointIndex1.reserve(1);
	rowGhostPointCoord.reserve(2);
	rowGhostPointCoord1.reserve(2);
	// Preallocate capacity for class member vectors
	this->ghostPointIndex.reserve((Nx + 1) * (Ny + 1));
	this->ghostPointIndex1.reserve((Nx + 1) * (Ny + 1));
	this->ghostPointCoord.reserve((Nx + 1) * (Ny + 1) / 4);
	this->ghostPointCoord1.reserve((Nx + 1) * (Ny + 1) / 4);
	int totalGhostPoints = -1;
	int totalGhostPoints1 = -1;

	for (int i = 0; i < Nx + 1; ++i) {
		for (int j = 0; j < Ny + 1; ++j) {
			Point curPoint = pointMatrix[i][j];
			// interface at current point's right direction
			int t = 0;
			int rowBegin = 0;
			int colBegin = 0;
			int rowSize = 2;
			int colSize = 2;
			int rowIndex = 0;
			int colIndex = 0;
			switch (curPoint.offset_x) {
			case 0:
				break;
			case -1:
				rowBegin = 0;
				rowSize = 1;
				break;
			case 1:
				rowBegin = 1;
				rowSize = 2;
				break;
			}
			switch (curPoint.offset_y) {
			case 0:
				break;
			case -1:
				colBegin = 0;
				colSize = 1;
				break;
			case 1:
				colBegin = 1;
				colSize = 2;
				break;
			}
			// Iterate in nine point scheme to detect ghost point: row direction
			for (int k = rowBegin; k < rowSize; ++k) {
				if (curPoint.point_material != pointMatrix[i + a[k]][j].point_material) {
					t = t + 1;
					int r = (a[k] == 1) ? 1 : 3;
					pointMatrix[i][j].point_type[1] = r;
					pointMatrix[i][j].point_interface_posi_x[r - 1] = sqrt(pow(R, 2) - pow(j * dy, 2));
					pointMatrix[i][j].point_interface_posi_y[r - 1] = j * dy;
					curGhostPoint = GhostPoint((i + a[k]) * dx, j * dy, i, j, r);
					ghostPoints.push_back(curGhostPoint);
					ghostPoints1.push_back(curGhostPoint);
					rowIndex = a[k];
					break;
				}
			}
			// Iterate in nine point scheme to detect ghost point: column direction
			for (int s = colBegin; s < colSize; ++s) {
				if (curPoint.point_material != pointMatrix[i][j + a[s]].point_material) {
					t = t + 1;
					int r = (a[s] == 1) ? 2 : 4;
					pointMatrix[i][j].point_type[1] = r;
					pointMatrix[i][j].point_interface_posi_x[r - 1] = i * dx;
					pointMatrix[i][j].point_interface_posi_y[r - 1] = sqrt(pow(R, 2) - pow(i * dx, 2));
					curGhostPoint = GhostPoint(i * dx, (j + a[s]) * dy, i, j, r);
					ghostPoints.push_back(curGhostPoint);
					ghostPoints1.push_back(curGhostPoint);
					colIndex = a[s];
					break;
				}
			}
			// Iterate in nine point scheme to detect ghost point: diagonal direction
			switch (t) {
			case 0:
				for (int k = rowBegin; k < rowSize; ++k) {
					for (int s = colBegin; s < colSize; ++s) {
						if (curPoint.point_material != pointMatrix[i + a[k]][j + a[s]].point_material) {
							t = -1;
							int r;
							// Get subtype of type -1 ghost point
							switch (a[k]) {
							case -1:
								switch (a[s]) {
								case -1:
									r = 3;
									break;
								case 1:
									r = 2;
									break;
								}
							case 1:
								switch (a[s]) {
								case -1:
									r = 4;
									break;
								case 1:
									r = 1;
									break;
								}
							}
							curGhostPoint = GhostPoint((i + a[k]) * dx, (j + a[s]) * dy, i, j, r);
							ghostPoints.push_back(curGhostPoint);
							pointMatrix[i][j].point_type[1] = curGhostPoint.orientation;
							break;
						}
					}
				}
				break;
			case 1:
				pointMatrix[i][j].point_type[1] = curGhostPoint.orientation;
				break;
			case 2:
				// Get subtype of type 2 ghost point
				int r;
				switch (rowIndex) {
				case -1:
					switch (colIndex) {
					case -1:
						r = 3;
						break;
					case 1:
						r = 2;
						break;
					}
				case 1:
					switch (colIndex) {
					case -1:
						r = 4;
						break;
					case 1:
						r = 1;
						break;
					}
				}
				curGhostPoint = GhostPoint((i + rowIndex) * dx, (j + colIndex) * dy, i, j, r);
				ghostPoints.push_back(curGhostPoint);
				pointMatrix[i][j].point_type[1] = curGhostPoint.orientation;
				break;
			}
			pointMatrix[i][j].point_type[0] = t;
			// Store ghost point information
			switch (t) {
			case -1:
				rowGhostPointIndex.push_back(totalGhostPoints + 1);
				rowGhostPointIndex1.push_back(-1);
				totalGhostPoints = totalGhostPoints + 1;
				rowGhostPointCoord.push_back(i);
				rowGhostPointCoord.push_back(j);
				ghostPointCoord.push_back(rowGhostPointCoord);
				parentGhostPoint.push_back(rowGhostPointCoord);
				rowGhostPointCoord.clear();
				break;
			case 0:
				rowGhostPointIndex.push_back(-1);
				rowGhostPointIndex1.push_back(-1);
				break;
			case 1:
				rowGhostPointIndex.push_back(totalGhostPoints + 1);
				rowGhostPointIndex1.push_back(totalGhostPoints1 + 1);
				totalGhostPoints = totalGhostPoints + 1;
				totalGhostPoints1 = totalGhostPoints1 + 1;
				rowGhostPointCoord.push_back(i);
				rowGhostPointCoord.push_back(j);
				ghostPointCoord.push_back(rowGhostPointCoord);
				ghostPointCoord1.push_back(rowGhostPointCoord);
				parentGhostPoint.push_back(rowGhostPointCoord);
				rowGhostPointCoord.clear();
				break;
			case 2:
				rowGhostPointIndex.push_back(totalGhostPoints + 1);
				rowGhostPointIndex1.push_back(totalGhostPoints1 + 1);
				totalGhostPoints = totalGhostPoints + 3;
				totalGhostPoints1 = totalGhostPoints1 + 2;
				rowGhostPointCoord.push_back(i);
				rowGhostPointCoord.push_back(j);
				ghostPointCoord.push_back(rowGhostPointCoord);
				ghostPointCoord.push_back(rowGhostPointCoord);
				ghostPointCoord.push_back(rowGhostPointCoord);
				ghostPointCoord1.push_back(rowGhostPointCoord);
				ghostPointCoord1.push_back(rowGhostPointCoord);
				parentGhostPoint.push_back(rowGhostPointCoord);
				rowGhostPointCoord.clear();
				break;
			}
		}
		ghostPointIndex.push_back(rowGhostPointIndex);
		ghostPointIndex1.push_back(rowGhostPointIndex1);
		rowGhostPointIndex.clear();
		rowGhostPointIndex1.clear();
	}
}

void Geometry::ShowGeometry() {
	std::cout << "Print geometry configuration (Number means material type) .............." << std::endl;
	for (auto rowIterator = this->pointMatrix.begin(); rowIterator != pointMatrix.end(); ++rowIterator) {
		for (auto colIterator = (*rowIterator).begin(); colIterator != (*rowIterator).end(); ++colIterator) {
			std::cout << std::setw(1) <<std::setfill('0') << (*colIterator).point_material << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------------------------------------------------------" << std::endl;
}

void Geometry::ShowGhostPoint() {
	std::cout << "Print ghost point configuration (Number means ghost point index) .............." << std::endl;
	for (auto rowIterator = this->ghostPointIndex.begin(); rowIterator != ghostPointIndex.end(); ++rowIterator) {
		for (auto colIterator = (*rowIterator).begin(); colIterator != (*rowIterator).end(); ++colIterator) {
			std::cout << std::setw(2) << std::setfill('0') << *colIterator << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------------------------------------------------------" << std::endl;

	std::cout << "Print interface (* means interface adjacent point) .............." << std::endl;
	for (auto rowIterator = this->ghostPointIndex.begin(); rowIterator != ghostPointIndex.end(); ++rowIterator) {
		for (auto colIterator = (*rowIterator).begin(); colIterator != (*rowIterator).end(); ++colIterator) {
			if (*colIterator == -1) std::cout << "o" << " ";
			else std::cout << "x" << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------------------------------------------------------" << std::endl;
}

void Geometry::ShowGhostPoint1() {
	std::cout << "Print ghost point configuration (Number means ghost point index1) .............." << std::endl;
	for (auto rowIterator = this->ghostPointIndex1.begin(); rowIterator != ghostPointIndex1.end(); ++rowIterator) {
		for (auto colIterator = (*rowIterator).begin(); colIterator != (*rowIterator).end(); ++colIterator) {
			std::cout << std::setw(2) << std::setfill('0') << *colIterator << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------------------------------------------------------" << std::endl;

	std::cout << "Print interface (* means interface adjacent point) .............." << std::endl;
	for (auto rowIterator = this->ghostPointIndex1.begin(); rowIterator != ghostPointIndex1.end(); ++rowIterator) {
		for (auto colIterator = (*rowIterator).begin(); colIterator != (*rowIterator).end(); ++colIterator) {
			if (*colIterator == -1) std::cout << "o" << " ";
			else std::cout << "x" << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "--------------------------------------------------------------------------" << std::endl;
}

Geometry::~Geometry() {

}

/*
-------------------------------------
 Ghost point1   | 0 | 1 | 2 | 3 | k |
-------------------------------------
Ghost point | 0 |   |   |   |   |   |
-------------------------------------
Ghost point | 1 |   |   |   |   |   |
-------------------------------------
Ghost point | 2 |   |   |   |   |   |
-------------------------------------
Ghost point | 3 |   |   |   |   |   |
-------------------------------------
Ghost point | i |   |   |   |   |   |
-------------------------------------
*/
void Geometry::ConfigureGhostPoint() {
	int x = -1;
	int y = -1;
	int gridIndex = 0;
	int rowIndex = 0;
	int colIndex = 0;
	int counts = 0;
	std::vector<double> rowghostPointGhostPointPortion;
	std::vector<double> rowghostPointGridPointPortion;
	//Preallocate memory for vectors
	rowghostPointGhostPointPortion.reserve(this->ghostPointCoord1.size() * 2);
	rowghostPointGridPointPortion.reserve((Nx + 1) * (Ny + 1) * 2);
	this->ghostPointGhostPointPortion.reserve(this->ghostPointCoord.size());
	this->ghostPointGridPointPortion.reserve(this->ghostPointCoord.size());

	for (int i = 0; i != this->ghostPointCoord.size(); ++i) {
		if (counts == 3) {
			i += 2;
			continue;
		}
		//Reset values of row vector as zero
		//TODO: More efficient method can be used here in the future
		rowghostPointGhostPointPortion = std::vector<double>(this->ghostPointCoord1.size() * 2, 0);
		rowghostPointGridPointPortion = std::vector<double>((Nx + 1) * (Ny + 1) * 2, 0);
		x = ghostPointCoord[i][0];
		y = ghostPointCoord[i][1];
		int k = 0;
		switch (((this->pointMatrix)[x][y]).point_type[0]) {
		case -1:
			counts = 0;
			// Get current ghost point index
			switch (this->ghostPoints[this->ghostPointIndex[x][y]].orientation) {
			case 1:
				rowIndex = 1;
				colIndex = 1;
				break;
			case 2:
				rowIndex = -1;
				colIndex = 1;
				break;
			case 3:
				rowIndex = -1;
				colIndex = -1;
				break;
			case 4:
				rowIndex = 1;
				colIndex = -1;
				break;
			}
			// 2D Taylor series approximation
			// (0,0)
			rowghostPointGridPointPortion[2 * (Ny * x + y)] = -4 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * x + y) + 1] = -4 / 3.0;
			// (rowIndex,0)
			rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + y)] = 4 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + y) + 1] = 4 / 3.0;
			// (0,colIndex)
			rowghostPointGridPointPortion[2 * (Ny * x + (y + colIndex))] = 4 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * x + (y + colIndex))+1] = 4 / 3.0;
			// (-rowIndex,colIndex)
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y + colIndex))] = -1 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y + colIndex)) + 1] = -1 / 3.0;
			// (rowIndex,-colIndex)
			rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + (y - colIndex))] = -1 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + (y - colIndex)) + 1] = -1 / 3.0;
			// (-rowIndex,-colIndex)
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y - colIndex))] = 1 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y - colIndex)) + 1] = 1 / 3.0;
			break;
		case 0:
			counts = 0;
			break;
		case 1:
			counts = 0;
			// Get index of current ghost point (exclude diagonal point)
			// index of current ghost point (start from 0)
			k = this->ghostPointIndex1[x][y];
			// Set ghost point portion of current ghost point
			rowghostPointGhostPointPortion[2 * k] = 1;
			rowghostPointGhostPointPortion[2 * k + 1] = 1;
			break;
		case 2:
			// Get index of current ghost point (exclude diagonal point)
		    k = this->ghostPointIndex1[x][y];
			// Set ghost point portion of current ghost point (row direction)
			rowghostPointGhostPointPortion[2 * k] = 1;
			rowghostPointGhostPointPortion[2 * k + 1] = 1;
			this->ghostPointGhostPointPortion.push_back(rowghostPointGhostPointPortion);
			this->ghostPointGridPointPortion.push_back(rowghostPointGridPointPortion);
			rowghostPointGhostPointPortion = std::vector<double>(this->ghostPointCoord1.size() * 2, 0);
			counts += 1;
			// Set ghost point portion of current ghost point (column direction)
			rowghostPointGhostPointPortion[2 * (k + 1)] = 1;
			rowghostPointGhostPointPortion[2 * (k + 1) + 1] = 1;
			this->ghostPointGhostPointPortion.push_back(rowghostPointGhostPointPortion);
			this->ghostPointGridPointPortion.push_back(rowghostPointGridPointPortion);
			rowghostPointGhostPointPortion = std::vector<double>(this->ghostPointCoord1.size() * 2, 0);
			counts += 1;
			// Set ghost and grid point portion of current ghost point (diagonal direction)
			switch (this->ghostPoints[this->ghostPointIndex[x][y] + 2].orientation) {
			case 1:
				rowIndex = 1;
				colIndex = 1;
				break;
			case 2:
				rowIndex = -1;
				colIndex = 1;
				break;
			case 3:
				rowIndex = -1;
				colIndex = -1;
				break;
			case 4:
				rowIndex = 1;
				colIndex = -1;
				break;
			}
			// 2D Taylor series approximation
			// (0,0)
			rowghostPointGridPointPortion[2 * (Ny * x + y)] = -4 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * x + y) + 1] = -4 / 3.0;
			// (rowIndex,0) ghost point
			rowghostPointGhostPointPortion[2 * k] = 4 / 3.0;
			rowghostPointGhostPointPortion[2 * k + 1] = 4 / 3.0;
			// (0,colIndex) ghost point
			rowghostPointGhostPointPortion[2 * (k + 1)] = 4 / 3.0;
			rowghostPointGhostPointPortion[2 * (k + 1) + 1] = 4 / 3.0;

			// (-rowIndex,colIndex)
			if (this->pointMatrix[x][y].point_material != this->pointMatrix[x - rowIndex][y + colIndex].point_material) {
				// ghost point case
				int s = this->ghostPointIndex1[x-rowIndex][y];
				switch (this->pointMatrix[x - rowIndex][y].point_type[0]) {
				case 1:
					// At column direction
					rowghostPointGhostPointPortion[2 * s] = -1 / 3.0;
					rowghostPointGhostPointPortion[2 * s + 1] = -1 / 3.0;
				case 2:
					// At column direction
					rowghostPointGhostPointPortion[2 * (s + 1)] = -1 / 3.0;
					rowghostPointGhostPointPortion[2 * (s + 1) + 1] = -1 / 3.0;
				}
			}
			else {
				// not ghost point case
				rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y + colIndex))] = -1 / 3.0;
				rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y + colIndex)) + 1] = -1 / 3.0;
			}

			// (rowIndex,-colIndex)
			if (this->pointMatrix[x][y].point_material != this->pointMatrix[x + rowIndex][y - colIndex].point_material) {
				// ghost point case
				int s = this->ghostPointIndex1[x - rowIndex][y];
				switch (this->pointMatrix[x - rowIndex][y].point_type[0]) {
				case 1:
					// At row direction
					rowghostPointGhostPointPortion[2 * s] = -1 / 3.0;
					rowghostPointGhostPointPortion[2 * s + 1] = -1 / 3.0;
				case 2:
					// At row direction
					rowghostPointGhostPointPortion[2 * s] = -1 / 3.0;
					rowghostPointGhostPointPortion[2 * s + 1] = -1 / 3.0;
				}
			}
			else {
				// not ghost point case
				rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + (y - colIndex))] = -1 / 3.0;
				rowghostPointGridPointPortion[2 * (Ny * (x + rowIndex) + (y - colIndex)) + 1] = -1 / 3.0;
			}

			// (-rowIndex,-colIndex)
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y - colIndex))] = 1 / 3.0;
			rowghostPointGridPointPortion[2 * (Ny * (x - rowIndex) + (y - colIndex)) + 1] = 1 / 3.0;
			break;
		}
		this->ghostPointGhostPointPortion.push_back(rowghostPointGhostPointPortion);
		this->ghostPointGridPointPortion.push_back(rowghostPointGridPointPortion);
		rowghostPointGhostPointPortion.clear();
		rowghostPointGridPointPortion.clear();
	}
}

void Geometry::ComputeFirstOrdDeriCentral() {
	int pointIndex = 0;
	int k = 0;
	std::vector<double> rowfirstOrdDeriCentralGridPointPortion_x;
	std::vector<double> rowfirstOrdDeriCentralGridPointPortion_y;
	std::vector<double> rowfirstOrdDeriCentralGhostPointPortion_x;
	std::vector<double> rowfirstOrdDeriCentralGhostPointPortion_y;
	// inner loop: x direction
	// outer loop: y direction
	// for point (x, y), use index = (Ny + 1) * x + y
	for (int i = 0; i < Nx + 1; ++i) {
		for (int j = 0; j < Ny + 1; ++j) {
			rowfirstOrdDeriCentralGridPointPortion_x = std::vector<double>((Nx + 1) * (Ny + 1) * 2, 0);
			rowfirstOrdDeriCentralGridPointPortion_y = std::vector<double>((Nx + 1) * (Ny + 1) * 2, 0);
			rowfirstOrdDeriCentralGhostPointPortion_x = std::vector<double>(this->ghostPointCoord1.size() * 2, 0);
			rowfirstOrdDeriCentralGhostPointPortion_y = std::vector<double>(this->ghostPointCoord1.size() * 2, 0);
			// Compute derivatives for type -1 and 0 points
			if (pointMatrix[i][j].point_type[0] == -1 || pointMatrix[i][j].point_type[0] == 0) {
				switch (pointMatrix[i][j].offset_x) {
				case 1:// Left boundary case
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j)] = -3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j) + 1] = -3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 2) + j)] = -1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 2) + j) + 1] = -1 / (2.0 * dx);
					firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
					firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
					rowfirstOrdDeriCentralGridPointPortion_x.clear();
					rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					break;
				case 0:// Center case
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -1 / (2.0 * dx);
					firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
					firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
					rowfirstOrdDeriCentralGridPointPortion_x.clear();
					rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					break;
				case -1:// Right boundary case
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j)] = 3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j) + 1] = 3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 2) + j)] = 1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 2) + j) + 1] = 1 / (2.0 * dx);
					firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
					firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
					rowfirstOrdDeriCentralGridPointPortion_x.clear();
					rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					break;
				}
				switch (pointMatrix[i][j].offset_y) {
				case 1:// Top boundary case
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j)] = -3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j) + 1] = -3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 2)] = -1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 2) + 1] = -1 / (2.0 * dy);
					firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
					firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
					rowfirstOrdDeriCentralGridPointPortion_y.clear();
					rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					break;
				case 0:
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -1 / (2.0 * dy);
					firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
					firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
					rowfirstOrdDeriCentralGridPointPortion_y.clear();
					rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					break;
				case -1:
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j)] = 3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j) + 1] = 3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 2)] = 1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 2) + 1] = 1 / (2.0 * dy);
					firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
					firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
					rowfirstOrdDeriCentralGridPointPortion_y.clear();
					rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					break;
				}
			}
			// Compute derivatives for type 1 points
			else if (pointMatrix[i][j].point_type[0] == 1) {
				switch (pointMatrix[i][j].offset_x) {// Compute x direction derivative
				case 1:// Left boundary case
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j)] = -3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j) + 1] = -3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 2) + j)] = -1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 2) + j) + 1] = -1 / (2.0 * dx);
					firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
					firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
					rowfirstOrdDeriCentralGridPointPortion_x.clear();
					rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					break;
				case 0:// Center case
					switch (pointMatrix[i][j].point_type[1]) {
					case 1:
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_x.push_back((*rowIterator) / (2.0 * dx));
						}
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -1 / (2.0 * dx);
						firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
						firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
						rowfirstOrdDeriCentralGridPointPortion_x.clear();
						rowfirstOrdDeriCentralGhostPointPortion_x.clear();
						break;
					case 3:
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_x.push_back(-(*rowIterator) / (2.0 * dx));
						}
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 1 / (2.0 * dx);
						firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
						firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
						rowfirstOrdDeriCentralGridPointPortion_x.clear();
						rowfirstOrdDeriCentralGhostPointPortion_x.clear();
						break;
					default:// No ghost point
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 1 / (2.0 * dx);
						firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
						firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
						rowfirstOrdDeriCentralGridPointPortion_x.clear();
						rowfirstOrdDeriCentralGhostPointPortion_x.clear();
						break;
					}
					break;
				case -1:// Right boundary case
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j)] = 3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * i + j) + 1] = 3 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -4 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 2) + j)] = 1 / (2.0 * dx);
					rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 2) + j) + 1] = 1 / (2.0 * dx);
					firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
					firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
					rowfirstOrdDeriCentralGridPointPortion_x.clear();
					rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					break;
				}
				switch (pointMatrix[i][j].offset_y) {// Compute y direction derivative
				case 1:// Top boundary case
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j)] = -3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j) + 1] = -3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 2)] = -1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 2) + 1] = -1 / (2.0 * dy);
					firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
					firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
					rowfirstOrdDeriCentralGridPointPortion_y.clear();
					rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					break;
				case 0:// Center case
					switch (pointMatrix[i][j].point_type[1]) {
					case 2:
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_y.push_back((*rowIterator) / (2.0 * dy));
						}
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -1 / (2.0 * dy);
						firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
						firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
						rowfirstOrdDeriCentralGridPointPortion_y.clear();
						rowfirstOrdDeriCentralGhostPointPortion_y.clear();
						break;
					case 4:
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_y.push_back(-(*rowIterator) / (2.0 * dy));
						}
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 1 / (2.0 * dy);
						firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
						firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
						rowfirstOrdDeriCentralGridPointPortion_y.clear();
						rowfirstOrdDeriCentralGhostPointPortion_y.clear();
						break;
					default:// No ghost point
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 1 / (2.0 * dy);
						firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
						firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
						rowfirstOrdDeriCentralGridPointPortion_y.clear();
						rowfirstOrdDeriCentralGhostPointPortion_y.clear();
						break;
					}
					break;
				case -1:// Bottom boundary case
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j)] = 3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j) + 1] = 3 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -4 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 2)] = 1 / (2.0 * dy);
					rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 2) + 1] = 1 / (2.0 * dy);
					firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
					firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
					rowfirstOrdDeriCentralGridPointPortion_y.clear();
					rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					break;
				}
			}
			// Compute derivative for type 2 points
			else if (pointMatrix[i][j].point_type[0] == 2) {
				// If pointMatrix[i][j].point_type[0] == 2, then normally offset_x == 0 and offset_y == 0	
				switch (pointMatrix[i][j].offset_x) {// Compute x direction derivative
				case -1:
					break;
				case 0:
					// Case when ghost point is at right direction
					if (pointMatrix[i][j].point_type[1] == 1 || pointMatrix[i][j].point_type[1] == 4) {
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_x.push_back((*rowIterator) / (2.0 * dx));
						}
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j)] = -1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i - 1) + j) + 1] = -1 / (2.0 * dx);
						firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
						firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
						rowfirstOrdDeriCentralGridPointPortion_x.clear();
						rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					}
					// Case when ghost point is at left direction
					else if (pointMatrix[i][j].point_type[1] == 2 || pointMatrix[i][j].point_type[1] == 3) {
						k = this->ghostPointIndex[i][j];
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_x.push_back(-(*rowIterator) / (2.0 * dx));
						}
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j)] = 1 / (2.0 * dx);
						rowfirstOrdDeriCentralGridPointPortion_x[2 * ((Ny + 1) * (i + 1) + j) + 1] = 1 / (2.0 * dx);
						firstOrdDeriCentralGridPointPortion_x.push_back(rowfirstOrdDeriCentralGridPointPortion_x);
						firstOrdDeriCentralGhostPointPortion_x.push_back(rowfirstOrdDeriCentralGhostPointPortion_x);
						rowfirstOrdDeriCentralGridPointPortion_x.clear();
						rowfirstOrdDeriCentralGhostPointPortion_x.clear();
					}
					break;
				case 1:
					break;
				}
				switch (pointMatrix[i][j].offset_y) {// Compute y direction derivative
				case -1:
					break;
				case 0:
					// Case when ghost point is at top direction
					if (pointMatrix[i][j].point_type[1] == 1 || pointMatrix[i][j].point_type[1] == 2) {
						k = this->ghostPointIndex[i][j] + 1;
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_y.push_back((*rowIterator) / (2.0 * dy));
						}
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1)] = -1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j - 1) + 1] = -1 / (2.0 * dy);
						firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
						firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
						rowfirstOrdDeriCentralGridPointPortion_y.clear();
						rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					}
					// Case when ghost point is at bottom direction
					else if (pointMatrix[i][j].point_type[1] == 3 || pointMatrix[i][j].point_type[1] == 4) {
						k = this->ghostPointIndex[i][j] + 1;
						for (auto rowIterator = (ghostPointGhostPointPortion[k]).begin(); rowIterator != (ghostPointGhostPointPortion[k]).end(); ++rowIterator) {
							rowfirstOrdDeriCentralGhostPointPortion_y.push_back(-(*rowIterator) / (2.0 * dy));
						}
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1)] = 1 / (2.0 * dy);
						rowfirstOrdDeriCentralGridPointPortion_y[2 * ((Ny + 1) * i + j + 1) + 1] = 1 / (2.0 * dy);
						firstOrdDeriCentralGridPointPortion_y.push_back(rowfirstOrdDeriCentralGridPointPortion_y);
						firstOrdDeriCentralGhostPointPortion_y.push_back(rowfirstOrdDeriCentralGhostPointPortion_y);
						rowfirstOrdDeriCentralGridPointPortion_y.clear();
						rowfirstOrdDeriCentralGhostPointPortion_y.clear();
					}
					break;
				case 1:
					break;
				}
			}
		}
	}
}

void Geometry::ComputeFirstOrdDeriForward() {

}

void Geometry::ComputeFirstOrdDeriBackward() {

}

void Geometry::ConfigureInterfacePoint(int x, int y, std::vector<std::vector<int>>& ninePointScheme) {
	//2D array containing information of nine point scheme: 
	for (int i = 0; i < ghostPointCoord1.size(); ++i) {
		x = ghostPointCoord1[i][0];
		y = ghostPointCoord1[i][1];
	}
}

void Geometry::NinePointInterpolation(int x, int y, double posi_x, double posi_y,
	std::vector<double>& dispGridPointPortion, std::vector<double>& deriGridPointPortion_x, std::vector<double>& deriGridPointPortion_y,
	std::vector<double>& deriGhostPointPortion_x, std::vector<double>& deriGhostPointPortion_y, bool dir) {
	// vectors should be initialized before using Geometry::NinePointInterpolation()
	std::vector<std::vector<int>> ninePointScheme(9, std::vector<int>(5, 0));
	int a[9] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };
	int b[9] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };
	// get neighbor point info first
	for (int i = 0; i < 9; ++i) {
		GetNeighborPointInfo(a[i], b[i], x, y, ninePointScheme[i]);
	}
	// nine point interpolation: based on Lagrange interpolation
	// Lagrange polynomial coefficient
	double coeff_x[3];
	double coeff_y[3];
	coeff_x[0] = (posi_x - (x)*dx) * (posi_x - (x + 1.) * dx) / ((-dx) * (-2 * dx));
	coeff_x[1] = (posi_x - (x - 1.) * dx) * (posi_x - (x + 1.) * dx) / ((dx) * (-dx));
	coeff_x[2] = (posi_x - (x - 1.) * dx) * (posi_x - (x)*dx) / ((dx) * (2 * dx));
	coeff_y[0] = (posi_y - (y)*dy) * (posi_y - (y + 1) * dy) / ((-dy) * (-2 * dy));
	coeff_y[1] = (posi_y - (y - 1.) * dy) * (posi_y - (y + 1.) * dy) / ((dy) * (-dy));
	coeff_y[2] = (posi_y - (y - 1.) * dy) * (posi_y - (y)*dy) / ((dy) * (2 * dy));
	
	// vector save infomation of three interpolated points
	std::vector<std::vector<double>> rowdispGridPointPortion_x(3, std::vector<double>((Nx + 1) * (Ny + 1) * 2, 0));
	std::vector<std::vector<double>> rowdispGridPointPortion_y(3, std::vector<double>((Nx + 1) * (Ny + 1) * 2, 0));
	std::vector<std::vector<double>> rowdispGhostPointPortion_x(3, std::vector<double>(this->ghostPointCoord1.size() * 2, 0));
	std::vector<std::vector<double>> rowdispGhostPointPortion_y(3, std::vector<double>(this->ghostPointCoord1.size() * 2, 0));
	int neighborPoint_x;
	int neighborPoint_y;
	int parient_x;
	int parient_y;
	int neighborPointIndex;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			// interpolate on x direction
			neighborPoint_x = ninePointScheme[3 * i + j][0];
			neighborPoint_y = ninePointScheme[3 * i + j][1];
			parient_x = ninePointScheme[3 * i + j][2];
			parient_y = ninePointScheme[3 * i + j][3];
			neighborPointIndex = ninePointScheme[3 * i + j][4];
			// case the point is not ghost point
			if (neighborPointIndex == -1) {
				rowdispGridPointPortion_x[i][Ny * neighborPoint_x + neighborPoint_y] += coeff_x[j];
				rowdispGridPointPortion_x[i][Ny * neighborPoint_x + neighborPoint_y + 1] += coeff_x[j];
			}
			// case the point is ghost point
			else {
				// add grid point portion
				for (int k = 0; k < ghostPointGridPointPortion[neighborPointIndex].size(); ++k) {
					rowdispGridPointPortion_x[i][k] += ghostPointGridPointPortion[neighborPointIndex][k] * coeff_y[j];
				}
				// add ghost point portion
				for (int k = 0; k < ghostPointGhostPointPortion[neighborPointIndex].size(); ++k) {
					rowdispGridPointPortion_x[i][k] += ghostPointGhostPointPortion[neighborPointIndex][k] * coeff_y[j];
				}
			}

			// interpolate on y direction
			neighborPoint_x = ninePointScheme[i + 3 * j][0];
			neighborPoint_y = ninePointScheme[i + 3 * j][1];
			parient_x = ninePointScheme[i + 3 * j][2];
			parient_y = ninePointScheme[i + 3 * j][3];
			neighborPointIndex = ninePointScheme[i + 3 * j][4];
			// case the point is not ghost point
			if (neighborPointIndex == -1) {
				rowdispGridPointPortion_y[i][Ny * neighborPoint_x + neighborPoint_y] += coeff_y[j];
				rowdispGridPointPortion_y[i][Ny * neighborPoint_x + neighborPoint_y + 1] += coeff_y[j];
			}
			// case the point is ghost point
			else {
				// add grid point portion
				for (int k = 0; k < ghostPointGridPointPortion[neighborPointIndex].size(); ++k) {
					rowdispGridPointPortion_y[i][k] += ghostPointGridPointPortion[neighborPointIndex][k] * coeff_y[j];
				}
				// add ghost point portion
				for (int k = 0; k < ghostPointGhostPointPortion[neighborPointIndex].size(); ++k) {
					rowdispGridPointPortion_y[i][k] += ghostPointGhostPointPortion[neighborPointIndex][k] * coeff_y[j];
				}
			}
		}
	}
	// total disp info
	for (int k = 0; k < dispGridPointPortion.size(); ++k) {
		dispGridPointPortion[k] += coeff_y[0] * rowdispGridPointPortion_x[0][k] 
			+ coeff_y[1] * rowdispGridPointPortion_x[1][k]
			+ coeff_y[2] * rowdispGridPointPortion_x[2][k];
	}
	// Lagrange polynomial coefficient for first order derivative
	double coeff_deri_x[3];
	double coeff_deri_y[3];
	coeff_deri_x[0] = (2 * posi_x - (x + x + 1.) * dx) / ((-dx) * (-2 * dx));
	coeff_deri_x[1] = (2 * posi_x - (x - 1. + x + 1.) * dx) / ((dx) * (-dx));
	coeff_deri_x[2] = (2 * posi_x - (x - 1. + x) * dx) / ((dx) * (2 * dx));
	coeff_deri_y[0] = (2 * posi_y - (y + y + 1.) * dy) / ((-dy) * (-2 * dy));
	coeff_deri_y[1] = (2 * posi_y - (y - 1. + y + 1.) * dy) / ((dy) * (-dy));
	coeff_deri_y[2] = (2 * posi_y - (y - 1. + y) * dy) / ((dy) * (2 * dy));
	// x derivative info
	for (int k = 0; k < rowdispGridPointPortion_y[0].size(); ++k) {
		deriGridPointPortion_x[k] += coeff_deri_x[0] * rowdispGridPointPortion_y[0][k]
			+ coeff_deri_x[1] * rowdispGridPointPortion_y[1][k]
			+ coeff_deri_x[2] * rowdispGridPointPortion_y[2][k];
	}
	for (int k = 0; k < rowdispGhostPointPortion_y[0].size(); ++k) {
		deriGhostPointPortion_x[k] += coeff_deri_x[0] * rowdispGhostPointPortion_y[0][k]
			+ coeff_deri_x[1] * rowdispGhostPointPortion_y[1][k]
			+ coeff_deri_x[2] * rowdispGhostPointPortion_y[2][k];
	}
	// y derivative info
	for (int k = 0; k < rowdispGridPointPortion_x[0].size(); ++k) {
		deriGridPointPortion_y[k] += coeff_deri_y[0] * rowdispGridPointPortion_x[0][k]
			+ coeff_deri_y[1] * rowdispGridPointPortion_x[1][k]
			+ coeff_deri_y[2] * rowdispGridPointPortion_x[2][k];
	}
	for (int k = 0; k < rowdispGhostPointPortion_x[0].size(); ++k) {
		deriGhostPointPortion_y[k] += coeff_deri_y[0] * rowdispGhostPointPortion_x[0][k]
			+ coeff_deri_y[1] * rowdispGhostPointPortion_x[1][k]
			+ coeff_deri_y[2] * rowdispGhostPointPortion_x[2][k];
	}
}

void Geometry::ThreePointInterpolation(std::vector<double>& dispGridPointPortion, std::vector<double>& dispGhostPointPortion,
	std::vector<double>& deriGridPointPortion, std::vector<double>& deriGhostPointPortion, bool dir) {

}

void Geometry::GetNeighborPointInfo(int offset_x, int offset_y, int central_x, int central_y, std::vector<int>& pointInfo) {
	//Get coordinates of current neighbor point
	//pointInfo[0]: x coordinate of current neighbor point
	//pointInfo[1]: y coordinate of current neighbor point
	//pointInfo[2]: x coordinate of parent point
	//pointInfo[3]: y coordinate of parent point
	//pointInfo[4]: Indicator of point type (-1: not ghost point; n: ghost point index)
	int x = central_x + offset_x;
	int y = central_y + offset_y;
	pointInfo[0] = x;
	pointInfo[1] = y;
	//Case that current neighbor point is not ghost point
	if (pointMatrix[central_x][central_y].point_material == pointMatrix[x][y].point_material) {
		pointInfo[2] = x;
		pointInfo[3] = y;
		pointInfo[4] = -1;
	}
	//Case that current neighbor point is ghost point
	//According to nine point scheme assumption, as long as the neighbor point belongs to
	//another material, then it must be a ghost point (normally it is true, but if the mesh
	//is too coarse, then this assumption might be violated especially for boundary point
	//example:
	//  o o o
	//  x o o
	//  x x x
	//In above example point at position (1,2) has one neighbor point at position (3,3)
	//which is not ghost point but it belongs to another material

	//TODO: In future the assumption should be prechecked when generating mesh
	//TODO: Adaptive finite difference method can be used to solve this problem
	else {
		switch (pointMatrix[x][y].point_type[0]) {
		case -1://Neighbor point is type -1 (its corresponding parent point is type 2)
			switch (pointMatrix[x][y].point_type[1]) {
			case 1://Neighbor point is subtype 1
				pointInfo[2] = x + 1;
				pointInfo[3] = y + 1;
				pointInfo[4] = ghostPointIndex[x + 1][y + 1] + 2;
				break;
			case 2://Neighbor point is subtype 2
				pointInfo[2] = x - 1;
				pointInfo[3] = y + 1;
				pointInfo[4] = ghostPointIndex[x - 1][y + 1] + 2;
				break;
			case 3://Neighbor point is subtype 3
				pointInfo[2] = x - 1;
				pointInfo[3] = y - 1;
				pointInfo[4] = ghostPointIndex[x - 1][y - 1] + 2;
				break;
			case 4://Neighbor point is subtype 4
				pointInfo[2] = x + 1;
				pointInfo[3] = y - 1;
				pointInfo[4] = ghostPointIndex[x + 1][y - 1] + 2;
				break;
			}
			break;
		case 0:
			//TODO: this case should not be runned, implement assert function in the future
			break;
		case 1://Neighbor point is type 1 (its corresponding parent point is type 1 or 2)
			switch (pointMatrix[x][y].point_type[1]) {
			case 1://Neighbor point is subtype 1
				pointInfo[2] = x + 1;
				pointInfo[3] = y;
				pointInfo[4] = (pointMatrix[x + 1][y].point_type[0] == 1) ? ghostPointIndex[x + 1][y] : ghostPointIndex[x + 1][y];
				break;
			case 2://Neighbor point is subtype 2
				pointInfo[2] = x;
				pointInfo[3] = y + 1;
				pointInfo[4] = (pointMatrix[x][y + 1].point_type[0] == 1) ? ghostPointIndex[x][y + 1] : ghostPointIndex[x][y + 1] + 1;
				break;
			case 3://Neighbor point is subtype 3
				pointInfo[2] = x - 1;
				pointInfo[3] = y;
				pointInfo[4] = (pointMatrix[x - 1][y].point_type[0] == 1) ? ghostPointIndex[x - 1][y] : ghostPointIndex[x - 1][y];
				break;
			case 4://Neighbor point is subtype 4
				pointInfo[2] = x;
				pointInfo[3] = y - 1;
				pointInfo[4] = (pointMatrix[x][y - 1].point_type[0] == 1) ? ghostPointIndex[x][y - 1] : ghostPointIndex[x][y - 1] + 1;
				break;
			}
			break;
		case 2://Neighbor point is type 2
			//In this case we only consider that offset_x = offset_y = 0 such that central point is the parent
			//point of current neighbor point (Boundary central point cases can be ignored as long as the mesh
			//is dense enough)
			//TODO: Add boundary central point processing code in the future
			if (offset_x == 0) {
				pointInfo[2] = x;
				pointInfo[3] = offset_y > 0 ? y - 1 : y + 1;
				pointInfo[4] = pointMatrix[x][pointInfo[3]].point_type[0] == 1 ? ghostPointIndex[x][pointInfo[3]] : ghostPointIndex[x][pointInfo[3]] + 1;
			}
			else if (offset_y == 0) {
				pointInfo[2] = offset_x > 0 ? x - 1 : x + 1;
				pointInfo[3] = y;
				pointInfo[4] = pointMatrix[pointInfo[2]][y].point_type[0] == 1 ? ghostPointIndex[pointInfo[2]][y] : ghostPointIndex[pointInfo[2]][y];
			}
			else {
				pointInfo[2] = offset_x > 0 ? x - 1 : x + 1;
				pointInfo[3] = offset_y > 0 ? y - 1 : y + 1;
				pointInfo[4] = ghostPointIndex[pointInfo[2]][pointInfo[3]] + 2;
			}
			break;
		default:
			//TODO: This case should not be runned, implement assert function in the future
			break;
		}
	}
}

// Compute fractional derivative according to point position
void Geometry::ComputeFracDeri(double posi_x, double posi_y, short materialType) {
	double curNlcLength = materials[materialType].nlc_range;
	double curAlpha = materials[materialType].alpha;
	double NlcLength[4];
}

// Compute fractional derivative according to point coordinates
void Geometry::ComputeFirstOrdPrecisionFracDeri(int x, int y, std::vector<double>& gridPointPortion, std::vector<double>& ghostPointPortion, bool dir) {
	double curNlcLength = materials[pointMatrix[x][y].point_material].nlc_range;
	double curAlpha = materials[pointMatrix[x][y].point_material].alpha;
	double posi_x = x * dx;
	double posi_y = y * dy;
	if (std::abs((curAlpha - (double)1)) < 1.e-9) {

	}
	else if (curAlpha > 0 && curAlpha < 1) {

	}
	// alpha cannot be either smaller than 0 or larger than 1
	else {
		exit(1);
	}
} 

void Geometry::ComputeSecOrdPrecisionFracDeri(int x, int y, std::vector<double>& gridPointPortion, std::vector<double>& ghostPointPortion, bool dir) {
	double lengthScale = materials[pointMatrix[x][y].point_material].nlc_range;
	double curAlpha = materials[pointMatrix[x][y].point_material].alpha;
	double posi_x = x * dx;
	double posi_y = y * dy;
	// vector stores nonlocal range on 4 directions:
	// curNlcRange[0]: right; curNlcRange[1]: top; curNlcRange[2]: left; curNlcRange[3]: bottom
	std::vector<double> curNlcRange = std::vector<double>(4, 0);
	// vector stores the farthest grid point coordinates on 4 directions
	// if the farthest point is not grid point, then it stores the second farthest grid point coordinates
	std::vector<std::vector<int>> curNlcGridPoints;
	// alpha==1, then fractional derivative is first order derivative
	std::vector<double> curNlcNonGridPoints = std::vector<double>(4, 0);
	if (std::abs((curAlpha - (double)1)) < 1.e-9) {
		switch (dir) {
		// x derivative
		case false:
			gridPointPortion = firstOrdDeriCentralGridPointPortion_x[(Ny + 1) * x + y];
			ghostPointPortion = firstOrdDeriCentralGhostPointPortion_x[(Ny + 1) * x + y];
			break;
		// y derivative
		case true:
			gridPointPortion = firstOrdDeriCentralGridPointPortion_y[(Ny + 1) * x + y];
			ghostPointPortion = firstOrdDeriCentralGhostPointPortion_y[(Ny + 1) * x + y];
			break;
		default:
			exit(1);
			break;
		}
	}
	else if (curAlpha > 0 && curAlpha < 1) {
		GetNlcInfo(curNlcRange, curNlcGridPoints, curNlcNonGridPoints, x, y, lengthScale);
	}
	// alpha cannot be either smaller than 0 or larger than 1
	else {
		exit(1);
	}
}

void Geometry::GetNlcInfo(std::vector<double>& nlcRange, std::vector<std::vector<int>>& nlcGridPoints, std::vector<double>& nlcNonGridPoints,
	int x, int y, double lengthScale, double pposi_x = -1, double pposi_y=-1) {
	double posi_x, posi_y;
	if (posi_x < 0) {
		posi_x = x * dx;
		posi_y = x * dy;
	}
	else {
		posi_x = pposi_x;
		posi_y = pposi_y;
		int originalNlcBoun[4] = { posi_y,posi_x,posi_y,posi_x };
	}
	double originalNlcBoun[4] = { posi_x + lengthScale,posi_y + lengthScale,posi_x - lengthScale,posi_y - lengthScale };
	double ExtremeNlcBoun[4] = { std::sqrt(std::pow(R,2) - std::pow(y,2)), std::sqrt(std::pow(R,2) - std::pow(x,2)),
		std::sqrt(std::pow(R,2) - std::pow(y,2)), std::sqrt(std::pow(R,2) - std::pow(x,2)) };
	nlcRange[0] = ExtremeNlcBoun[0] < originalNlcBoun[0] ? ExtremeNlcBoun[0] : originalNlcBoun[0];
	nlcRange[1] = ExtremeNlcBoun[1] < originalNlcBoun[1] ? ExtremeNlcBoun[1] : originalNlcBoun[1];
	nlcRange[2] = ExtremeNlcBoun[2] > originalNlcBoun[2] ? ExtremeNlcBoun[2] : originalNlcBoun[2];
	nlcRange[3] = ExtremeNlcBoun[3] > originalNlcBoun[3] ? ExtremeNlcBoun[3] : originalNlcBoun[3];

	int floor_x, floor_y, ceil_x, ceil_y;
	for (int i = 0; i < 2; ++i) {
		// get nlc information on x direction
		floor_x = std::floor(originalNlcBoun[2 * i] / dx);
		ceil_x = std::floor(originalNlcBoun[2 * i] / dx);
	}
}