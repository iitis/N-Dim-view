#include "K3HyperCloud.h"
#include "iostream"

K3HyperCloud::K3HyperCloud(int new_n, int new_m) : n(new_n), m(new_m) { // m points in n-dimensional space
	K3RawData.resize(m);
	for (int i = 0; i < m; i++) {
		K3RawData[i].setN(n);
	}
}

// K3HyperCloud::K3HyperCloud(int n, int m) {
//	typedef std::vector<double> K3HyperPoint;
//	typedef std::vector<K3HyperPoint> HyperPoints;
//	K3RawData=new public Eigen::HyperPoints()
// }
//
