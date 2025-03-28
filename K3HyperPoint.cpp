#include "K3HyperPoint.h"

K3HyperPoint::K3HyperPoint(int n) {
	setN(n);
}

void K3HyperPoint::setN(int n) {
	coords.resize(n);
}


