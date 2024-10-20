#include "tdd_calculation.h"

void TDD(Profile p, Results r)
{
	Eigen::ArrayXd dD(p.h.rows());

	dD = 4 * r.misfit_density; // TODo serve capire come gestire diff eq.

	r.D = dD;
}
