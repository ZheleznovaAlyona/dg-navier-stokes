#include "SLAE.h"

namespace slae
{

	void SLAE::initialize(int size, int unzero_matrix_elements)
	{
		n = size;
		A.initialize(n, unzero_matrix_elements);
		b.initialize(n);
	}

	void SLAE::reinitialize()
	{
		A.reinitialize();
		b.make_zero();
	}

}
