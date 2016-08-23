#include <string>
namespace testingparameters
{

	class Testing_parameters
	{
	public:

		bool use_LU;
		int test;
		int solver;
		void Testing_parameters::initialize(std::string file_name);
	};

}