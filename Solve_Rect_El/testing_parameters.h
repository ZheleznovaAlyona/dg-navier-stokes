#include <string>
namespace testingparameters
{

	class Testing_parameters
	{
	public:

		static bool use_LU;
		static int test;
		static int solver;
		static void Testing_parameters::initialize(std::string file_name);
	private:
		Testing_parameters();
		~Testing_parameters();
		Testing_parameters(Testing_parameters& parameters);
	};

}