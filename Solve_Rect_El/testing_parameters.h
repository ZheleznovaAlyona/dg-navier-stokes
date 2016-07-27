namespace testingparameters
{

	class Testing_parameters
	{
	public:

		static bool use_LU;
		static int test;
		static int solver;

		static const Testing_parameters& instance();
		
	private:
		Testing_parameters();
		Testing_parameters(bool use_LU_in, int test_in, int solver_in);
		~Testing_parameters();
		Testing_parameters(Testing_parameters& parameters);
		Testing_parameters& operator=(Testing_parameters&);
	};

}