namespace integration
{
	class Gauss_integration
	{
	public:
		double gauss_points[2][9];//����� ������
		double gauss_weights[9];// ���� ������
		double gauss_points_1[3];//����� ������
		double gauss_weights_1[3];// ���� ������

		void initialize();
	};

}