#include "ScriptJMMV.h"
#include <chrono>

int main() {
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	vvi in_1 = Smolyak_Elem_Isotrop();

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count();
	std::cout << "Time Duration of Smolyak_Elem_Isotrop: " << duration << std::endl;
	std::cout << "Mat length: " << static_cast<int>(in_1.size()) << std::endl;

	t1 = std::chrono::high_resolution_clock::now();
	vvd in_2 = Smolyak_Grid(in_1);
	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count();
	std::cout << "Time Duration of Smolyak_Grid: " << duration << std::endl;

	/*
	t1 = std::chrono::high_resolution_clock::now();
	vd in_3 = Smolyak_Polynomial(in_2, in_1);
	//double *in_3 = Smolyak_Polynomial(in_2, in_1);

	t2 = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds> (t2 - t1).count();
	std::cout << "Duration of Smolyak_Polynomial: " << duration << std::endl;
	*/
	//delete[] in_3;

	return 0;
}
