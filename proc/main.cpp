#include <gatb/gatb_core.hpp>

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how instantiate the Model class and how to get info       */
/* about it.                                                                    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
	Kmer<>::ModelCanonical model;
	std::cout << "----------   MODEL   ----------"               << std::endl;
	std::cout << "span:             " << model.getSpan()         << std::endl;
	std::cout << "kmer size:        " << model.getKmerSize()     << std::endl;
	std::cout << "kmer memory size: " << model.getMemorySize()   << std::endl;
}


