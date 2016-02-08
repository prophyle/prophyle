#include <gatb/gatb_core.hpp>
#include <unordered_map>
#include <boost/unordered_set.hpp>
#include <unordered_set>

#include <gatb/tools/math/Integer.hpp>

//using boost::unordered_set;

/********************************************************************************/
/*                              Kmer management                                 */
/*                                                                              */
/* This snippet shows how instantiate the Model class and how to get info       */
/* about it.                                                                    */
/*                                                                              */
/********************************************************************************/
int main (int argc, char* argv[])
{
    // We define a sequence of nucleotides
    //const char* seq = "CATTGATAGTGGATGGT";
    const char* seq = "TTTTTTTCCACFGTCCCCCCCCCC";
    std::cout << "Initial sequence: " << seq << std::endl;

    // We configure a data object with a sequence (in ASCII format)
    Data data ((char*)seq);

    // We declare a kmer model with kmer of size 5.
    // Note that we want "direct" kmers, not the min(forward,revcomp) default behavior.
    Kmer<>::ModelCanonical model (5);

    // We declare an iterator on a given sequence.
    Kmer<>::ModelCanonical::Iterator it (model);

	Kmer<>::ModelCanonical::Kmer kmer;
    //std::unordered_set< uint64_t > s;
    boost::unordered_set< uint64_t > s;

    // We configure the iterator with our sequence
    it.setData (data);

    // We iterate the kmers.
    for (it.first(); !it.isDone(); it.next())
    {
		kmer = model.codeSeed (
			 model.toString(it->value()).c_str() , Data::ASCII
		);

        //std::cout << "kmer " << model.toString(it->value()) << ",  value " << it->value() << std::endl;
        s.insert(kmer.value().getVal());
    }

     for ( auto itt = s.begin(); itt != s.end(); ++itt ){
     	kmer.set(*itt);
     	cout << *itt << " " << model.toString( kmer.value() )  << endl;
    }

    for(auto const &element: s){
     	cout << element  << endl;
    }


}