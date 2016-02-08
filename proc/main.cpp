#include <gatb/gatb_core.hpp>
#include <unordered_map>
#include <boost/unordered_set.hpp>
#include <unordered_set>
#include "kseq.h"

#include <gatb/tools/math/Integer.hpp>

const int32_t k=20;
const int32_t fasta_line_length=60;
const int32_t max_contig_length=10000000;

KSEQ_INIT(gzFile, gzread);


typedef uint64_t nkmer_t;
typedef std::set<nkmer_t> ordered_set_t;
typedef std::unordered_set<nkmer_t> unordered_set_t;

char complement(char x){
	switch(x){
		case 'A':
		case 'a':
			return 'T';
		case 'C':
		case 'c':
			return 'G';
		case 'G':
		case 'g':
			return 'C';
		case 'T':
		case 't':
			return 'A';
		default:
			return 'N';
	}
}

struct contig_t{
	int32_t k;

	char *seq_buffer;
	char *r_ext;
	char *l_ext;

	contig_t(const char *base_kmer){
		k=strlen(base_kmer);
		seq_buffer=new char[k+2*max_contig_length+1]();

		reinit(base_kmer);
	}

	int32_t reinit(const char *base_kmer){
		assert(strlen(base_kmer)==k);

		l_ext = r_ext = &seq_buffer[max_contig_length];
		*r_ext='\0';

		for(int32_t i=0;i<k;i++){
			r_extend(base_kmer[i]);
		}
	}

	int32_t r_extend(char c){
		*r_ext=c;
		++r_ext;
		*r_ext='\0';
	}

	int32_t l_extend(char c){
		--l_ext;
		*l_ext=complement(c);
	}

	~contig_t(){
		delete[] seq_buffer;
	}

	int32_t print_to_fasta(FILE* fasta_file, char* contig_name, char *comment){
		if (comment==nullptr){
			fprintf(fasta_file,">%s\n",contig_name);
		} else {
			fprintf(fasta_file,">%s %s\n",contig_name,comment);
		}

		int32_t buffer_len=r_ext-l_ext;
		char print_buffer[fasta_line_length+1]={0};

		int32_t j=0;

		for (char *p=l_ext;p<r_ext;p+=fasta_line_length){
			strncpy(print_buffer,p,fasta_line_length);
			fprintf(fasta_file, "%s\n", print_buffer);
		}		

		delete[] seq_buffer;
	}
};


/*
	TODO: test if kmer is correct
*/
int ordered_kmers_from_fasta(const string &fasta_fn, ordered_set_t &ordered_set, int32_t k){
	ordered_set.clear();

 	gzFile fp;
    kseq_t *seq;
    int64_t l;

    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);

    char buffer[100]={};

    Kmer<>::ModelCanonical model (k);


	for(int32_t seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
 		//seq->seq.s;
	    Data data (seq->seq.s);
	
	    Kmer<>::ModelCanonical::Iterator it (model);
	    it.setData (data);

	    for (it.first(); !it.isDone(); it.next())
	    {
	    	std::cout << it->value() << " " << model.toString( it->value() ) << std::endl;
	    }

    }

    kseq_destroy(seq);
    gzclose(fp);

	return 0;
}

int unsort_kmers(const ordered_set_t &ordered_set, unordered_set_t &unordered_set){
	unordered_set.clear();

    for(nkmer_t const &element: ordered_set){
     	unordered_set.insert(element);
    }
	return 0;
}

int ordered_kmers_to_fasta(const string &fasta_fn, ordered_set_t &ordered_set, int32_t k){
	return 0;
}



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
	ordered_set_t a;
	ordered_kmers_from_fasta("test.fa", a, 5);
	return 0;


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