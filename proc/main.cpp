#include "kseq.h"

#include <cstdio>
#include <climits>

#include <vector>
#include <unordered_set>
#include <boost/unordered_set.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <gatb/gatb_core.hpp>
#include <gatb/tools/math/Integer.hpp>

const int32_t k=22;
const int32_t fasta_line_length=60;
const int32_t max_contig_length=1000000;

static const uint8_t nt4_nt256[] = "ACGTN";


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

	char *r_ext_border;
	char *l_ext_border;

	contig_t(uint32_t k){
		this->k=k;
		seq_buffer=new char[k+2*max_contig_length+1]();
		l_ext_border=seq_buffer;
		r_ext_border=seq_buffer+k+2*max_contig_length;
	}

	int32_t reinit(const char *base_kmer){
		assert(strlen(base_kmer)==k);

		l_ext = r_ext = &seq_buffer[max_contig_length];
		*r_ext='\0';

		for(int32_t i=0;i<k;i++){
			r_extend(base_kmer[i]);
		}
		return 0;
	}

	int32_t r_extend(char c){
		*r_ext=c;
		++r_ext;
		*r_ext='\0';
		return 0;
	}

	int32_t l_extend(char c){
		--l_ext;
		*l_ext=complement(c);
		return 0;
	}

	~contig_t(){
		delete[] seq_buffer;
	}

	bool is_full(){
		return (r_ext>=r_ext_border) || (l_ext<=l_ext_border);
	}

	int32_t print_to_fasta(FILE* fasta_file, char* contig_name, char *comment=nullptr) const {
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

		return 0;
	}
};


/*
	TODO: test if kmer is correct
*/
template <typename T>
int kmers_from_fasta(const string &fasta_fn, T &set, int32_t k){
	set.clear();

 	gzFile fp;
    kseq_t *seq;
    int64_t l;

    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);

    char buffer[100]={};

    Kmer<>::ModelCanonical model (k);


	for(int32_t seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
		cout << "start" << endl;
	    Data data (seq->seq.s);	
		cout << "data ok" << endl;
	    Kmer<>::ModelCanonical::Iterator it (model);
	    it.setData (data);

		cout << "starting iterator" << endl;
	    for (it.first(); !it.isDone(); it.next())
	    {
	    	set.insert(it->value().getVal());
	    	//std::cout << it->value() << " " << model.toString( it->value() ) << std::endl;
	    }
		cout << "iterator finished" << endl;
    }

    kseq_destroy(seq);
    gzclose(fp);

	return 0;
}

template <typename T1, typename T2>
int32_t find_intersection(const std::vector<T1> &sets, T2 &intersection){
	assert(sets.len()>0);

	/* 
		1) Find the smallest set from sets.
	*/

	int32_t min=INT32_MAX;
	int32_t i_min=-1;

	cout << "searching min" << endl;

	for(int32_t i=0;i<sets.size();i++){
		if (sets[i].size()<min){
			min=sets[i].size();
			i_min=i;
			cout << "new min" << i << endl;
		}
	}

	assert(i_min != Int32.MaxValue && i_min!=-1);

	/*
		2) Take it as the intersection.
	*/

	cout << "2" << endl;

	intersection.clear();
	cout << "2.1" << endl;
	std::copy(
		sets[i_min].cbegin(), sets[i_min].cend(),
		std::inserter(intersection,intersection.end())
	);

	cout << "first intersectino size " << intersection.size() << endl;

	/*
		3) Remove elements from intersection present in other sets.
	*/

	for(const T1 &current_set : sets) {

		for(boost::unordered_set<nkmer_t>::iterator it = intersection.begin(); it !=intersection.end();){

			if(current_set.find(*it) == current_set.cend()){
				it=intersection.erase(it);
			}
			else{
				++it;
			}
		}
	}

	return 0;
}


template <typename T1, typename T2>
int32_t remove_subset(std::vector<T1> &sets, const T2 &subset){

	for(int32_t i=0;i<sets.size();i++){

		T1 &current_set = sets[i];

		for(const nkmer_t &nkmer : subset){
			current_set.erase(nkmer);
		}
	}

}


template <typename T>
int assemble(const string &fasta_fn, T &set, int32_t k){
    Kmer<>::ModelCanonical model (k);


    cout << "assembling, size: " << set.size() << endl;

	FILE *file=fopen(fasta_fn.c_str(),"w+");

	char kmer_str[k+1];

	contig_t contig(k);

	const std::vector<char> nucls = {'A','C','G','T'};

	int i=0;
	while(set.size()>0){
		i++;
		//printf("writing contig %d\n",i);
		const nkmer_t central_kmer=*(set.begin());
		set.erase(central_kmer);

		std::string central_kmer_string=model.toString(central_kmer);
		contig.reinit(central_kmer_string.c_str());

		strncpy(kmer_str,central_kmer_string.c_str(),k);
		kmer_str[k]='\0';

		//printf("central k-mer: %s\n",central_kmer_string.c_str());

		bool extending = true;
		while (extending){
			
			for(int32_t i=0;i<k;i++){
				kmer_str[i]=kmer_str[i+1];
			}
			kmer_str[k]='\0';


			extending=false;
			for(const char &c : nucls){
				kmer_str[k-1]=c;


				nkmer_t nkmer = model.codeSeed (
 					kmer_str , Data::ASCII
				).value().getVal();


				if(set.count( nkmer )){
					contig.r_extend(c);
					extending=true;
					set.erase(nkmer);
					break;
				}

			}

			if(contig.is_full()){
				extending=false;
			}
		}

		contig.print_to_fasta(file,"conting xx");
	}

	fclose(file);

    cout << "...finished " << endl;

	return 0;

}


int main (int argc, char* argv[])
{

	std::string intersection_fn;
	std::vector<std::string> input_fns;
	std::vector<std::string> output_fns;

    try
    {
        namespace po = boost::program_options;
        
        po::positional_options_description pos;
        pos.add("input-file", -1);
        
        po::options_description vol("Command-line parameters");

        vol.add_options()
               ("input,i", po::value<std::vector<std::string>>(&input_fns)->required(), "Input files.")
               ("output,o", po::value<std::vector<std::string>>(&output_fns)->required(), "Output files.")
               ("intersection,x", po::value<std::string>(&intersection_fn)->required(), "Intersection file.")
               ;


		po::variables_map vm;
        try
        {
            po::store(po::command_line_parser(argc, argv).options(vol).positional(pos).run(),vm); // can throw
            po::notify(vm); // throws on error, so do after help in case there are any problems
            
            if (input_fns.size() != output_fns.size()) {
                fprintf(stderr,"There must be equal number of input and output files.\n");
		        return EXIT_FAILURE;
            }
            
        }
        catch(po::error& e)
        {
            std::cout << vol << "\n";
            fprintf(stderr,"%s.\n",e.what());
	        return EXIT_FAILURE;
        }
        
    }
    catch(std::exception& e)
    {
        fprintf(stderr,"Unhandled Exception: %s.\n",e.what());
        return EXIT_FAILURE;
    }

	std::vector< boost::unordered_set<nkmer_t> > full_sets(output_fns.size());

    for(int32_t i=0;i<input_fns.size();i++){
    	cout << "Loading " << input_fns[i] << endl;
    	kmers_from_fasta(input_fns[i],full_sets[i],k);
    }

    boost::unordered_set<nkmer_t> intersection;

	cout << "Intersection" << endl;

    find_intersection(full_sets, intersection);
    cout << "intersection size " <<  intersection.size() << endl;

	cout << "Removing subsets" << endl;

    remove_subset(full_sets, intersection);

    for(int32_t i=0;i<input_fns.size();i++){
    	cout << "Assembling " << input_fns[i] << endl;
		assemble(output_fns[i],full_sets[i],k);
    }

   	cout << "Assembling intersection" << endl;
	assemble(intersection_fn,intersection,k);


    return 0;
}