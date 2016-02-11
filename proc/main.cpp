#include "kseq.h"

#include <zlib.h>

#include <cstdio>
#include <limits>
#include <iostream>

#include <vector>
//#include <unordered_set>
#include <boost/unordered_set.hpp>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

const int32_t k=22;
const int32_t fasta_line_length=60;
const int32_t max_contig_length=1000000;

static const uint8_t nt4_nt256[] = "ACGTN";

static const uint8_t nt256_nt4[] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };


KSEQ_INIT(gzFile, gzread);


typedef __uint128_t nkmer_t;

//typedef std::set<nkmer_t> ordered_set_t;
//typedef std::unordered_set<nkmer_t> unordered_set_t;

typedef std::set<nkmer_t> set_t;


template<typename _nkmer_T>
int32_t encode_forward(const char *kmers, const int32_t k, _nkmer_T &nkmer){
	nkmer=0;
	for (int32_t i=0;i<k;i++){
		uint8_t nt4 = nt256_nt4[kmers[i]];
		if (nt4==4){
			return -1;
		}

		nkmer <<= 2;
		nkmer |= nt4;
	}
	return 0;
}

template<typename _nkmer_T>
int32_t encode_reverse(const char *kmers, const int32_t k, _nkmer_T &nkmer){
	nkmer=0;
	for (int32_t i=0;i<k;i++){
		uint8_t nt4 = nt256_nt4[kmers[k-i-1]];
		if (nt4==4){
			return -1;
		}

		//complement
		nt4 = 3-nt4;

		nkmer <<= 2;
		nkmer |= nt4;
	}
	return 0;
}

template<typename _nkmer_T>
int32_t encode_canonical(const char *kmers, const int32_t k, _nkmer_T &nkmer){
	_nkmer_T nkmer_f;
	_nkmer_T nkmer_r;

	int32_t error_code;

	error_code=encode_forward(kmers, k, nkmer_f);
	if(error_code!=0){
		return error_code;
	}

	error_code=encode_reverse(kmers, k, nkmer_r);
	if(error_code!=0){
		return error_code;
	}

	nkmer=std::min(nkmer_f,nkmer_r);

	return 0;
}

template<typename _nkmer_T>
int32_t decode_kmer(_nkmer_T nkmer, int32_t k, std::string &kmer){
	kmer.resize(k);
	for(int32_t i=0;i<k;i++){
		uint8_t nt4 = nkmer & 0x3;
		nkmer >>=2;
		kmer[k-i-1]=nt4_nt256[nt4];
	}

	return 0;
}

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
		uint8_t nt4 = nt256_nt4[c];

		if (nt4==4){
			return -1;
		}

		*r_ext=nt4_nt256[nt4];
		++r_ext;
		*r_ext='\0';
		return 0;
	}

	int32_t l_extend(char c){
		uint8_t nt4 = nt256_nt4[c];

		if (nt4==4){
			return -1;
		}

		--l_ext;
		*l_ext=nt4_nt256[3-nt4];
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

//template<typename _nkmer_T, typename _set_T>
template<typename _set_T>
int kmers_from_fasta(const std::string &fasta_fn, _set_T &set, int32_t k){

	set.clear();

 	gzFile fp;
    kseq_t *seq;
    int64_t l;

    fp = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);


    char buffer[100]={};

	typename _set_T::value_type nkmer;

	for(int32_t seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
		std::cout << "kmers from fasta" << std::endl;

		std::cout << "starting iterator" << std::endl;
		for(char *kmer=seq->seq.s; kmer < (seq->seq.s) + (seq->seq.l) -k +1 ;kmer++){
			encode_canonical(kmer, k, nkmer);
			set.insert(nkmer);
		}
		std::cout << "iterator finished" << std::endl;

    }

    kseq_destroy(seq);
    gzclose(fp);

	return 0;
}

template<typename _set_T>
int32_t find_intersection(const std::vector<_set_T> &sets, _set_T &intersection){
	assert(sets.size()>0);

	/* 
		1) Find the smallest set from sets.
	*/

	int32_t min=std::numeric_limits<int32_t>::max();
	int32_t i_min=-1;

	std::cout << "searching min" << std::endl;

	for(int32_t i=0;i<sets.size();i++){
		if (sets[i].size()<min){
			min=sets[i].size();
			i_min=i;
			std::cout << "new min" << i << std::endl;
		}
	}

	assert(i_min != std::numeric_limits<int32_t>::max() && i_min!=-1);

	/*
		2) Take it as the intersection.
	*/

	std::cout << "2" << std::endl;

	intersection.clear();
	std::cout << "2.1" << std::endl;
	std::copy(
		sets[i_min].cbegin(), sets[i_min].cend(),
		std::inserter(intersection,intersection.end())
	);

	std::cout << "first intersectino size " << intersection.size() << std::endl;

	/*
		3) Remove elements from intersection present in other sets.
	*/

	for(const _set_T &current_set : sets) {

		for(auto it = intersection.begin(); it !=intersection.end();){

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


template<typename _set_T, typename _subset_T>
int32_t remove_subset(std::vector<_set_T> &sets, const _subset_T &subset){

	for(int32_t i=0;i<sets.size();i++){

		_set_T &current_set = sets[i];

		for(const auto &nkmer : subset){
			current_set.erase(nkmer);
		}
	}

}


template<typename _set_T>
int assemble(const std::string &fasta_fn, _set_T &set, int32_t k){

    std::cout << "assembling, size: " << set.size() << std::endl;

	FILE *file=fopen(fasta_fn.c_str(),"w+");

	char kmer_str[k+1];

	contig_t contig(k);

	const std::vector<char> nucls = {'A','C','G','T'};

	int i=0;
	while(set.size()>0){
		i++;
		//printf("writing contig %d\n",i);

		const auto central_nkmer=*(set.begin());
		set.erase(central_nkmer);

		std::string central_kmer_string;
		decode_kmer(central_nkmer,k,central_kmer_string);
		contig.reinit(central_kmer_string.c_str());

		strncpy(kmer_str,central_kmer_string.c_str(),k);
		kmer_str[k]='\0';

		//printf("central k-mer: %s\n",central_kmer_string.c_str());

		typename _set_T::value_type nkmer;

		bool extending = true;
		while (extending){
			
			for(int32_t i=0;i<k;i++){
				kmer_str[i]=kmer_str[i+1];
			}
			kmer_str[k]='\0';


			extending=false;
			for(const char &c : nucls){
				kmer_str[k-1]=c;

				encode_canonical(kmer_str, k, nkmer);

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

    std::cout << "...finished " << std::endl;

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
    	std::cout << "Loading " << input_fns[i] << std::endl;
    	kmers_from_fasta(input_fns[i],full_sets[i],k);
    }

    boost::unordered_set<nkmer_t> intersection;

	std::cout << "Intersection" << std::endl;

    find_intersection(full_sets, intersection);
    std::cout << "intersection size " <<  intersection.size() << std::endl;

	std::cout << "Removing subsets" << std::endl;

    remove_subset(full_sets, intersection);

    for(int32_t i=0;i<input_fns.size();i++){
    	std::cout << "Assembling " << input_fns[i] << std::endl;
		assemble(output_fns[i],full_sets[i],k);
    }

   	std::cout << "Assembling intersection" << std::endl;
	assemble(intersection_fn,intersection,k);


    return 0;
}