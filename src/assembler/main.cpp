#include "kseq.h"

#include <zlib.h>

#include <cinttypes>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <set>
#include <cassert>
#include <sstream>
#include <unordered_set>
#include <getopt.h>

//typedef __uint128_t nkmer_t;
typedef uint64_t nkmer_t;
typedef std::set<nkmer_t> set_t;

const int32_t fasta_line_length=60;
const int32_t max_contig_length=10000000;
const int32_t max_allowed_kmer_length=sizeof(nkmer_t)*4;

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


KSEQ_INIT(gzFile, gzread)


void print_help(int k){
	std::cerr <<
		"\n" <<
		"Program: assembler (for k-mer propagation)\n" <<
		"Contact: Karel Brinda <karel.brinda@gmail.com>\n" <<
		"\n" <<
		"Usage:   assembler [options]\n" <<
		"\n"
		"Command-line parameters:\n" <<
		" -i arg  Input FASTA files.\n" <<
		" -o arg  Output FASTA files. They will contain the same \n" <<
		"            k-mers as input file except those from \n" <<
		"            intersection.\n" <<
		" -x arg  Intersection FASTA file.\n" <<
		" -k arg  K-mer size. [" << k << "]\n" << std::endl;
}


template<typename _nkmer_T>
int32_t encode_forward(const char *kmers, const int32_t k, _nkmer_T &nkmer){
	nkmer=0;
	for (int32_t i=0;i<k;i++){
		uint8_t nt4 = nt256_nt4[static_cast<int32_t>(kmers[i])];
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
		uint8_t nt4 = nt256_nt4[static_cast<int32_t>(kmers[k-i-1])];
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

void reverse_complement_in_place(std::string &kmer){
	//std::cerr << "before reverse complementing " << kmer << std::endl;
    std::reverse(kmer.begin(), kmer.end());
	for (int32_t i=0;i<static_cast<int32_t>(kmer.size());i++){
		char nt4=nt256_nt4[static_cast<int32_t>(kmer[i])];
		if (nt4<4){
			nt4=3-nt4;
		}
		kmer[i]=nt4_nt256[static_cast<int32_t>(nt4)];
	}
	//std::cerr << "after reverse complementing " << kmer << std::endl;
}

template<typename _set_T>
void debug_print_kmer_set(_set_T &set, int k){
	std::string kmer;
	for(auto x: set){
		decode_kmer(x, k, kmer);
		std::cerr << x << " " << kmer << ";  ";
	}
	std::cerr<<std::endl;
}


struct contig_t{
	int32_t k;

	/* contig buffer */
	char *seq_buffer;

	/* the first position of the contig */
	char *l_ext;

	/* the last position of the contig +1 (semiopen ) */
	char *r_ext;

	/* min possible value of l_ext */
	char *l_ext_border;

	/* max possible value of l_ext */
	char *r_ext_border;

	contig_t(uint32_t _k){
		this->k=_k;
		seq_buffer=new char[k+2*max_contig_length+1]();
		l_ext_border=seq_buffer;
		r_ext_border=seq_buffer+2*max_contig_length;
	}

	int32_t new_contig(const char *base_kmer){
		assert(static_cast<int32_t>(strlen(base_kmer))==k);

		l_ext = r_ext = &seq_buffer[max_contig_length];
		*r_ext='\0';

		for(int32_t i=0;i<k;i++){
			r_extend(base_kmer[i]);
		}
		return 0;
	}

	int32_t r_extend(char c){
		uint8_t nt4 = nt256_nt4[static_cast<int32_t>(c)];

		if (nt4==4){
			return -1;
		}

		*r_ext=nt4_nt256[nt4];
		++r_ext;
		*r_ext='\0';
		return 0;
	}

	int32_t l_extend(char c){
		uint8_t nt4 = nt256_nt4[static_cast<int32_t>(c)];

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

	int32_t print_to_fasta(FILE* fasta_file, const char* contig_name, const char *comment=nullptr) const {
		if (comment==nullptr){
			fprintf(fasta_file,">%s\n",contig_name);
		} else {
			fprintf(fasta_file,">%s %s\n",contig_name,comment);
		}

		char print_buffer[fasta_line_length+1]={0};

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

	std::cout << "Loading " << fasta_fn << std::endl;

	set.clear();

	gzFile fp;
	kseq_t *seq;
	int64_t l;

	fp = gzopen(fasta_fn.c_str(), "r");
	seq = kseq_init(fp);

	typename _set_T::value_type nkmer;

	for(int32_t seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
		//std::cout << "kmers from fasta" << std::endl;

		//std::cout << "starting iterator" << std::endl;
		for(char *kmer=seq->seq.s; kmer < (seq->seq.s) + (seq->seq.l) -k +1 ;kmer++){
			int c = encode_canonical(kmer, k, nkmer);
			if (c==0){
				set.insert(nkmer);
			}
			else{
				//std::cout << "problem" <<std::endl;
			}
		}
	}

	std::cout << "  " << fasta_fn << " loaded - #kmers: " << set.size()<<std::endl;
	//std::cout << "iterator finished" << std::endl;

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

	//std::cout << "searching min" << std::endl;

	for(int32_t i=0;i<static_cast<int32_t>(sets.size());i++){
		if (static_cast<int32_t>(sets[i].size())<min){
			min=sets[i].size();
			i_min=i;
			//std::cout << "new min" << i << std::endl;
		}
	}

	assert(i_min != std::numeric_limits<int32_t>::max() && i_min!=-1);

	/*
		2) Take it as the intersection.
	*/

	//std::cout << "2" << std::endl;

	intersection.clear();
	//std::cout << "2.1" << std::endl;
	std::copy(
		sets[i_min].cbegin(), sets[i_min].cend(),
		std::inserter(intersection,intersection.end())
	);

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

	for(int32_t i=0;i<static_cast<int32_t>(sets.size());i++){

		_set_T &current_set = sets[i];

		for(const auto &nkmer : subset){
			current_set.erase(nkmer);
		}
	}

	return 0;
}


template<typename _set_T>
int assemble(const std::string &fasta_fn, _set_T &set, int32_t k){

	std::cout << "Assembling " << fasta_fn << " from " << set.size() << " kmers" << std::endl;

	FILE *file=fopen(fasta_fn.c_str(),"w+");

	char kmer_str[max_allowed_kmer_length+1];

	contig_t contig(k);

	const std::vector<char> nucls = {'A','C','G','T'};

	//int32_t i=0;
	int32_t contig_id=0;
	while(set.size()>0){

		const auto central_nkmer=*(set.begin());
		set.erase(central_nkmer);

		std::string central_kmer_string;
		decode_kmer(central_nkmer,k,central_kmer_string);
		contig.new_contig(central_kmer_string.c_str());

		typename _set_T::value_type nkmer;


		for (int direction=0;direction<2;direction++){

			//std::cerr << "direction " << direction << std::endl;

			if (direction==0){
				// forward
			}
			else{
				// reverse
				reverse_complement_in_place(central_kmer_string);
			}

			strncpy(kmer_str,central_kmer_string.c_str(),k);
			kmer_str[k]='\0';

			bool extending = true;

			//std::cerr << "central k-mer: " << central_kmer_string << std::endl;


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
						//std::cerr << "extending " << c << std::endl;
						//debug_print_kmer_set(set,k);
						//std::cerr << std::string(contig.l_ext) << c << std::endl;
						//std::cerr << std::endl;
						if(direction==0){
							contig.r_extend(c);
						}
						else{
							contig.l_extend(c);
						}
						set.erase(nkmer);

						if(!contig.is_full()){
							extending=true;
						}
						break;
					}
				}
			}
		}

		//std::cerr << "====================" << std::endl;

		std::stringstream ss;
		ss<<"contig_"<<contig_id;
		const std::string contig_name(ss.str());
		contig.print_to_fasta(file,contig_name.c_str());
		contig_id++;
	}

	fclose(file);

	std::cout << "   assembly finished (" << contig_id << " contigs)" << std::endl;

	return 0;

}


int main (int argc, char* argv[])
{
	int32_t k=22;

	std::string intersection_fn;
	std::vector<std::string> input_fns;
	std::vector<std::string> output_fns;

	if (argc<2){
		print_help(k);
		exit(1);
	}

	int c;
	while ((c = getopt(argc, (char *const *)argv, "hi:o:x:k:")) >= 0) {
		switch (c) {
			case 'h': {
				print_help(k);
				exit(0);
				break;
			}
			case 'i': {
				input_fns.push_back(std::string(optarg));
				break;
			}
			case 'o': {
				output_fns.push_back(std::string(optarg));
				break;
			}
			case 'x': {
				intersection_fn=std::string(optarg);
				break;
			}
			case 'k': {
				k = atoi(optarg);
				break;
			}
			case '?': {
				std::cerr<<"Unknown error"<<std::endl;
				exit(1);
				break;
			}
		}
	}

	if(k < 2 || max_allowed_kmer_length<k){
		fprintf(stderr, "K-mer size must satisfy 1 < k <= %" PRId32 ".\n", max_allowed_kmer_length);
		return EXIT_FAILURE;
	}

	std::vector< std::unordered_set<nkmer_t> > full_sets(output_fns.size());

	std::cout << "=====================" << std::endl;
	std::cout << "1) Loading references" << std::endl;
	std::cout << "=====================" << std::endl;

	for(int32_t i=0;i<static_cast<int32_t>(input_fns.size());i++){
		kmers_from_fasta(input_fns[i],full_sets[i],k);
		//debug_print_kmer_set(full_sets[i],k);
	}


	std::cout << "===============" << std::endl;
	std::cout << "2) Intersecting" << std::endl;
	std::cout << "===============" << std::endl;


	std::unordered_set<nkmer_t> intersection;

	std::cout << "Computing intersection" << std::endl;

	find_intersection(full_sets, intersection);
	std::cout << "   intersection size: " <<  intersection.size() << std::endl;

	std::cout << "Removing this intersection from all kmer sets" << std::endl;

	int32_t intersection_size  = intersection.size();

	std::vector<int32_t> old_sizes;	
	for (auto const &x : full_sets){
		old_sizes.insert(old_sizes.end(),x.size());
	}

	remove_subset(full_sets, intersection);

	std::vector<int32_t> new_sizes;	
	for (auto const &x : full_sets){
		new_sizes.insert(new_sizes.end(),x.size());
	}

	//for (auto a=old_sizes.begin(),auto b=new_sizes.begin();a<old_sizes.end() && b<new_sizes.end();++a,++b){
	for (int32_t i=0;i<static_cast<int32_t>(old_sizes.size());i++){
		assert(old_sizes[i]==new_sizes[i]+intersection_size);
		std::cout << old_sizes[i] << " " << new_sizes[i] << " ...inter:" << intersection_size << std::endl;
	}


	std::cout << "=============" << std::endl;
	std::cout << "3) Assembling" << std::endl;
	std::cout << "=============" << std::endl;


	for(int32_t i=0;i<static_cast<int32_t>(input_fns.size());i++){
		assemble(output_fns[i],full_sets[i],k);
	}

	assemble(intersection_fn,intersection,k);


	return 0;
}
