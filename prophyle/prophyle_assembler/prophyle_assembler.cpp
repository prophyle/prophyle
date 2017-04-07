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
//const int32_t default_k=31;

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


void print_help(){
	std::cerr <<
		"\n" <<
		"Program:  prophyle_assembler (greedy assembler for ProPhyle)\n" <<
		"Contact:  Karel Brinda <karel.brinda@gmail.com>\n" <<
		"\n" <<
		"Usage:    prophyle_assembler [options]\n" <<
		"\n" <<
		"Examples: prophyle_assembler -k 15 -i f1.fa -i f2.fa -x fx.fa\n" <<
		"             - compute intersection of f1 and f2\n" <<
		"          prophyle_assembler -k 15 -i f1.fa -i f2.fa -x fx.fa -o g1.fa -o g2.fa\n" <<
		"             - compute intersection of f1 and f2, and subtract it from them\n" <<
		"          prophyle_assembler -k 15 -i f1.fa -o g1.fa\n" <<
		"             - re-assemble f1 to g1\n" <<
		"\n" <<
		"Command-line parameters:\n" <<
		" -k INT   K-mer size.\n" <<
		" -i FILE  Input FASTA file (can be used multiple times).\n" <<
		" -o FILE  Output FASTA file (if used, must be used as many times as -i).\n" <<
		" -x FILE  Compute intersection, subtract it, save it.\n" <<
		" -s FILE  Output file with k-mer statistics.\n" <<
		//" -k INT   K-mer size. [" << default_k << "]\n" <<
		" -S       Silent mode.\n" <<
		"\n" <<
		"Note that '-' can be used for standard input/output. \n" <<
		std::endl;
}

void test_file(FILE *fo, std::string fn){
	if(fo==nullptr){
		std::cerr << "Error: file '" << fn << "' could not be open." << std::endl;
		exit(1);
	}
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
void debug_print_kmer_set(_set_T &set, int k, bool verbose){
	std::string kmer;
	for(auto x: set){
		decode_kmer(x, k, kmer);
		if(verbose){
			std::cerr << x << " " << kmer << ";  ";
		}
	}
	if(verbose){
		std::cerr<<std::endl;
	}
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
int kmers_from_fasta(const std::string &fasta_fn, _set_T &set, int32_t k, FILE* fstats,bool verbose){

	if (verbose){
		std::cerr << "Loading " << fasta_fn << std::endl;
	}

	set.clear();

	kseq_t *seq;
	int64_t l;

	FILE *instream = nullptr;
	if(fasta_fn=="-"){
		instream = stdin;
	}
	else {
		instream = fopen(fasta_fn.c_str(), "r");
		test_file(instream, fasta_fn);
	}
	gzFile fp = gzdopen(fileno(instream), "r");
	seq = kseq_init(fp);

	typename _set_T::value_type nkmer;

	for(int32_t seqid=0;(l = kseq_read(seq)) >= 0;seqid++) {
		//std::cerr << "kmers from fasta" << std::endl;

		//std::cerr << "starting iterator" << std::endl;
		for(char *kmer=seq->seq.s; kmer < (seq->seq.s) + (seq->seq.l) -k +1 ;kmer++){
			int c = encode_canonical(kmer, k, nkmer);
			if (c==0){
				set.insert(nkmer);
			}
			else{
				//std::cerr << "problem" <<std::endl;
			}
		}
	}

	if(fstats){
		fprintf(fstats,"%s\t%lu\n",fasta_fn.c_str(),set.size());
	}
	//std::cerr << "iterator finished" << std::endl;

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

	//std::cerr << "searching min" << std::endl;

	for(int32_t i=0;i<static_cast<int32_t>(sets.size());i++){
		if (static_cast<int32_t>(sets[i].size())<min){
			min=sets[i].size();
			i_min=i;
			//std::cerr << "new min" << i << std::endl;
		}
	}

	assert(i_min != std::numeric_limits<int32_t>::max() && i_min!=-1);

	/*
		2) Take it as the intersection.
	*/

	//std::cerr << "2" << std::endl;

	intersection.clear();
	//std::cerr << "2.1" << std::endl;
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
int assemble(const std::string &fasta_fn, _set_T &set, int32_t k, FILE* fstats, bool verbose){
	if(fstats){
		fprintf(fstats,"%s\t%lu\n",fasta_fn.c_str(),set.size());
	}

	FILE *file=nullptr;
	if(fasta_fn=="-"){
		file=stdout;
	}
	else{
		file=fopen(fasta_fn.c_str(),"w+");
		test_file(file, fasta_fn);
	}
	char kmer_str[max_allowed_kmer_length+1];
	contig_t contig(k);
	const std::vector<char> nucls = {'A','C','G','T'};

	//int32_t i=0;
	int32_t contig_id=1;
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
		ss<<"c"<<contig_id;
		const std::string contig_name(ss.str());
		contig.print_to_fasta(file,contig_name.c_str());
		contig_id++;
	}

	fclose(file);

	if(verbose){
		std::cerr << "   assembly finished (" << contig_id << " contigs)" << std::endl;
	}

	return 0;

}


int main (int argc, char* argv[])
{
	int32_t k=-1;

	std::string intersection_fn;
	std::vector<std::string> in_fns;
	std::vector<std::string> out_fns;
	std::string stats_fn;
	FILE *fstats=nullptr;

	if (argc<2){
		print_help();
		exit(1);
	}

	bool compute_intersection=false;
	bool compute_output=false;
	bool verbose=true;
	int32_t no_sets=0;

	int c;
	while ((c = getopt(argc, (char *const *)argv, "hSi:o:x:s:k:")) >= 0) {
		switch (c) {
			case 'h': {
				print_help();
				exit(0);
				break;
			}
			case 'i': {
				in_fns.push_back(std::string(optarg));
				no_sets+=1;
				break;
			}
			case 'o': {
				out_fns.push_back(std::string(optarg));
				compute_output=true;
				break;
			}
			case 'x': {
				intersection_fn=std::string(optarg);
				compute_intersection=true;
				break;
			}
			case 's': {
				stats_fn=std::string(optarg);
				if(stats_fn=="-"){
					fstats = stdin;
				}
				else {
					fstats = fopen(stats_fn.c_str(), "w+");
					test_file(fstats, stats_fn);
				}

				break;
			}
			case 'S': {
				verbose=false;

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

	if (k == -1){
		print_help();
		std::cerr << "K-mer size (-k) is required." << std::endl;
		return EXIT_FAILURE;
	}


	if (k <= 0 || max_allowed_kmer_length<k){
		std::cerr << "K-mer size must satisfy 1 <= k <= " << max_allowed_kmer_length << "." << std::endl;
		return EXIT_FAILURE;
	}

	if (compute_output && (static_cast<int32_t>(out_fns.size())!=no_sets)){
		std::cerr << "If -o is used, it must be used as many times as -i (" << no_sets << "!=" << out_fns.size() << ")." << std::endl;
		return EXIT_FAILURE;
	}

	if(fstats){
		fprintf(fstats,"# cmd: %s",argv[0]);

		for (int32_t i=1;i<argc;i++){
			fprintf(fstats," %s",argv[i]);
		}
		fprintf(fstats,"\n");
	}

	std::vector< std::unordered_set<nkmer_t> > full_sets(no_sets);

	if(verbose){
		std::cerr << "=====================" << std::endl;
		std::cerr << "1) Loading references" << std::endl;
		std::cerr << "=====================" << std::endl;
	}


	std::vector<int32_t> in_sizes;
	std::vector<int32_t> out_sizes;

	for(int32_t i=0;i<no_sets;i++){
		kmers_from_fasta(in_fns[i],full_sets[i],k,fstats,verbose);
		//debug_print_kmer_set(full_sets[i],k);
		in_sizes.insert(in_sizes.end(),full_sets[i].size());
	}

	if(verbose){
		std::cerr << "===============" << std::endl;
		std::cerr << "2) Intersecting" << std::endl;
		std::cerr << "===============" << std::endl;
	}


	std::unordered_set<nkmer_t> intersection;

	int32_t intersection_size = 0;

	if(compute_intersection){
		if (verbose){
			std::cerr << "2.1) Computing intersection" << std::endl;
		}

		find_intersection(full_sets, intersection);
		intersection_size  = intersection.size();
		if (verbose){
			std::cerr << "   intersection size: " <<  intersection_size << std::endl;
		}
		if(compute_output){
			if (verbose){
				std::cerr << "2.2) Removing this intersection from all kmer sets" << std::endl;
			}
			remove_subset(full_sets, intersection);
		}
	}

	if(compute_output){
		for (int32_t i=0;i<no_sets;i++){
			out_sizes.insert(out_sizes.end(),full_sets[i].size());
			assert(in_sizes[i]==out_sizes[i]+intersection_size);
			if (verbose){
				std::cerr << in_sizes[i] << " " << out_sizes[i] << " ...inter:" << intersection_size << std::endl;
			}
		}
	}

	if(verbose){
		std::cerr << "=============" << std::endl;
		std::cerr << "3) Assembling" << std::endl;
		std::cerr << "=============" << std::endl;
	}

	if(compute_output){
		for(int32_t i=0;i<static_cast<int32_t>(in_fns.size());i++){
			assemble(out_fns[i], full_sets[i], k, fstats, verbose);
		}
	}
	if(compute_intersection){
		assemble(intersection_fn, intersection, k, fstats, verbose);
	}

	if (fstats){
		fclose(fstats);
	}

	return 0;
}
