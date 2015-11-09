#include <cstdio>
#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <unordered_set>
#include <zlib.h>
#include "kseq.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

namespace logging = boost::log;

KSEQ_INIT(gzFile, gzread)

using namespace std;

void log_init()
{
    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::warning
        //logging::trivial::severity >= logging::trivial::debug
    );
}

char nucl_complement(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
	return 'x';
}

int compare_kmers(const char *kmer1,const char *kmer2,int length){
	for (int i=0;i<length;i++){
		if (kmer1[i]<kmer2[i]){
			return i+1;
		}else{
			if (kmer1[i]>kmer2[i]){
				return -i-1;
			}
		}
	}
	return 0;
}

bool is_kmer_valid(const char *kmer,int length){
	for (int i=0;i<length;i++){
		if (kmer[i]!='A' && kmer[i]!='C' && kmer[i]!='G' && kmer[i]!='T'){
			return false;
		}
	}
	return true;
}

void reverse_compl(const char *kmer, char *reverse_kmer, int length){
	char buffer[length];
	memcpy(buffer,kmer,length);
	for (int i=0;i<length;i++){
		reverse_kmer[length-i-1]=nucl_complement(buffer[i]);
	}
}

void reverse_compl(string &kmer){
	int l=kmer.length();
	char buffer[l];
	reverse_compl(kmer.c_str(),buffer,l);
	kmer=string(buffer,l);	
}

string reverse_compl(const char *kmer, int length){
	char buffer[length];
	reverse_compl(kmer,buffer,length);
	return string(buffer,length);
}

string canonical_kmer(const char *kmer, int length){
	char buffer[length];
	reverse_compl(kmer,buffer,length);
	if (compare_kmers(kmer,buffer,length)>0){
		return string(*kmer,length);
	}
	else {		
		return string(buffer,length);
	}
}

class KmerSet {
	private:
		int k;
		unordered_set<string> kmers_set;
	public:
	
	KmerSet(int k)
	:k(k)
	{}
	
	void insert(string kmer)
	{
	    BOOST_LOG_TRIVIAL(debug) << "Inserting kmer " << kmer;
	    BOOST_LOG_TRIVIAL(trace) << "Set of kmers: " << this->debug_internal_state();
		this->kmers_set.insert(kmer);	
	    BOOST_LOG_TRIVIAL(trace) << "Set of kmers: " << this->debug_internal_state();
	}

	bool pop(string &kmer)
	{
	    BOOST_LOG_TRIVIAL(debug) << "Poping kmer " << kmer;
	    BOOST_LOG_TRIVIAL(trace) << "Set of kmers: " << this->debug_internal_state();
		auto search = kmers_set.find(kmer);
		if (search != kmers_set.end()){
			kmers_set.erase(kmer);
			reverse_compl(kmer);
			kmers_set.erase(kmer);
			reverse_compl(kmer);
		    BOOST_LOG_TRIVIAL(trace) << "Set of kmers: " << this->debug_internal_state();
			return true;
		}
		else{
			return false;
		}
	}	

	string pop()
	{
		string kmer=*(kmers_set.begin());
		//this->debug();
		this->pop(kmer);
		//this->debug();
		return kmer;		
	}
	
	int size(){
		return this->kmers_set.size();
	}
	
	string debug_internal_state(){
		return boost::algorithm::join(this->kmers_set, ", ");
	}
};

void extend_contig_to_right(KmerSet &kmers_set,vector<char> &contig,int k) {
	BOOST_LOG_TRIVIAL(debug) << "Extending contig: " << string(contig.begin(),contig.end());

	string candidate_kmer(k,'x');
	bool extending=true;
	while(extending){
		for(int ci=contig.size()-k+1,i=0;i<k-1;i++,ci++){
			candidate_kmer[i]=contig[ci];
		}
		
		extending=false;
		for(char base: {'A','C','G','T'}){
			candidate_kmer[k-1]=base;
			BOOST_LOG_TRIVIAL(debug) << "Candidate kmer for extension: " << candidate_kmer;
			
			if (kmers_set.pop(candidate_kmer)){
				BOOST_LOG_TRIVIAL(debug) << " ...yes";
				extending=true;
				contig.push_back(base);
				BOOST_LOG_TRIVIAL(debug) << "Current state of contig: " << string(contig.begin(),contig.end());
				break;
			} else {
				BOOST_LOG_TRIVIAL(debug) << " ...no";
			}
		}
	}
}

int main(){
	log_init();
	
	int k=2;
	
	gzFile fp;
	kseq_t *seq;
	int l;
	
	KmerSet kmers_set(k);
	
	fp = gzopen("tests/Borrelia_garinii.fa", "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
	    BOOST_LOG_TRIVIAL(info) << "Loading sequence: " << seq->name.s;
		int positions=strlen(seq->seq.s)-k+1;
		for (int i=0;i<positions;i++){
			if (is_kmer_valid(&(seq->seq.s[i]),k)){
				for(int ii=0;ii<k;ii++){
					kmers_set.insert(string(&(seq->seq.s[i]),k));
					kmers_set.insert(reverse_compl(&(seq->seq.s[i]),k));
				}
			}
		}
	}
	
	
	
	while (kmers_set.size() > 0){
		// initial kmer
		string initial_kmer(kmers_set.pop());
		vector<char> contig(initial_kmer.begin(),initial_kmer.end());		
	    BOOST_LOG_TRIVIAL(debug) << "Initial kmer: " << string(contig.begin(),contig.end());
		
		//extending
		extend_contig_to_right(kmers_set,contig,k);

		//printing
		printf(">");
		int i=0;
		for(char x: contig){
			if (i%40==0){
				printf("\n");
			}
			printf("%c",x);
			i++;
		}
		printf("\n");
	}

	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
