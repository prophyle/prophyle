#include <cstdio>
#include <iostream>
#include <cassert>
#include <string>
#include <unordered_set>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

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
		this->kmers_set.insert(kmer);	
	}

	bool pop(string &kmer)
	{
		auto search = kmers_set.find(kmer);
		if (search != kmers_set.end()){
			kmers_set.erase(kmer);
			reverse_compl(kmer);
			kmers_set.erase(kmer);
			reverse_compl(kmer);
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
	
	void debug(){
		printf("Set of k-mers:\n");
		for ( auto it = this->kmers_set.begin(); it != this->kmers_set.end(); ++it ){
			printf("%s ", it->c_str());
		}
		printf("\n\n");
		
	}
};

int main(){
	int k=2;
	printf("pseudoassembler\n");
	
	gzFile fp;
	kseq_t *seq;
	int l;
	
	KmerSet kmers_set(k);
	
	fp = gzopen("tests/Borrelia_garinii.fa", "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		printf("name: %s\n", seq->name.s);
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
	
	string contig;
	
	while (kmers_set.size() > 0){
		// initial kmer
		contig=kmers_set.pop();		
		
		//extending
		string candidate_kmer;

		bool extending=true;
		while(extending){
			candidate_kmer=contig.substr(contig.size()-k,k-1)+"?";

			//printf("candidate for extension %s\n", candidate_kmer.c_str());
			
			extending=false;
			for(char x: {'A','C','G','T'}){
				candidate_kmer[k-1]=x;

				if (kmers_set.pop(candidate_kmer)){
					//printf("yes %c\n", x);
					extending=true;
					contig=contig+string(1,x);
					break;
				} else {
					//printf("no %c\n", x);
				}
			}
		}
	    std::cout << ">\n" << contig << "\n";
		
			
	}

	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
