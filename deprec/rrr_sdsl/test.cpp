#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <time.h>

using namespace std;
using namespace sdsl;

template<size_t size>
void check_rrr_vector(bit_vector &b, int M, vector<int> random_indexes) {
  rrr_vector<size> rrrb(b);
  //rrr_vector<127>::rank_1_type rank_rrrb(&rrrb);

  clock_t t = clock();
  int sum_rrrb = 0;
  for(int i = 0; i < M; ++i) {
    sum_rrrb += rrrb[random_indexes[i]];;
  }

  float time = (float)(clock() - t) / CLOCKS_PER_SEC;
  cout << endl;
  cout << "block size " << size << endl;
  cout << "sum rrrb = " << sum_rrrb << endl;
  cout << "sum rrrb counting time " << time << " sec" << endl;
  cout << "one element access time in rrrvector " << time / M << endl;
  cout<< "size of rrrb in MB: " << size_in_mega_bytes(rrrb)<< endl;

  // rrr_vector<i>::select_1_type select_rrrb(&rrrb);
  // cout<<"position of first one in b: "<<select_rrrb(1)<<endl;
  // cout << "size of select_rrrb in MB = " << size_in_mega_bytes(select_rrrb) << endl;
  cout<<endl;
}

void check_compressed_vectors(bit_vector &b) {
  int M = 10000000;
  vector<int> random_indexes(M, 0);
  for(int i = 0; i < M; ++i) {
    random_indexes[i] = rand() % b.size();
  }

  clock_t t = clock();
  int sum_indexes = 0;
  for(int i = 0; i < M; ++i) {
    sum_indexes += random_indexes[i];
  }
  float time = (float)(clock() - t) / CLOCKS_PER_SEC;
  cout << "sum indexes = " << sum_indexes << endl;
  cout << "sum indexes counting time " << time << " sec" << endl;
  cout << "one element access time in usual vector " << time / M << endl;

  t = clock();
  int sum_b = 0;
  for(int i = 0; i < M; ++i) {
    sum_b += b[random_indexes[i]];
  }
  time = (float)(clock() - t) / CLOCKS_PER_SEC;
  cout << "sum b = " << sum_b << endl;
  cout << "sum b counting time " << time << " sec" << endl;
  cout << "one element access time in bitvector " << time / M << endl;

  //rank_support_v<> rb(&b);

  cout<< "size of b in MB: " << size_in_mega_bytes(b)<< endl;
  //cout<< "size of rb in MB: " << size_in_mega_bytes(rb)<< endl;

  check_rrr_vector<7>(b, M, random_indexes);
  check_rrr_vector<15>(b, M, random_indexes);
  check_rrr_vector<31>(b, M, random_indexes);
  check_rrr_vector<63>(b, M, random_indexes);
  check_rrr_vector<127>(b, M, random_indexes);
  check_rrr_vector<255>(b, M, random_indexes);
}

int main()
{
    srand(time(NULL));

    uint64_t N = 10000000000;
    bit_vector b(N, 0);
    for(uint64_t i = 0; i < N; ++i) {
      if (rand() % 100 == 3) {
		    b[i] = 1;
	    }
    }
    check_compressed_vectors(b);
    return 0;

    FILE *fp;
	  fp = fopen("klcp_big.bin", "rb");
    uint64_t seq_len;
	  fread(&seq_len, sizeof(uint64_t), 1, fp);
    uint64_t capacity = (seq_len + 7) / 8;
    cout << "seq_len " << seq_len << endl;
    cout << "capacity " << capacity << endl;
    bit_vector klcp(seq_len, 0);
    cout << "bitvector created" << std::endl;
    FILE *log;
    uint64_t ones = 0;
    //log = fopen("klcp_log.txt", "w");
    for(uint64_t i = 0; i < capacity; ++i) {
      //fprintf(log, "%llu ", i);
      int8_t element;
      int a = fread(&element, sizeof(int8_t), 1, fp);
      if (a != sizeof(int8_t)) {
        cerr << "can not read" << std::endl;
      }
      //cout << (int)element << endl;
      for(int j = 0; j < 8; ++j) {
        klcp[i * 8 + j] = ((element & (1 << (j % 8))) > 0);
        ones += klcp[i * 8 + j];
      }
    }
    cout << "ones count: " << ones << endl;
    // for(int i = 0; i < seq_len; ++i) {
    //   cout << klcp[i] << endl;
    // }
    check_compressed_vectors(klcp);

    //
    // t = clock();
    // int sum_rank_rrrb = 0;
    // for(int i = 0; i < N; ++i) {
    //   sum_rank_rrrb += rank_rrrb(i);
    // }
    //
    // cout << "sum rank rrrb counting time " << (float)(clock() - t) / CLOCKS_PER_SEC << " sec" << endl;
    // cout << "sum b = " << sum_b << ", sum rank rrrb = " << sum_rank_rrrb << endl;
    // cout << "sum rrrb = " << sum_rrrb << endl;

    //cout<< "size of rank_rrrb in MB: " << size_in_mega_bytes(rank_rrrb)<< endl;
}
