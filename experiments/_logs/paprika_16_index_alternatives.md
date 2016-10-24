Table of Contents
=================

* [01_1from10_bacteria_bin_k31](#01_1from10_bacteria_bin_k31)
  * [01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log](#01_1from10_bacteria_bin_k3111_kmer_propagationlog)
  * [01_1from10_bacteria_bin_k31/1.2_merging_fasta.log](#01_1from10_bacteria_bin_k3112_merging_fastalog)
  * [01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log](#01_1from10_bacteria_bin_k3121_bwa_fa2paclog)
  * [01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log](#01_1from10_bacteria_bin_k3122_bwa_pac2bwtgenlog)
  * [01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log](#01_1from10_bacteria_bin_k3123_bwa_bwtupdatelog)
  * [01_1from10_bacteria_bin_k31/2.5_klcp_sa.log](#01_1from10_bacteria_bin_k3125_klcp_salog)
  * [01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log](#01_1from10_bacteria_bin_k3131a_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log](#01_1from10_bacteria_bin_k3131b_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log](#01_1from10_bacteria_bin_k3132a_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log](#01_1from10_bacteria_bin_k3132b_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/4.1_read_assignment.log](#01_1from10_bacteria_bin_k3141_read_assignmentlog)
  * [01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log](#01_1from10_bacteria_bin_k3142_read_assignment_simlcalog)
  * [01_1from10_bacteria_bin_k31/5.1_contigs_stats.log](#01_1from10_bacteria_bin_k3151_contigs_statslog)
  * [01_1from10_bacteria_bin_k31/5.2_index_size.log](#01_1from10_bacteria_bin_k3152_index_sizelog)
* [02_bacteria_orig_k31](#02_bacteria_orig_k31)
  * [02_bacteria_orig_k31/1.1_kmer_propagation.log](#02_bacteria_orig_k3111_kmer_propagationlog)
  * [02_bacteria_orig_k31/1.2_merging_fasta.log](#02_bacteria_orig_k3112_merging_fastalog)
  * [02_bacteria_orig_k31/2.1_bwa_fa2pac.log](#02_bacteria_orig_k3121_bwa_fa2paclog)
  * [02_bacteria_orig_k31/2.2_bwa_pac2bwtgen.log](#02_bacteria_orig_k3122_bwa_pac2bwtgenlog)
  * [02_bacteria_orig_k31/2.3_bwa_bwtupdate.log](#02_bacteria_orig_k3123_bwa_bwtupdatelog)
  * [02_bacteria_orig_k31/2.5_klcp_sa.log](#02_bacteria_orig_k3125_klcp_salog)
  * [02_bacteria_orig_k31/3.1a_matching_rolling.log](#02_bacteria_orig_k3131a_matching_rollinglog)
  * [02_bacteria_orig_k31/3.1b_matching_rolling.log](#02_bacteria_orig_k3131b_matching_rollinglog)
  * [02_bacteria_orig_k31/3.2a_matching_restarted.log](#02_bacteria_orig_k3132a_matching_restartedlog)
  * [02_bacteria_orig_k31/3.2b_matching_restarted.log](#02_bacteria_orig_k3132b_matching_restartedlog)
  * [02_bacteria_orig_k31/4.1_read_assignment.log](#02_bacteria_orig_k3141_read_assignmentlog)
  * [02_bacteria_orig_k31/4.2_read_assignment_simlca.log](#02_bacteria_orig_k3142_read_assignment_simlcalog)
  * [02_bacteria_orig_k31/5.1_contigs_stats.log](#02_bacteria_orig_k3151_contigs_statslog)
  * [02_bacteria_orig_k31/5.2_index_size.log](#02_bacteria_orig_k3152_index_sizelog)
* [03_hmp_orig_k31](#03_hmp_orig_k31)
  * [03_hmp_orig_k31/1.1_kmer_propagation.log](#03_hmp_orig_k3111_kmer_propagationlog)
  * [03_hmp_orig_k31/1.2_merging_fasta.log](#03_hmp_orig_k3112_merging_fastalog)
  * [03_hmp_orig_k31/2.1_bwa_fa2pac.log](#03_hmp_orig_k3121_bwa_fa2paclog)
  * [03_hmp_orig_k31/2.2_bwa_pac2bwtgen.log](#03_hmp_orig_k3122_bwa_pac2bwtgenlog)
  * [03_hmp_orig_k31/2.3_bwa_bwtupdate.log](#03_hmp_orig_k3123_bwa_bwtupdatelog)
  * [03_hmp_orig_k31/2.5_klcp_sa.log](#03_hmp_orig_k3125_klcp_salog)
  * [03_hmp_orig_k31/3.1a_matching_rolling.log](#03_hmp_orig_k3131a_matching_rollinglog)
  * [03_hmp_orig_k31/3.1b_matching_rolling.log](#03_hmp_orig_k3131b_matching_rollinglog)
  * [03_hmp_orig_k31/3.2a_matching_restarted.log](#03_hmp_orig_k3132a_matching_restartedlog)
  * [03_hmp_orig_k31/3.2b_matching_restarted.log](#03_hmp_orig_k3132b_matching_restartedlog)
  * [03_hmp_orig_k31/4.1_read_assignment.log](#03_hmp_orig_k3141_read_assignmentlog)
  * [03_hmp_orig_k31/4.2_read_assignment_simlca.log](#03_hmp_orig_k3142_read_assignment_simlcalog)
  * [03_hmp_orig_k31/5.1_contigs_stats.log](#03_hmp_orig_k3151_contigs_statslog)
  * [03_hmp_orig_k31/5.2_index_size.log](#03_hmp_orig_k3152_index_sizelog)
* [04_bacteria_orig_k31_lowcomp_masked](#04_bacteria_orig_k31_lowcomp_masked)
  * [04_bacteria_orig_k31_lowcomp_masked/1.1_kmer_propagation.log](#04_bacteria_orig_k31_lowcomp_masked11_kmer_propagationlog)
  * [04_bacteria_orig_k31_lowcomp_masked/1.2_merging_fasta.log](#04_bacteria_orig_k31_lowcomp_masked12_merging_fastalog)
  * [04_bacteria_orig_k31_lowcomp_masked/2.1_bwa_fa2pac.log](#04_bacteria_orig_k31_lowcomp_masked21_bwa_fa2paclog)
  * [04_bacteria_orig_k31_lowcomp_masked/2.2_bwa_pac2bwtgen.log](#04_bacteria_orig_k31_lowcomp_masked22_bwa_pac2bwtgenlog)
  * [04_bacteria_orig_k31_lowcomp_masked/2.3_bwa_bwtupdate.log](#04_bacteria_orig_k31_lowcomp_masked23_bwa_bwtupdatelog)
  * [04_bacteria_orig_k31_lowcomp_masked/2.5_klcp_sa.log](#04_bacteria_orig_k31_lowcomp_masked25_klcp_salog)
  * [04_bacteria_orig_k31_lowcomp_masked/3.1a_matching_rolling.log](#04_bacteria_orig_k31_lowcomp_masked31a_matching_rollinglog)
  * [04_bacteria_orig_k31_lowcomp_masked/3.1b_matching_rolling.log](#04_bacteria_orig_k31_lowcomp_masked31b_matching_rollinglog)
  * [04_bacteria_orig_k31_lowcomp_masked/3.2a_matching_restarted.log](#04_bacteria_orig_k31_lowcomp_masked32a_matching_restartedlog)
  * [04_bacteria_orig_k31_lowcomp_masked/3.2b_matching_restarted.log](#04_bacteria_orig_k31_lowcomp_masked32b_matching_restartedlog)
  * [04_bacteria_orig_k31_lowcomp_masked/4.1_read_assignment.log](#04_bacteria_orig_k31_lowcomp_masked41_read_assignmentlog)
  * [04_bacteria_orig_k31_lowcomp_masked/4.2_read_assignment_simlca.log](#04_bacteria_orig_k31_lowcomp_masked42_read_assignment_simlcalog)
  * [04_bacteria_orig_k31_lowcomp_masked/5.1_contigs_stats.log](#04_bacteria_orig_k31_lowcomp_masked51_contigs_statslog)
  * [04_bacteria_orig_k31_lowcomp_masked/5.2_index_size.log](#04_bacteria_orig_k31_lowcomp_masked52_index_sizelog)
* [05_bacteria_orig_k32](#05_bacteria_orig_k32)
  * [05_bacteria_orig_k32/1.1_kmer_propagation.log](#05_bacteria_orig_k3211_kmer_propagationlog)
  * [05_bacteria_orig_k32/1.2_merging_fasta.log](#05_bacteria_orig_k3212_merging_fastalog)
  * [05_bacteria_orig_k32/2.1_bwa_fa2pac.log](#05_bacteria_orig_k3221_bwa_fa2paclog)
  * [05_bacteria_orig_k32/2.2_bwa_pac2bwtgen.log](#05_bacteria_orig_k3222_bwa_pac2bwtgenlog)
  * [05_bacteria_orig_k32/2.3_bwa_bwtupdate.log](#05_bacteria_orig_k3223_bwa_bwtupdatelog)
  * [05_bacteria_orig_k32/2.5_klcp_sa.log](#05_bacteria_orig_k3225_klcp_salog)
  * [05_bacteria_orig_k32/3.1a_matching_rolling.log](#05_bacteria_orig_k3231a_matching_rollinglog)
  * [05_bacteria_orig_k32/3.1b_matching_rolling.log](#05_bacteria_orig_k3231b_matching_rollinglog)
  * [05_bacteria_orig_k32/3.2a_matching_restarted.log](#05_bacteria_orig_k3232a_matching_restartedlog)
  * [05_bacteria_orig_k32/3.2b_matching_restarted.log](#05_bacteria_orig_k3232b_matching_restartedlog)
  * [05_bacteria_orig_k32/4.1_read_assignment.log](#05_bacteria_orig_k3241_read_assignmentlog)
  * [05_bacteria_orig_k32/4.2_read_assignment_simlca.log](#05_bacteria_orig_k3242_read_assignment_simlcalog)
  * [05_bacteria_orig_k32/5.1_contigs_stats.log](#05_bacteria_orig_k3251_contigs_statslog)
  * [05_bacteria_orig_k32/5.2_index_size.log](#05_bacteria_orig_k3252_index_sizelog)
* [06_bacteria_orig_k31_nondel](#06_bacteria_orig_k31_nondel)
  * [06_bacteria_orig_k31_nondel/1.1_kmer_propagation.log](#06_bacteria_orig_k31_nondel11_kmer_propagationlog)
  * [06_bacteria_orig_k31_nondel/1.2_merging_fasta.log](#06_bacteria_orig_k31_nondel12_merging_fastalog)
  * [06_bacteria_orig_k31_nondel/2.1_bwa_fa2pac.log](#06_bacteria_orig_k31_nondel21_bwa_fa2paclog)
  * [06_bacteria_orig_k31_nondel/2.2_bwa_pac2bwtgen.log](#06_bacteria_orig_k31_nondel22_bwa_pac2bwtgenlog)
  * [06_bacteria_orig_k31_nondel/2.3_bwa_bwtupdate.log](#06_bacteria_orig_k31_nondel23_bwa_bwtupdatelog)
  * [06_bacteria_orig_k31_nondel/2.5_klcp_sa.log](#06_bacteria_orig_k31_nondel25_klcp_salog)
  * [06_bacteria_orig_k31_nondel/3.1a_matching_rolling.log](#06_bacteria_orig_k31_nondel31a_matching_rollinglog)
  * [06_bacteria_orig_k31_nondel/3.1b_matching_rolling.log](#06_bacteria_orig_k31_nondel31b_matching_rollinglog)
  * [06_bacteria_orig_k31_nondel/3.2a_matching_restarted.log](#06_bacteria_orig_k31_nondel32a_matching_restartedlog)
  * [06_bacteria_orig_k31_nondel/3.2b_matching_restarted.log](#06_bacteria_orig_k31_nondel32b_matching_restartedlog)
  * [06_bacteria_orig_k31_nondel/4.1_read_assignment.log](#06_bacteria_orig_k31_nondel41_read_assignmentlog)
  * [06_bacteria_orig_k31_nondel/4.2_read_assignment_simlca.log](#06_bacteria_orig_k31_nondel42_read_assignment_simlcalog)
  * [06_bacteria_orig_k31_nondel/5.1_contigs_stats.log](#06_bacteria_orig_k31_nondel51_contigs_statslog)
  * [06_bacteria_orig_k31_nondel/5.2_index_size.log](#06_bacteria_orig_k31_nondel52_index_sizelog)
* [07_bacteria_orig_k31_nonprop](#07_bacteria_orig_k31_nonprop)
  * [07_bacteria_orig_k31_nonprop/1.1_kmer_propagation.log](#07_bacteria_orig_k31_nonprop11_kmer_propagationlog)
  * [07_bacteria_orig_k31_nonprop/1.2_merging_fasta.log](#07_bacteria_orig_k31_nonprop12_merging_fastalog)
  * [07_bacteria_orig_k31_nonprop/2.1_bwa_fa2pac.log](#07_bacteria_orig_k31_nonprop21_bwa_fa2paclog)
  * [07_bacteria_orig_k31_nonprop/2.2_bwa_pac2bwtgen.log](#07_bacteria_orig_k31_nonprop22_bwa_pac2bwtgenlog)
  * [07_bacteria_orig_k31_nonprop/2.3_bwa_bwtupdate.log](#07_bacteria_orig_k31_nonprop23_bwa_bwtupdatelog)
  * [07_bacteria_orig_k31_nonprop/2.5_klcp_sa.log](#07_bacteria_orig_k31_nonprop25_klcp_salog)
  * [07_bacteria_orig_k31_nonprop/3.1a_matching_rolling.log](#07_bacteria_orig_k31_nonprop31a_matching_rollinglog)
  * [07_bacteria_orig_k31_nonprop/3.1b_matching_rolling.log](#07_bacteria_orig_k31_nonprop31b_matching_rollinglog)
  * [07_bacteria_orig_k31_nonprop/3.2a_matching_restarted.log](#07_bacteria_orig_k31_nonprop32a_matching_restartedlog)
  * [07_bacteria_orig_k31_nonprop/3.2b_matching_restarted.log](#07_bacteria_orig_k31_nonprop32b_matching_restartedlog)
  * [07_bacteria_orig_k31_nonprop/4.1_read_assignment.log](#07_bacteria_orig_k31_nonprop41_read_assignmentlog)
  * [07_bacteria_orig_k31_nonprop/4.2_read_assignment_simlca.log](#07_bacteria_orig_k31_nonprop42_read_assignment_simlcalog)
  * [07_bacteria_orig_k31_nonprop/5.1_contigs_stats.log](#07_bacteria_orig_k31_nonprop51_contigs_statslog)
  * [07_bacteria_orig_k31_nonprop/5.2_index_size.log](#07_bacteria_orig_k31_nonprop52_index_sizelog)
* [08_bacteria_orig_k31_nonprop_reasm](#08_bacteria_orig_k31_nonprop_reasm)
  * [08_bacteria_orig_k31_nonprop_reasm/1.1_kmer_propagation.log](#08_bacteria_orig_k31_nonprop_reasm11_kmer_propagationlog)
  * [08_bacteria_orig_k31_nonprop_reasm/1.2_merging_fasta.log](#08_bacteria_orig_k31_nonprop_reasm12_merging_fastalog)
  * [08_bacteria_orig_k31_nonprop_reasm/2.1_bwa_fa2pac.log](#08_bacteria_orig_k31_nonprop_reasm21_bwa_fa2paclog)
  * [08_bacteria_orig_k31_nonprop_reasm/2.2_bwa_pac2bwtgen.log](#08_bacteria_orig_k31_nonprop_reasm22_bwa_pac2bwtgenlog)
  * [08_bacteria_orig_k31_nonprop_reasm/2.3_bwa_bwtupdate.log](#08_bacteria_orig_k31_nonprop_reasm23_bwa_bwtupdatelog)
  * [08_bacteria_orig_k31_nonprop_reasm/2.5_klcp_sa.log](#08_bacteria_orig_k31_nonprop_reasm25_klcp_salog)
  * [08_bacteria_orig_k31_nonprop_reasm/3.1a_matching_rolling.log](#08_bacteria_orig_k31_nonprop_reasm31a_matching_rollinglog)
  * [08_bacteria_orig_k31_nonprop_reasm/3.1b_matching_rolling.log](#08_bacteria_orig_k31_nonprop_reasm31b_matching_rollinglog)
  * [08_bacteria_orig_k31_nonprop_reasm/3.2a_matching_restarted.log](#08_bacteria_orig_k31_nonprop_reasm32a_matching_restartedlog)
  * [08_bacteria_orig_k31_nonprop_reasm/3.2b_matching_restarted.log](#08_bacteria_orig_k31_nonprop_reasm32b_matching_restartedlog)
  * [08_bacteria_orig_k31_nonprop_reasm/4.1_read_assignment.log](#08_bacteria_orig_k31_nonprop_reasm41_read_assignmentlog)
  * [08_bacteria_orig_k31_nonprop_reasm/4.2_read_assignment_simlca.log](#08_bacteria_orig_k31_nonprop_reasm42_read_assignment_simlcalog)
  * [08_bacteria_orig_k31_nonprop_reasm/5.1_contigs_stats.log](#08_bacteria_orig_k31_nonprop_reasm51_contigs_statslog)
  * [08_bacteria_orig_k31_nonprop_reasm/5.2_index_size.log](#08_bacteria_orig_k31_nonprop_reasm52_index_sizelog)

***
## 01_1from10_bacteria_bin_k31

### 01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log
* Sat Oct 15 22:32:18 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:37:16
System time:   0:01:04
Elapsed time:  0:02:38
CPU usage:        1452%

Memory peak:      1.30 GB
```

* 2002544inputs+4217000outputs (62major+36663544minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/1.2_merging_fasta.log
* Sat Oct 15 22:34:57 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:32
System time:   0:00:01
Elapsed time:  0:00:35
CPU usage:          96%

Memory peak:      0.01 GB
```

* 32inputs+1826744outputs (0major+7283minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log
* Sat Oct 15 22:35:33 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:07
System time:   0:00:00
Elapsed time:  0:00:12
CPU usage:          64%

Memory peak:      0.62 GB
```

* 864inputs+1037184outputs (3major+120668minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log
* Sat Oct 15 22:35:46 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     0:10:57
System time:   0:00:00
Elapsed time:  0:10:57
CPU usage:         100%

Memory peak:      0.69 GB
```

* 0inputs+858048outputs (0major+21935minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log
* Sat Oct 15 22:46:44 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:03
System time:   0:00:01
Elapsed time:  0:00:11
CPU usage:          43%

Memory peak:      1.23 GB
```

* 48inputs+1716096outputs (0major+2025minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.5_klcp_sa.log
* Sat Oct 15 22:46:55 CEST 2016
* jobs: 24
* ../../bin/exk index -s -k 31 index.fa

```
User time:     0:14:36
System time:   0:00:02
Elapsed time:  0:10:39
CPU usage:         137%

Memory peak:      1.43 GB
```

* 48inputs+1287072outputs (0major+63963minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log
* Sat Oct 15 22:57:35 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:17
System time:   0:00:01
Elapsed time:  0:01:19
CPU usage:          99%

Memory peak:      1.68 GB
```

* 538704inputs+171288outputs (0major+299032minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log
* bwt_loading	0.43s
* sa_loading	0.21s
* bns_loading	1.32s
* klcp_loading	0.10s
* matching_time	76.51s
* reads	1000000
* kmers	70000000
* rpm	784244
* kpm	54897050


### 01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log
* Sat Oct 15 22:58:54 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:57
System time:   0:00:01
Elapsed time:  0:01:59
CPU usage:         100%

Memory peak:      1.48 GB
```

* 0inputs+171288outputs (0major+346472minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.37s
* sa_loading	0.17s
* bns_loading	1.41s
* matching_time	116.76s
* reads	1000000
* kmers	70000000
* rpm	513857
* kpm	35969967


### 01_1from10_bacteria_bin_k31/4.1_read_assignment.log
* Sat Oct 15 22:58:54 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a

```
User time:     0:00:29
System time:   0:00:00
Elapsed time:  0:00:29
CPU usage:         100%

Memory peak:      0.05 GB
```

* 328inputs+0outputs (1major+21801minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log
* Sat Oct 15 22:58:54 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a -t

```
User time:     0:00:40
System time:   0:00:00
Elapsed time:  0:00:40
CPU usage:         100%

Memory peak:      0.05 GB
```

* 56inputs+0outputs (0major+22506minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/5.1_contigs_stats.log
* Number of contigs: 2192947
* Total length: 878639828
* Average length: 400.66623953976085

* Median length: 65
* Number of k-mer occurencies: 812851418


### 01_1from10_bacteria_bin_k31/5.2_index_size.log
* 210M	index.fa.31.bit.klcp
* 88M	index.fa.ann
* 838M	index.fa.bwt
* 419M	index.fa.pac
* 419M	index.fa.sa

***
## 02_bacteria_orig_k31

### 02_bacteria_orig_k31/1.1_kmer_propagation.log
* Sat Oct 15 23:04:08 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     5:06:04
System time:   0:18:16
Elapsed time:  0:20:09
CPU usage:        1609%

Memory peak:     10.17 GB
```

* 23457312inputs+36108504outputs (6major+392448998minor)pagefaults 0swaps


### 02_bacteria_orig_k31/1.2_merging_fasta.log
* Sat Oct 15 23:24:18 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:05:10
System time:   0:00:15
Elapsed time:  0:06:41
CPU usage:          81%

Memory peak:      0.01 GB
```

* 8249944inputs+16816776outputs (7major+8432minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.1_bwa_fa2pac.log
* Sat Oct 15 23:30:59 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:04
System time:   0:00:14
Elapsed time:  0:02:10
CPU usage:          60%

Memory peak:      6.07 GB
```

* 568inputs+9860920outputs (2major+2799317minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Sat Oct 15 23:33:09 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:15:20
System time:   0:00:35
Elapsed time:  2:15:56
CPU usage:         100%

Memory peak:      4.42 GB
```

* 0inputs+7806112outputs (0major+50361340minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.3_bwa_bwtupdate.log
* Sun Oct 16 01:49:05 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:33
System time:   0:00:36
Elapsed time:  0:02:11
CPU usage:          53%

Memory peak:     11.17 GB
```

* 408inputs+15612216outputs (0major+5738358minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.5_klcp_sa.log
* Sun Oct 16 01:51:17 CEST 2016
* jobs: 24
* ../../bin/exk index -s -k 31 index.fa

```
User time:     2:15:06
System time:   0:02:13
Elapsed time:  1:41:21
CPU usage:         136%

Memory peak:     13.03 GB
```

* 376inputs+11709168outputs (0major+159166378minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1a_matching_rolling.log
* Sun Oct 16 03:32:39 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:38
System time:   0:00:18
Elapsed time:  0:01:58
CPU usage:          99%

Memory peak:     14.22 GB
```

* 492640inputs+317960outputs (2major+7083479minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	3.25s
* sa_loading	3.92s
* bns_loading	16.48s
* klcp_loading	4.55s
* matching_time	88.10s
* reads	1000000
* kmers	70000000
* rpm	681031
* kpm	47672162


### 02_bacteria_orig_k31/3.2a_matching_restarted.log
* Sun Oct 16 03:34:37 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:07:55
System time:   0:00:40
Elapsed time:  0:08:37
CPU usage:         100%

Memory peak:     12.37 GB
```

* 0inputs+317960outputs (0major+39968526minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	10.33s
* sa_loading	1.45s
* bns_loading	15.27s
* matching_time	489.12s
* reads	1000000
* kmers	70000000
* rpm	122670
* kpm	8586926


### 02_bacteria_orig_k31/4.1_read_assignment.log
* Sun Oct 16 03:34:37 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:01:50
System time:   0:00:00
Elapsed time:  0:01:56
CPU usage:          95%

Memory peak:      0.07 GB
```

* 38472inputs+0outputs (162major+133747minor)pagefaults 0swaps


### 02_bacteria_orig_k31/4.2_read_assignment_simlca.log
* Sun Oct 16 03:34:37 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:04:08
System time:   0:00:00
Elapsed time:  0:04:14
CPU usage:          98%

Memory peak:      0.07 GB
```

* 38384inputs+0outputs (154major+233710minor)pagefaults 0swaps


### 02_bacteria_orig_k31/5.1_contigs_stats.log
* Number of contigs: 24170781
* Total length: 7993453766
* Average length: 330.70730176240477

* Median length: 70
* Number of k-mer occurencies: 7268330336


### 02_bacteria_orig_k31/5.2_index_size.log
* 1.9G	index.fa.31.bit.klcp
* 1004M	index.fa.ann
* 7.5G	index.fa.bwt
* 3.8G	index.fa.pac
* 3.8G	index.fa.sa

***
## 03_hmp_orig_k31

### 03_hmp_orig_k31/1.1_kmer_propagation.log
* Sun Oct 16 03:53:14 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     3:02:33
System time:   0:09:21
Elapsed time:  0:12:49
CPU usage:        1495%

Memory peak:      6.32 GB
```

* 15537176inputs+22930936outputs (52major+219198608minor)pagefaults 0swaps


### 03_hmp_orig_k31/1.2_merging_fasta.log
* Sun Oct 16 04:06:04 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:03:07
System time:   0:00:10
Elapsed time:  0:03:33
CPU usage:          92%

Memory peak:      0.01 GB
```

* 1233120inputs+11151128outputs (23major+10181minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.1_bwa_fa2pac.log
* Sun Oct 16 04:09:38 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:40
System time:   0:00:05
Elapsed time:  0:01:16
CPU usage:          61%

Memory peak:      3.56 GB
```

* 944inputs+6189472outputs (4major+879640minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.2_bwa_pac2bwtgen.log
* Sun Oct 16 04:10:54 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     1:17:48
System time:   0:00:06
Elapsed time:  1:17:55
CPU usage:         100%

Memory peak:      3.06 GB
```

* 0inputs+5277880outputs (0major+5154946minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.3_bwa_bwtupdate.log
* Sun Oct 16 05:28:49 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:22
System time:   0:00:17
Elapsed time:  0:01:22
CPU usage:          49%

Memory peak:      7.55 GB
```

* 312inputs+10555752outputs (0major+3642252minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.5_klcp_sa.log
* Sun Oct 16 05:30:12 CEST 2016
* jobs: 24
* ../../bin/exk index -s -k 31 index.fa

```
User time:     1:35:04
System time:   0:06:28
Elapsed time:  1:08:59
CPU usage:         147%

Memory peak:      8.81 GB
```

* 256inputs+7916824outputs (0major+386203872minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1a_matching_rolling.log
* Sun Oct 16 06:39:11 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:19
System time:   0:00:24
Elapsed time:  0:02:44
CPU usage:         100%

Memory peak:      9.43 GB
```

* 0inputs+238816outputs (0major+12184606minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1b_matching_rolling.log
* bwt_loading	8.34s
* sa_loading	5.73s
* bns_loading	7.46s
* klcp_loading	2.99s
* matching_time	138.85s
* reads	1000000
* kmers	70000000
* rpm	432132
* kpm	30249229


### 03_hmp_orig_k31/3.2a_matching_restarted.log
* Sun Oct 16 06:41:55 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:59
System time:   0:00:11
Elapsed time:  0:05:11
CPU usage:         100%

Memory peak:      8.17 GB
```

* 0inputs+238816outputs (0major+4602648minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.2b_matching_restarted.log
* bwt_loading	2.89s
* sa_loading	2.52s
* bns_loading	7.36s
* matching_time	297.86s
* reads	1000000
* kmers	70000000
* rpm	201440
* kpm	14100784


### 03_hmp_orig_k31/4.1_read_assignment.log
* Sun Oct 16 06:41:55 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a

```
User time:     0:01:06
System time:   0:00:00
Elapsed time:  0:01:09
CPU usage:          97%

Memory peak:      0.09 GB
```

* 9992inputs+0outputs (46major+94838minor)pagefaults 0swaps


### 03_hmp_orig_k31/4.2_read_assignment_simlca.log
* Sun Oct 16 06:41:55 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a -t

```
User time:     0:02:11
System time:   0:00:00
Elapsed time:  0:02:13
CPU usage:          99%

Memory peak:      0.09 GB
```

* 15464inputs+0outputs (48major+310529minor)pagefaults 0swaps


### 03_hmp_orig_k31/5.1_contigs_stats.log
* Number of contigs: 10761113
* Total length: 5404543228
* Average length: 502.2290192473585

* Median length: 70
* Number of k-mer occurencies: 5081709838


### 03_hmp_orig_k31/5.2_index_size.log
* 1.3G	index.fa.31.bit.klcp
* 446M	index.fa.ann
* 5.1G	index.fa.bwt
* 2.6G	index.fa.pac
* 2.6G	index.fa.sa

***
## 04_bacteria_orig_k31_lowcomp_masked

### 04_bacteria_orig_k31_lowcomp_masked/1.1_kmer_propagation.log
* Sun Oct 16 14:14:18 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     5:03:49
System time:   0:15:09
Elapsed time:  0:19:19
CPU usage:        1651%

Memory peak:     10.08 GB
```

* 30566184inputs+35891456outputs (726major+343527118minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/1.2_merging_fasta.log
* Sun Oct 16 14:33:37 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:05:29
System time:   0:00:15
Elapsed time:  0:07:00
CPU usage:          82%

Memory peak:      0.01 GB
```

* 9280992inputs+16622872outputs (22major+62606minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/2.1_bwa_fa2pac.log
* Sun Oct 16 14:40:38 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:05
System time:   0:00:17
Elapsed time:  0:02:18
CPU usage:          60%

Memory peak:      6.29 GB
```

* 1184inputs+9971856outputs (3major+2951828minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/2.2_bwa_pac2bwtgen.log
* Sun Oct 16 14:42:56 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:07:45
System time:   0:00:24
Elapsed time:  2:08:10
CPU usage:         100%

Memory peak:      4.34 GB
```

* 0inputs+7653312outputs (0major+21089127minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/2.3_bwa_bwtupdate.log
* Sun Oct 16 16:51:07 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:32
System time:   0:00:30
Elapsed time:  0:02:13
CPU usage:          47%

Memory peak:     10.95 GB
```

* 344inputs+15306624outputs (0major+4554974minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/2.5_klcp_sa.log
* Sun Oct 16 16:53:20 CEST 2016
* jobs: 24
* ../../bin/exk index -s -k 31 index.fa

```
User time:     3:34:23
System time:   0:42:11
Elapsed time:  2:47:24
CPU usage:         153%

Memory peak:     12.77 GB
```

* 512inputs+11479968outputs (1major+2091935554minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/3.1a_matching_rolling.log
* Sun Oct 16 19:40:44 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:46
System time:   0:00:34
Elapsed time:  0:02:26
CPU usage:          96%

Memory peak:     14.10 GB
```

* 2293512inputs+314336outputs (0major+14517274minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/3.1b_matching_rolling.log
* bwt_loading	8.51s
* sa_loading	8.70s
* bns_loading	18.62s
* klcp_loading	5.75s
* matching_time	98.44s
* reads	1000000
* kmers	70000000
* rpm	609483
* kpm	42663780


### 04_bacteria_orig_k31_lowcomp_masked/3.2a_matching_restarted.log
* Sun Oct 16 19:43:10 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:08:26
System time:   0:00:24
Elapsed time:  0:08:51
CPU usage:         100%

Memory peak:     12.27 GB
```

* 24inputs+314336outputs (1major+13975849minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/3.2b_matching_restarted.log
* bwt_loading	3.12s
* sa_loading	5.19s
* bns_loading	17.61s
* matching_time	504.20s
* reads	1000000
* kmers	70000000
* rpm	119001
* kpm	8330074


### 04_bacteria_orig_k31_lowcomp_masked/4.1_read_assignment.log
* Sun Oct 16 19:43:10 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:01:49
System time:   0:00:00
Elapsed time:  0:01:56
CPU usage:          94%

Memory peak:      0.07 GB
```

* 38904inputs+0outputs (180major+100569minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/4.2_read_assignment_simlca.log
* Sun Oct 16 19:43:10 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:03:52
System time:   0:00:00
Elapsed time:  0:03:59
CPU usage:          97%

Memory peak:      0.07 GB
```

* 38784inputs+0outputs (168major+242967minor)pagefaults 0swaps


### 04_bacteria_orig_k31_lowcomp_masked/5.1_contigs_stats.log
* Number of contigs: 27254016
* Total length: 7836989407
* Average length: 287.55356300517326

* Median length: 76.0
* Number of k-mer occurencies: 7019368927


### 04_bacteria_orig_k31_lowcomp_masked/5.2_index_size.log
* 1.9G	index.fa.31.bit.klcp
* 1.2G	index.fa.ann
* 7.3G	index.fa.bwt
* 3.7G	index.fa.pac
* 3.7G	index.fa.sa

***
## 05_bacteria_orig_k32

### 05_bacteria_orig_k32/1.1_kmer_propagation.log
* Mon Oct 17 16:52:51 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     5:04:29
System time:   0:17:37
Elapsed time:  0:18:54
CPU usage:        1704%

Memory peak:     10.17 GB
```

* 26833296inputs+36130040outputs (5major+393148592minor)pagefaults 0swaps


### 05_bacteria_orig_k32/1.2_merging_fasta.log
* Mon Oct 17 17:11:45 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:05:20
System time:   0:00:16
Elapsed time:  0:07:01
CPU usage:          80%

Memory peak:      0.01 GB
```

* 9395488inputs+16840120outputs (21major+62172minor)pagefaults 0swaps


### 05_bacteria_orig_k32/2.1_bwa_fa2pac.log
* Mon Oct 17 17:18:47 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:04
System time:   0:00:12
Elapsed time:  0:02:14
CPU usage:          58%

Memory peak:      6.02 GB
```

* 1312inputs+9837136outputs (3major+2021937minor)pagefaults 0swaps


### 05_bacteria_orig_k32/2.2_bwa_pac2bwtgen.log
* Mon Oct 17 17:21:01 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:14:52
System time:   0:00:32
Elapsed time:  2:15:25
CPU usage:         100%

Memory peak:      4.43 GB
```

* 0inputs+7828704outputs (0major+35280620minor)pagefaults 0swaps


### 05_bacteria_orig_k32/2.3_bwa_bwtupdate.log
* Mon Oct 17 19:36:26 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:33
System time:   0:00:24
Elapsed time:  0:02:24
CPU usage:          40%

Memory peak:     11.20 GB
```

* 416inputs+15657400outputs (0major+2277408minor)pagefaults 0swaps


### 05_bacteria_orig_k32/2.5_klcp_sa.log
* Mon Oct 17 19:38:51 CEST 2016
* jobs: 24
* ../../bin/prophyle-index index -s -k 32 index.fa

```
User time:     2:26:57
System time:   0:08:33
Elapsed time:  1:52:55
CPU usage:         138%

Memory peak:     13.07 GB
```

* 336inputs+11743056outputs (1major+503057853minor)pagefaults 0swaps


### 05_bacteria_orig_k32/3.1a_matching_rolling.log
* Mon Oct 17 21:31:46 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.1b_matching_rolling.log -k 32 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:36
System time:   0:00:16
Elapsed time:  0:01:54
CPU usage:          99%

Memory peak:     14.24 GB
```

* 2507800inputs+706248outputs (1major+4275658minor)pagefaults 0swaps


### 05_bacteria_orig_k32/3.1b_matching_rolling.log
* bwt_loading	3.45s
* sa_loading	3.45s
* bns_loading	16.11s
* klcp_loading	3.72s
* matching_time	85.85s
* reads	1000000
* kmers	69000000
* rpm	698879
* kpm	48222640


### 05_bacteria_orig_k32/3.2a_matching_restarted.log
* Mon Oct 17 21:33:41 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.2b_matching_restarted.log -k 32 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:08:23
System time:   0:00:14
Elapsed time:  0:08:38
CPU usage:         100%

Memory peak:     12.37 GB
```

* 144inputs+706248outputs (2major+7590358minor)pagefaults 0swaps


### 05_bacteria_orig_k32/3.2b_matching_restarted.log
* bwt_loading	3.14s
* sa_loading	2.05s
* bns_loading	15.96s
* matching_time	496.16s
* reads	1000000
* kmers	69000000
* rpm	120929
* kpm	8344090


### 05_bacteria_orig_k32/4.1_read_assignment.log
* Mon Oct 17 21:33:41 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 32 -f sam -a

```
User time:     0:01:51
System time:   0:00:00
Elapsed time:  0:01:58
CPU usage:          94%

Memory peak:      0.07 GB
```

* 38616inputs+0outputs (197major+66720minor)pagefaults 0swaps


### 05_bacteria_orig_k32/4.2_read_assignment_simlca.log
* Mon Oct 17 21:33:41 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 32 -f sam -a -t

```
User time:     0:04:31
System time:   0:00:00
Elapsed time:  0:04:38
CPU usage:          98%

Memory peak:      0.07 GB
```

* 41992inputs+0outputs (176major+232326minor)pagefaults 0swaps


### 05_bacteria_orig_k32/5.1_contigs_stats.log
* Number of contigs: 23624544
* Total length: 8016586123
* Average length: 339.3329464052301

* Median length: 72.0
* Number of k-mer occurencies: 7284225259


### 05_bacteria_orig_k32/5.2_index_size.log
* 1.9G	index.fa.32.bit.klcp
* 981M	index.fa.ann
* 7.5G	index.fa.bwt
* 3.8G	index.fa.pac
* 3.8G	index.fa.sa

***
## 06_bacteria_orig_k31_nondel

### 06_bacteria_orig_k31_nondel/1.1_kmer_propagation.log
* Tue Oct 18 11:41:09 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     6:57:12
System time:   0:24:00
Elapsed time:  0:20:29
CPU usage:        2154%

Memory peak:     10.17 GB
```

* 21621248inputs+19096840outputs (32major+506794640minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/1.2_merging_fasta.log
* Tue Oct 18 12:01:38 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py --nondel index

```
User time:     0:06:08
System time:   0:00:16
Elapsed time:  0:06:46
CPU usage:          95%

Memory peak:      0.01 GB
```

* 368inputs+19156264outputs (5major+62158minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/2.1_bwa_fa2pac.log
* Tue Oct 18 12:08:25 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:05
System time:   0:00:13
Elapsed time:  0:02:07
CPU usage:          62%

Memory peak:      5.04 GB
```

* 952inputs+9829864outputs (3major+1517614minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/2.2_bwa_pac2bwtgen.log
* Tue Oct 18 12:10:32 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     3:09:40
System time:   0:00:57
Elapsed time:  3:10:38
CPU usage:         100%

Memory peak:      5.23 GB
```

* 0inputs+9305024outputs (0major+97100480minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/2.3_bwa_bwtupdate.log
* Tue Oct 18 15:21:10 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:39
System time:   0:00:30
Elapsed time:  0:02:47
CPU usage:          42%

Memory peak:     13.31 GB
```

* 424inputs+18610040outputs (0major+3534280minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/2.5_klcp_sa.log
* Tue Oct 18 15:23:58 CEST 2016
* jobs: 24
* ../../bin/prophyle-index index -s -k 31 index.fa

```
User time:     3:21:18
System time:   0:06:56
Elapsed time:  2:32:54
CPU usage:         136%

Memory peak:     15.53 GB
```

* 488inputs+13957536outputs (0major+454150498minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/3.1a_matching_rolling.log
* Tue Oct 18 17:56:53 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:20
System time:   0:00:45
Elapsed time:  0:03:19
CPU usage:          93%

Memory peak:     15.94 GB
```

* 4596696inputs+888424outputs (1major+19236632minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/3.1b_matching_rolling.log
* bwt_loading	14.99s
* sa_loading	11.63s
* bns_loading	4.50s
* klcp_loading	7.39s
* matching_time	146.73s
* reads	1000000
* kmers	70000000
* rpm	408926
* kpm	28624847


### 06_bacteria_orig_k31_nondel/3.2a_matching_restarted.log
* Tue Oct 18 18:00:13 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:11:50
System time:   0:00:16
Elapsed time:  0:12:10
CPU usage:         100%

Memory peak:     13.73 GB
```

* 104inputs+888424outputs (1major+6244731minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/3.2b_matching_restarted.log
* bwt_loading	3.78s
* sa_loading	4.08s
* bns_loading	4.20s
* matching_time	714.67s
* reads	1000000
* kmers	70000000
* rpm	83955
* kpm	5876869


### 06_bacteria_orig_k31_nondel/4.1_read_assignment.log
* Tue Oct 18 18:00:13 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:03:05
System time:   0:00:01
Elapsed time:  0:03:07
CPU usage:         100%

Memory peak:      0.07 GB
```

* 32inputs+0outputs (0major+265725minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/4.2_read_assignment_simlca.log
* Tue Oct 18 18:00:13 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:07:11
System time:   0:00:01
Elapsed time:  0:07:14
CPU usage:         100%

Memory peak:      0.07 GB
```

* 344inputs+0outputs (1major+511755minor)pagefaults 0swaps


### 06_bacteria_orig_k31_nondel/5.1_contigs_stats.log
* Number of contigs: 6204367
* Total length: 9528336564
* Average length: 1535.7467673978667

* Median length: 62
* Number of k-mer occurencies: 9342205554


### 06_bacteria_orig_k31_nondel/5.2_index_size.log
* 2.3G	index.fa.31.bit.klcp
* 257M	index.fa.ann
* 8.9G	index.fa.bwt
* 4.5G	index.fa.pac
* 4.5G	index.fa.sa

***
## 07_bacteria_orig_k31_nonprop

### 07_bacteria_orig_k31_nonprop/1.1_kmer_propagation.log
* Tue Oct 18 20:37:20 CEST 2016
* jobs: 24
* make -f Makefile.generated V=1

```
User time:     0:00:08
System time:   0:00:21
Elapsed time:  0:01:29
CPU usage:          34%

Memory peak:      0.07 GB
```

* 5054800inputs+17937056outputs (0major+1628310minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/1.2_merging_fasta.log
* Tue Oct 18 23:03:06 CEST 2016
* jobs: 1
* ../../bin/create_final_fasta.py --nondel index

```
User time:     0:04:38
System time:   0:00:16
Elapsed time:  0:05:16
CPU usage:          93%

Memory peak:      0.01 GB
```

* 0inputs+17926800outputs (0major+6879minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/2.1_bwa_fa2pac.log
* Tue Oct 18 23:08:22 CEST 2016
* jobs: 1
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:57
System time:   0:00:18
Elapsed time:  0:02:01
CPU usage:          62%

Memory peak:      4.23 GB
```

* 264inputs+8838472outputs (0major+1714321minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/2.2_bwa_pac2bwtgen.log
* Tue Oct 18 23:10:24 CEST 2016
* jobs: 1
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     3:00:25
System time:   0:00:24
Elapsed time:  3:00:50
CPU usage:         100%

Memory peak:      4.97 GB
```

* 0inputs+8836624outputs (0major+30978757minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/2.3_bwa_bwtupdate.log
* Wed Oct 19 02:11:14 CEST 2016
* jobs: 1
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:37
System time:   0:00:30
Elapsed time:  0:02:51
CPU usage:          40%

Memory peak:     12.64 GB
```

* 472inputs+17673248outputs (0major+4606134minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/2.5_klcp_sa.log
* Wed Oct 19 02:14:06 CEST 2016
* jobs: 1
* ../../bin/prophyle-index index -s -k 31 index.fa

```
User time:     3:46:34
System time:   0:09:27
Elapsed time:  2:57:24
CPU usage:         133%

Memory peak:     14.75 GB
```

* 520inputs+13254936outputs (1major+692607598minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/3.1a_matching_rolling.log
* Wed Oct 19 05:12:02 CEST 2016
* jobs: 1
* ../../bin/prophyle-index match -b -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:02
System time:   0:00:51
Elapsed time:  0:04:55
CPU usage:          99%

Memory peak:     14.90 GB
```

* 524016inputs+904392outputs (0major+37075297minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/3.1b_matching_rolling.log
* bwt_loading	11.38s
* sa_loading	5.10s
* bns_loading	0.01s
* klcp_loading	5.46s
* matching_time	271.00s
* reads	1000000
* kmers	70000000
* rpm	221405
* kpm	15498361


### 07_bacteria_orig_k31_nonprop/3.2a_matching_restarted.log
* Wed Oct 19 05:16:58 CEST 2016
* jobs: 1
* ../../bin/prophyle-index match -b -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:12:53
System time:   0:00:49
Elapsed time:  0:13:43
CPU usage:         100%

Memory peak:     12.80 GB
```

* 14752inputs+904392outputs (1major+40356726minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/3.2b_matching_restarted.log
* bwt_loading	16.28s
* sa_loading	10.04s
* bns_loading	0.01s
* matching_time	795.91s
* reads	1000000
* kmers	70000000
* rpm	75386
* kpm	5277007


### 07_bacteria_orig_k31_nonprop/4.1_read_assignment.log
* Wed Oct 19 05:30:41 CEST 2016
* jobs: 1
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:03:14
System time:   0:00:00
Elapsed time:  0:03:18
CPU usage:          99%

Memory peak:      0.07 GB
```

* 80240inputs+0outputs (223major+110252minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/4.2_read_assignment_simlca.log
* Wed Oct 19 05:34:00 CEST 2016
* jobs: 1
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:06:54
System time:   0:00:00
Elapsed time:  0:06:55
CPU usage:         100%

Memory peak:      0.07 GB
```

* 0inputs+0outputs (0major+300930minor)pagefaults 0swaps


### 07_bacteria_orig_k31_nonprop/5.1_contigs_stats.log
* Number of contigs: 4875
* Total length: 9048702447
* Average length: 1856144.0916923077

* Median length: 1207381
* Number of k-mer occurencies: 9048556197


### 07_bacteria_orig_k31_nonprop/5.2_index_size.log
* 2.2G	index.fa.31.bit.klcp
* 576K	index.fa.ann
* 8.5G	index.fa.bwt
* 4.3G	index.fa.pac
* 4.3G	index.fa.sa

***
## 08_bacteria_orig_k31_nonprop_reasm

### 08_bacteria_orig_k31_nonprop_reasm/1.1_kmer_propagation.log
* Wed Oct 19 10:35:29 CEST 2016
* jobs: 24
* make -f Makefile.generated V=1

```
User time:     3:54:19
System time:   0:06:05
Elapsed time:  0:10:32
CPU usage:        2279%

Memory peak:      0.59 GB
```

* 13040568inputs+17472328outputs (25major+189838266minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/1.2_merging_fasta.log
* Wed Oct 19 10:46:02 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py --nondel index

```
User time:     0:05:20
System time:   0:00:15
Elapsed time:  0:06:00
CPU usage:          93%

Memory peak:      0.01 GB
```

* 0inputs+17481544outputs (0major+14223minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/2.1_bwa_fa2pac.log
* Wed Oct 19 10:52:03 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:57
System time:   0:00:14
Elapsed time:  0:02:03
CPU usage:          59%

Memory peak:      4.25 GB
```

* 904inputs+8705488outputs (3major+1781452minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/2.2_bwa_pac2bwtgen.log
* Wed Oct 19 10:54:06 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     3:14:24
System time:   0:01:46
Elapsed time:  3:16:11
CPU usage:         100%

Memory peak:      4.83 GB
```

* 0inputs+8567976outputs (0major+161918088minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/2.3_bwa_bwtupdate.log
* Wed Oct 19 14:10:18 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:36
System time:   0:00:36
Elapsed time:  0:02:13
CPU usage:          55%

Memory peak:     12.26 GB
```

* 424inputs+17135944outputs (1major+5832984minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/2.5_klcp_sa.log
* Wed Oct 19 14:12:31 CEST 2016
* jobs: 24
* ../../bin/prophyle-index index -s -k 31 index.fa

```
User time:     4:06:38
System time:   0:22:02
Elapsed time:  3:21:50
CPU usage:         133%

Memory peak:     14.30 GB
```

* 1216inputs+12851968outputs (5major+1420739563minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/3.1a_matching_rolling.log
* Wed Oct 19 17:34:21 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:56
System time:   0:00:28
Elapsed time:  0:02:27
CPU usage:          98%

Memory peak:     14.52 GB
```

* 1066008inputs+836608outputs (0major+8047911minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/3.1b_matching_rolling.log
* bwt_loading	8.72s
* sa_loading	9.25s
* bns_loading	1.13s
* klcp_loading	5.10s
* matching_time	120.06s
* reads	1000000
* kmers	70000000
* rpm	499737
* kpm	34981605


### 08_bacteria_orig_k31_nonprop_reasm/3.2a_matching_restarted.log
* Wed Oct 19 17:36:49 CEST 2016
* jobs: 24
* ../../bin/prophyle-index match -b -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:39
System time:   0:00:53
Elapsed time:  0:10:33
CPU usage:         100%

Memory peak:     12.48 GB
```

* 744inputs+836608outputs (0major+57822361minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/3.2b_matching_restarted.log
* bwt_loading	11.74s
* sa_loading	9.80s
* bns_loading	1.16s
* matching_time	609.81s
* reads	1000000
* kmers	70000000
* rpm	98391
* kpm	6887342


### 08_bacteria_orig_k31_nonprop_reasm/4.1_read_assignment.log
* Wed Oct 19 17:36:49 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:02:45
System time:   0:00:01
Elapsed time:  0:02:48
CPU usage:          99%

Memory peak:      0.07 GB
```

* 10800inputs+0outputs (66major+499991minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/4.2_read_assignment_simlca.log
* Wed Oct 19 17:36:49 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:06:24
System time:   0:00:00
Elapsed time:  0:06:27
CPU usage:         100%

Memory peak:      0.07 GB
```

* 5976inputs+0outputs (54major+420276minor)pagefaults 0swaps


### 08_bacteria_orig_k31_nonprop_reasm/5.1_contigs_stats.log
* Number of contigs: 1646440
* Total length: 8773601662
* Average length: 5328.831698695367

* Median length: 91.0
* Number of k-mer occurencies: 8724208462


### 08_bacteria_orig_k31_nonprop_reasm/5.2_index_size.log
* 2.1G	index.fa.31.bit.klcp
* 68M	index.fa.ann
* 8.2G	index.fa.bwt
* 4.1G	index.fa.pac
* 4.1G	index.fa.sa

