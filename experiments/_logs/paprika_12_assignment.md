Table of Contents
=================

* [01_1from10_bacteria_bin_k31](#01_1from10_bacteria_bin_k31)
  * [01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log](#01_1from10_bacteria_bin_k3111_kmer_propagationlog)
  * [01_1from10_bacteria_bin_k31/1.2_merging_fasta.log](#01_1from10_bacteria_bin_k3112_merging_fastalog)
  * [01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log](#01_1from10_bacteria_bin_k3121_bwa_fa2paclog)
  * [01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log](#01_1from10_bacteria_bin_k3122_bwa_pac2bwtgenlog)
  * [01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log](#01_1from10_bacteria_bin_k3123_bwa_bwtupdatelog)
  * [01_1from10_bacteria_bin_k31/2.4_bwa_bwt2sa.log](#01_1from10_bacteria_bin_k3124_bwa_bwt2salog)
  * [01_1from10_bacteria_bin_k31/2.5_klcp.log](#01_1from10_bacteria_bin_k3125_klcplog)
  * [01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log](#01_1from10_bacteria_bin_k3131a_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log](#01_1from10_bacteria_bin_k3131b_matching_rollinglog)
  * [01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log](#01_1from10_bacteria_bin_k3132a_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log](#01_1from10_bacteria_bin_k3132b_matching_restartedlog)
  * [01_1from10_bacteria_bin_k31/3.3a_matching_rolling_skipping.log](#01_1from10_bacteria_bin_k3133a_matching_rolling_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.3b_matching_rolling_skipping.log](#01_1from10_bacteria_bin_k3133b_matching_rolling_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.4a_matching_restarted_skipping.log](#01_1from10_bacteria_bin_k3134a_matching_restarted_skippinglog)
  * [01_1from10_bacteria_bin_k31/3.4b_matching_restarted_skipping.log](#01_1from10_bacteria_bin_k3134b_matching_restarted_skippinglog)
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
  * [02_bacteria_orig_k31/2.4_bwa_bwt2sa.log](#02_bacteria_orig_k3124_bwa_bwt2salog)
  * [02_bacteria_orig_k31/2.5_klcp.log](#02_bacteria_orig_k3125_klcplog)
  * [02_bacteria_orig_k31/3.1a_matching_rolling.log](#02_bacteria_orig_k3131a_matching_rollinglog)
  * [02_bacteria_orig_k31/3.1b_matching_rolling.log](#02_bacteria_orig_k3131b_matching_rollinglog)
  * [02_bacteria_orig_k31/3.2a_matching_restarted.log](#02_bacteria_orig_k3132a_matching_restartedlog)
  * [02_bacteria_orig_k31/3.2b_matching_restarted.log](#02_bacteria_orig_k3132b_matching_restartedlog)
  * [02_bacteria_orig_k31/3.3a_matching_rolling_skipping.log](#02_bacteria_orig_k3133a_matching_rolling_skippinglog)
  * [02_bacteria_orig_k31/3.3b_matching_rolling_skipping.log](#02_bacteria_orig_k3133b_matching_rolling_skippinglog)
  * [02_bacteria_orig_k31/3.4a_matching_restarted_skipping.log](#02_bacteria_orig_k3134a_matching_restarted_skippinglog)
  * [02_bacteria_orig_k31/3.4b_matching_restarted_skipping.log](#02_bacteria_orig_k3134b_matching_restarted_skippinglog)
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
  * [03_hmp_orig_k31/2.4_bwa_bwt2sa.log](#03_hmp_orig_k3124_bwa_bwt2salog)
  * [03_hmp_orig_k31/2.5_klcp.log](#03_hmp_orig_k3125_klcplog)
  * [03_hmp_orig_k31/3.1a_matching_rolling.log](#03_hmp_orig_k3131a_matching_rollinglog)
  * [03_hmp_orig_k31/3.1b_matching_rolling.log](#03_hmp_orig_k3131b_matching_rollinglog)
  * [03_hmp_orig_k31/3.2a_matching_restarted.log](#03_hmp_orig_k3132a_matching_restartedlog)
  * [03_hmp_orig_k31/3.2b_matching_restarted.log](#03_hmp_orig_k3132b_matching_restartedlog)
  * [03_hmp_orig_k31/3.3a_matching_rolling_skipping.log](#03_hmp_orig_k3133a_matching_rolling_skippinglog)
  * [03_hmp_orig_k31/3.3b_matching_rolling_skipping.log](#03_hmp_orig_k3133b_matching_rolling_skippinglog)
  * [03_hmp_orig_k31/3.4a_matching_restarted_skipping.log](#03_hmp_orig_k3134a_matching_restarted_skippinglog)
  * [03_hmp_orig_k31/3.4b_matching_restarted_skipping.log](#03_hmp_orig_k3134b_matching_restarted_skippinglog)
  * [03_hmp_orig_k31/4.1_read_assignment.log](#03_hmp_orig_k3141_read_assignmentlog)
  * [03_hmp_orig_k31/4.2_read_assignment_simlca.log](#03_hmp_orig_k3142_read_assignment_simlcalog)
  * [03_hmp_orig_k31/5.1_contigs_stats.log](#03_hmp_orig_k3151_contigs_statslog)
  * [03_hmp_orig_k31/5.2_index_size.log](#03_hmp_orig_k3152_index_sizelog)

***
## 01_1from10_bacteria_bin_k31

### 01_1from10_bacteria_bin_k31/1.1_kmer_propagation.log
* Fri Jul 22 19:35:45 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     0:38:29
System time:   0:01:27
Elapsed time:  0:02:54
CPU usage:        1377%

Memory peak:      1.12 GB
```

* 0inputs+5007784outputs (0major+41134683minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/1.2_merging_fasta.log
* Fri Jul 22 19:38:39 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:00:58
System time:   0:00:01
Elapsed time:  0:01:01
CPU usage:          98%

Memory peak:      0.01 GB
```

* 0inputs+2443112outputs (0major+4543minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.1_bwa_fa2pac.log
* Fri Jul 22 19:39:41 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:00:12
System time:   0:00:02
Elapsed time:  0:00:22
CPU usage:          64%

Memory peak:      1.32 GB
```

* 64inputs+1741464outputs (0major+555955minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul 22 19:40:03 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     0:13:19
System time:   0:00:00
Elapsed time:  0:13:19
CPU usage:         100%

Memory peak:      0.79 GB
```

* 0inputs+1042696outputs (0major+71377minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.3_bwa_bwtupdate.log
* Fri Jul 22 19:53:23 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:04
System time:   0:00:02
Elapsed time:  0:00:15
CPU usage:          45%

Memory peak:      1.49 GB
```

* 64inputs+2085384outputs (0major+334022minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.4_bwa_bwt2sa.log
* Fri Jul 22 19:53:39 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:04:17
System time:   0:00:03
Elapsed time:  0:04:26
CPU usage:          98%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+2757392minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/2.5_klcp.log
* Fri Jul 22 19:58:06 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     0:11:39
System time:   0:00:04
Elapsed time:  0:11:56
CPU usage:          98%

Memory peak:      1.49 GB
```

* 40inputs+1042696outputs (0major+2932504minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1a_matching_rolling.log
* Fri Jul 22 20:10:02 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:23
System time:   0:00:02
Elapsed time:  0:01:30
CPU usage:          94%

Memory peak:      2.75 GB
```

* 538728inputs+173160outputs (0major+327969minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.1b_matching_rolling.log
* bwt_loading	0.54s
* sa_loading	0.27s
* bns_loading	5.26s
* klcp_loading	0.21s
* matching_time	79.24s
* reads	1000000
* kmers	70000000
* rpm	757207
* kpm	53004474


### 01_1from10_bacteria_bin_k31/3.2a_matching_restarted.log
* Fri Jul 22 20:11:33 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:02
System time:   0:00:01
Elapsed time:  0:02:04
CPU usage:         100%

Memory peak:      2.26 GB
```

* 0inputs+173160outputs (0major+411177minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.2b_matching_restarted.log
* bwt_loading	0.41s
* sa_loading	0.20s
* bns_loading	5.26s
* matching_time	118.23s
* reads	1000000
* kmers	70000000
* rpm	507505
* kpm	35525333


### 01_1from10_bacteria_bin_k31/3.3a_matching_rolling_skipping.log
* Fri Jul 22 20:13:38 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:01:26
System time:   0:00:01
Elapsed time:  0:01:28
CPU usage:         100%

Memory peak:      2.75 GB
```

* 0inputs+173160outputs (0major+424662minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	0.43s
* sa_loading	0.19s
* bns_loading	4.99s
* klcp_loading	0.18s
* matching_time	82.19s
* reads	1000000
* kmers	70000000
* rpm	730053
* kpm	51103705


### 01_1from10_bacteria_bin_k31/3.4a_matching_restarted_skipping.log
* Fri Jul 22 20:15:06 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:00
System time:   0:00:01
Elapsed time:  0:02:01
CPU usage:         100%

Memory peak:      2.26 GB
```

* 0inputs+173160outputs (0major+379439minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	0.41s
* sa_loading	0.19s
* bns_loading	4.99s
* matching_time	115.78s
* reads	1000000
* kmers	70000000
* rpm	518221
* kpm	36275497


### 01_1from10_bacteria_bin_k31/4.1_read_assignment.log
* Fri Jul 22 20:11:33 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a

```
User time:     0:00:29
System time:   0:00:00
Elapsed time:  0:00:30
CPU usage:          98%

Memory peak:      0.06 GB
```

* 6880inputs+0outputs (54major+20299minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/4.2_read_assignment_simlca.log
* Fri Jul 22 20:11:33 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/1from10.bacteria.bin.newick -k 31 -f sam -a -t

```
User time:     0:00:46
System time:   0:00:00
Elapsed time:  0:00:46
CPU usage:          98%

Memory peak:      0.05 GB
```

* 5648inputs+0outputs (54major+22371minor)pagefaults 0swaps


### 01_1from10_bacteria_bin_k31/5.1_contigs_stats.log
* Number of contigs: 8495372
* Total length: 1067712578
* Average length: 125.68167444580413

* Median length: 34.0
* Number of k-mer occurencies: 812851418


### 01_1from10_bacteria_bin_k31/5.2_index_size.log
* 510M	index.fa.31.bit.klcp
* 342M	index.fa.ann
* 1019M	index.fa.bwt
* 510M	index.fa.pac
* 510M	index.fa.sa

***
## 02_bacteria_orig_k31

### 02_bacteria_orig_k31/1.1_kmer_propagation.log
* Fri Jul 22 20:17:09 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     4:59:08
System time:   0:16:03
Elapsed time:  0:22:41
CPU usage:        1389%

Memory peak:      9.31 GB
```

* 32611512inputs+43545296outputs (46major+376635908minor)pagefaults 0swaps


### 02_bacteria_orig_k31/1.2_merging_fasta.log
* Fri Jul 22 20:39:50 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:12:08
System time:   0:00:27
Elapsed time:  0:17:37
CPU usage:          71%

Memory peak:      0.01 GB
```

* 16204336inputs+23879936outputs (21major+27143minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.1_bwa_fa2pac.log
* Fri Jul 22 20:57:29 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:02:09
System time:   0:00:40
Elapsed time:  0:04:43
CPU usage:          60%

Memory peak:     13.94 GB
```

* 1688208inputs+18043024outputs (5major+9165965minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.2_bwa_pac2bwtgen.log
* Fri Jul 22 21:02:12 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     2:42:13
System time:   0:00:07
Elapsed time:  2:42:32
CPU usage:         100%

Memory peak:      5.54 GB
```

* 2145672inputs+9887296outputs (2major+422671minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.3_bwa_bwtupdate.log
* Fri Jul 22 23:44:45 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:42
System time:   0:00:28
Elapsed time:  0:02:43
CPU usage:          44%

Memory peak:     14.15 GB
```

* 360inputs+19774584outputs (0major+3352956minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.4_bwa_bwt2sa.log
* Fri Jul 22 23:47:28 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:50:12
System time:   0:01:44
Elapsed time:  0:52:49
CPU usage:          98%

Memory peak:     14.15 GB
```

* 840inputs+9887296outputs (4major+125206730minor)pagefaults 0swaps


### 02_bacteria_orig_k31/2.5_klcp.log
* Sat Jul 23 00:40:17 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     1:59:46
System time:   0:03:38
Elapsed time:  2:04:34
CPU usage:          99%

Memory peak:     14.15 GB
```

* 5372496inputs+9887296outputs (3major+187924523minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1a_matching_rolling.log
* Sat Jul 23 02:44:51 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:14
System time:   0:01:01
Elapsed time:  0:04:39
CPU usage:          91%

Memory peak:     25.84 GB
```

* 18290592inputs+360824outputs (1major+23959885minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.1b_matching_rolling.log
* bwt_loading	8.95s
* sa_loading	11.76s
* bns_loading	65.83s
* klcp_loading	17.63s
* matching_time	149.40s
* reads	1000000
* kmers	70000000
* rpm	401613
* kpm	28112905


### 02_bacteria_orig_k31/3.2a_matching_restarted.log
* Sat Jul 23 02:49:31 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:29
System time:   0:00:16
Elapsed time:  0:10:40
CPU usage:          92%

Memory peak:     21.13 GB
```

* 17584280inputs+360824outputs (1major+2597057minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.2b_matching_restarted.log
* bwt_loading	8.49s
* sa_loading	1.94s
* bns_loading	62.83s
* matching_time	512.14s
* reads	1000000
* kmers	70000000
* rpm	117156
* kpm	8200908


### 02_bacteria_orig_k31/3.3a_matching_rolling_skipping.log
* Sat Jul 23 03:00:12 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:03:10
System time:   0:01:06
Elapsed time:  0:04:35
CPU usage:          93%

Memory peak:     25.84 GB
```

* 10539440inputs+360824outputs (2major+22781379minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	8.83s
* sa_loading	10.65s
* bns_loading	67.63s
* klcp_loading	16.39s
* matching_time	151.12s
* reads	1000000
* kmers	70000000
* rpm	397040
* kpm	27792812


### 02_bacteria_orig_k31/3.4a_matching_restarted_skipping.log
* Sat Jul 23 03:04:47 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:09:04
System time:   0:00:39
Elapsed time:  0:10:45
CPU usage:          90%

Memory peak:     21.12 GB
```

* 20054872inputs+360824outputs (3major+16854823minor)pagefaults 0swaps


### 02_bacteria_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	8.04s
* sa_loading	12.40s
* bns_loading	65.59s
* matching_time	496.15s
* reads	1000000
* kmers	70000000
* rpm	120932
* kpm	8465241


### 02_bacteria_orig_k31/4.1_read_assignment.log
* Sat Jul 23 02:49:31 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a

```
User time:     0:02:05
System time:   0:00:00
Elapsed time:  0:02:53
CPU usage:          72%

Memory peak:      0.07 GB
```

* 51752inputs+0outputs (198major+88646minor)pagefaults 0swaps


### 02_bacteria_orig_k31/4.2_read_assignment_simlca.log
* Sat Jul 23 02:49:31 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/bacteria.orig.newick -k 31 -f sam -a -t

```
User time:     0:05:29
System time:   0:00:00
Elapsed time:  0:06:17
CPU usage:          87%

Memory peak:      0.07 GB
```

* 44432inputs+0outputs (201major+111886minor)pagefaults 0swaps


### 02_bacteria_orig_k31/5.1_contigs_stats.log
* Number of contigs: 95208424
* Total length: 10124583056
* Average length: 106.34125249253154

* Median length: 34.0
* Number of k-mer occurencies: 7268330336


### 02_bacteria_orig_k31/5.2_index_size.log
* 4.8G	index.fa.31.bit.klcp
* 3.9G	index.fa.ann
* 9.5G	index.fa.bwt
* 4.8G	index.fa.pac
* 4.8G	index.fa.sa

***
## 03_hmp_orig_k31

### 03_hmp_orig_k31/1.1_kmer_propagation.log
* Sat Jul 23 03:15:37 CEST 2016
* jobs: 24
* make -f Makefile.generated

```
User time:     2:53:08
System time:   0:04:24
Elapsed time:  0:12:04
CPU usage:        1470%

Memory peak:      5.81 GB
```

* 11552456inputs+26263224outputs (9major+147070940minor)pagefaults 0swaps


### 03_hmp_orig_k31/1.2_merging_fasta.log
* Sat Jul 23 03:27:42 CEST 2016
* jobs: 24
* ../../bin/create_final_fasta.py index

```
User time:     0:05:27
System time:   0:00:12
Elapsed time:  0:05:57
CPU usage:          95%

Memory peak:      0.01 GB
```

* 56inputs+14285504outputs (2major+25802minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.1_bwa_fa2pac.log
* Sat Jul 23 03:33:39 CEST 2016
* jobs: 24
* ../../bin/bwa fa2pac index.fa index.fa

```
User time:     0:01:07
System time:   0:00:11
Elapsed time:  0:02:18
CPU usage:          57%

Memory peak:      7.06 GB
```

* 1104inputs+9814544outputs (3major+2687492minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.2_bwa_pac2bwtgen.log
* Sat Jul 23 03:35:57 CEST 2016
* jobs: 24
* ../../bin/bwa pac2bwtgen -b 50000000 index.fa.pac index.fa.bwt

```
User time:     1:31:47
System time:   0:00:03
Elapsed time:  1:31:51
CPU usage:         100%

Memory peak:      3.56 GB
```

* 0inputs+6201448outputs (0major+223255minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.3_bwa_bwtupdate.log
* Sat Jul 23 05:07:49 CEST 2016
* jobs: 24
* ../../bin/bwa bwtupdate index.fa.bwt

```
User time:     0:00:25
System time:   0:00:17
Elapsed time:  0:01:52
CPU usage:          38%

Memory peak:      8.87 GB
```

* 272inputs+12402896outputs (0major+2441202minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.4_bwa_bwt2sa.log
* Sat Jul 23 05:09:42 CEST 2016
* jobs: 24
* ../../bin/bwa bwt2sa index.fa.bwt index.fa.sa

```
User time:     0:32:54
System time:   0:01:01
Elapsed time:  0:34:30
CPU usage:          98%

Memory peak:      8.87 GB
```

* 168inputs+6201448outputs (0major+68580216minor)pagefaults 0swaps


### 03_hmp_orig_k31/2.5_klcp.log
* Sat Jul 23 05:44:12 CEST 2016
* jobs: 24
* ../../bin/exk index -k 31 index.fa

```
User time:     1:10:18
System time:   0:01:17
Elapsed time:  1:12:10
CPU usage:          99%

Memory peak:      8.87 GB
```

* 304inputs+6201448outputs (1major+106897669minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1a_matching_rolling.log
* Sat Jul 23 06:56:23 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.1b_matching_rolling.log -k 31 -u index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:45
System time:   0:00:55
Elapsed time:  0:04:14
CPU usage:          87%

Memory peak:     15.02 GB
```

* 15795784inputs+258864outputs (0major+32974271minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.1b_matching_rolling.log
* bwt_loading	13.22s
* sa_loading	9.58s
* bns_loading	29.27s
* klcp_loading	12.20s
* matching_time	155.29s
* reads	1000000
* kmers	70000000
* rpm	386368
* kpm	27045780


### 03_hmp_orig_k31/3.2a_matching_restarted.log
* Sat Jul 23 07:00:38 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.2b_matching_restarted.log -k 31 index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:05:19
System time:   0:00:27
Elapsed time:  0:05:47
CPU usage:         100%

Memory peak:     12.06 GB
```

* 248inputs+258864outputs (3major+11850864minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.2b_matching_restarted.log
* bwt_loading	5.63s
* sa_loading	7.82s
* bns_loading	30.20s
* matching_time	302.49s
* reads	1000000
* kmers	70000000
* rpm	198353
* kpm	13884694


### 03_hmp_orig_k31/3.3a_matching_rolling_skipping.log
* Sat Jul 23 07:06:25 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.3b_matching_rolling_skipping.log -k 31 -u -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:02:36
System time:   0:00:43
Elapsed time:  0:03:20
CPU usage:         100%

Memory peak:     15.01 GB
```

* 0inputs+258864outputs (0major+24915295minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.3b_matching_rolling_skipping.log
* bwt_loading	10.93s
* sa_loading	7.21s
* bns_loading	28.95s
* klcp_loading	7.24s
* matching_time	144.63s
* reads	1000000
* kmers	70000000
* rpm	414861
* kpm	29040236


### 03_hmp_orig_k31/3.4a_matching_restarted_skipping.log
* Sat Jul 23 07:09:46 CEST 2016
* jobs: 24
* ../../bin/exk match -l 3.4b_matching_restarted_skipping.log -k 31 -s index.fa ../../reads/simulation_bacteria.1000000.fq

```
User time:     0:04:51
System time:   0:00:07
Elapsed time:  0:04:59
CPU usage:         100%

Memory peak:     12.06 GB
```

* 136inputs+258864outputs (1major+676194minor)pagefaults 0swaps


### 03_hmp_orig_k31/3.4b_matching_restarted_skipping.log
* bwt_loading	3.53s
* sa_loading	1.70s
* bns_loading	26.14s
* matching_time	267.74s
* reads	1000000
* kmers	70000000
* rpm	224094
* kpm	15686609


### 03_hmp_orig_k31/4.1_read_assignment.log
* Sat Jul 23 07:00:38 CEST 2016
* jobs: 24
* ../../bin/assignment.py -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a

```
User time:     0:01:15
System time:   0:00:00
Elapsed time:  0:01:18
CPU usage:          96%

Memory peak:      0.09 GB
```

* 45872inputs+0outputs (196major+168331minor)pagefaults 0swaps


### 03_hmp_orig_k31/4.2_read_assignment_simlca.log
* Sat Jul 23 07:00:38 CEST 2016
* jobs: 24
* ../../bin/assignment.py -l -i kmers_rolling.txt -n ../../trees/hmp.orig.newick -k 31 -f sam -a -t

```
User time:     0:02:44
System time:   0:00:00
Elapsed time:  0:02:46
CPU usage:          99%

Memory peak:      0.09 GB
```

* 46496inputs+0outputs (200major+280736minor)pagefaults 0swaps


### 03_hmp_orig_k31/5.1_contigs_stats.log
* Number of contigs: 42285628
* Total length: 6350278678
* Average length: 150.17581571686722

* Median length: 34.0
* Number of k-mer occurencies: 5081709838


### 03_hmp_orig_k31/5.2_index_size.log
* 3.0G	index.fa.31.bit.klcp
* 1.8G	index.fa.ann
* 6.0G	index.fa.bwt
* 3.0G	index.fa.pac
* 3.0G	index.fa.sa

