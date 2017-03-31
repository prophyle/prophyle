To check whether kmers mathcing works correctly,
just run 
./exk -u -v -k 2 test/ref.fa test/reads.fq > test/current.txt
from upper directory
and check that test/current.txt is the same as test/expected.txt
