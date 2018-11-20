samtools view ../../projects/cnag_rnaseq_buffer_optimization/data/alignments/c_elegans/CCB68ANXX_4_N703-S503-NX-xt.sorted.bam |\
 grep AAAAAAAAAAA | head -n 20 | awk '{print $1}' > ./test/read_names.txt

samtools view ../../projects/cnag_rnaseq_buffer_optimization/data/alignments/c_elegans/CCB68ANXX_4_N703-S503-NX-xt.sorted.bam |\
 grep -f ./test/read_names.txt | head -n 40 > ./test/paired_ends.sam

awk 'NR%2==0{print ">"$1"\n"$10}' ./test/paired_ends.sam > ./test/paired_ends_1.fa
awk 'NR%2==1{print ">"$1"\n"$10}' ./test/paired_ends.sam > ./test/paired_ends_2.fa

awk 'NR%2==0{print "@"$1"\n"$10"\n+\n"$11}' ./test/paired_ends.sam > ./test/paired_ends_1.fq
awk 'NR%2==1{print "@"$1"\n"$10"\n+\n"$11}' ./test/paired_ends.sam > ./test/paired_ends_2.fq

awk 'NR%2==0{print $0}' ./test/single_ends.sam > ./test/single_ends.sam
awk '{print ">"$1"\n"$10}' ./test/single_ends.sam > ./test/single_ends.fa
awk '{print "@"$1"\n"$10"\n+\n"$11}' ./test/single_ends.sam > ./test/single_ends_1.fq
