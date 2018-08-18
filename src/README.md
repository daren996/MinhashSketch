### Compile

$ nvcc -std=c++11 -I /export/project/hondius/opt/cub -rdc=true main.cu Utils.cu Hash.cu SpookyV2_d.cu -o minhash

### Run

$ ./minhash ../testing_files/sequence_clip1.fasta ../testing_files/sequence_clip2.fasta all -e --k=5 --m=10 --t=10

$ ./minhash ../testing_files/sequence_abstract1.fasta ../testing_files/sequence_abstract2.fasta all -e --k=5 --m=10 --t=10
