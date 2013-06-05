strain=$1

curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/$strain/ -l -s | grep 'fna'
