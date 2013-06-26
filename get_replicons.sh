strain=$1

curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/$strain/ -l -s | grep 'fna'
curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria_DRAFT/$strain/ -l -s | grep 'fna'
