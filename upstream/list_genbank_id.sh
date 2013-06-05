if [ "$#" -lt 2 ]; then
	echo 'list_genbank_id GENUS SPECIES'
	exit 65
fi

genus=$1
species=$2

curl ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ -l -s | grep  $genus"_"$species

