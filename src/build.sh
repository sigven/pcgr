cp /Users/sigven/research/software/vcf2tsv/vcf2tsv.py pcgr/
tar czvfhL pcgr.tgz pcgr/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/pcgr:$TAG --rm=true .
