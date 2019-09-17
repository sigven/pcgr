echo "Build the Docker Image"
#git clone --branch grch38 --single-branch https://github.com/konradjk/loftee.git
#cd loftee
#tar czf loftee_1.0.3.tgz *
#cd ..
#mv loftee/loftee_1.0.3.tgz .
#rm -rf loftee/
cp /Users/sigven/research/software/vcf2tsv/vcf2tsv.py pcgr/
tar czvfhL pcgr.tgz pcgr/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/pcgr:$TAG --rm=true .
