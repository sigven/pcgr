#cp -LR ../../../../ncgc/R/pcgrr R/
#rm -rf R/pcgrr/.git*
#rm -rf R/pcgrr/.R*
#rm -rf R/pcgrr/rlogging.log
tar czvfh pcgr.tgz pcgr/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/pcgr:$TAG --rm=true .

