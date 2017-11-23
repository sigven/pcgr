tar czvfh pcgr.tgz pcgr/
echo "Build the Docker Image"
TAG=`date "+%Y%m%d"`
docker build -t sigven/pcgr:$TAG --rm=true .

