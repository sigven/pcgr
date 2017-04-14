############################################################
# Dockerfile to build Personal Cancer Genome Reporter (PCGR)
# Main software components:
#  1. Variant Effect Predictor (VEP)
#  2. vcfanno (0.0.11)
#  3. custom scripts (pcgr.tgz)
############################################################

# use Debian:testing as base image
FROM debian:jessie-backports

# set non interactive
ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGE_BIO="tabix samtools libhts1 bedtools"
ENV PACKAGE_DEV="perl debconf-utils build-essential gfortran python-dev python-pip gcc-multilib autoconf zlib1g-dev git libncurses5-dev libblas-dev liblapack-dev cpanminus libcurl4-gnutls-dev libssh2-1-dev libxml2-dev vim libssl-dev openssl libcairo2-dev"
ENV PYTHON_MODULES="numpy cython scipy transvar bx-python pyvcf cyvcf cyvcf2 biopython crossmap pandas"

ENV R_BASE_VERSION 3.3.3

USER root
WORKDIR /

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		apt-utils wget curl nano unzip ca-certificates ed less locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
#deb http://ftp.debian.org/debian jessie-backports main
RUN echo "deb http://http.debian.net/debian jessie-backports main" > /etc/apt/sources.list.d/debian-jessie-backports.list \
	&& echo 'APT::Default-Release "jessie-backports";' > /etc/apt/apt.conf.d/default
#RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
	#&& echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update && apt-get install -y --no-install-recommends \
		littler \
    #r-cran-littler \
		r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}* \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# Install pandoc (for HTML reports)
RUN wget https://github.com/jgm/pandoc/releases/download/1.19.1/pandoc-1.19.1-1-amd64.deb && \
      dpkg -i pandoc* && \
      rm pandoc* && \
      apt-get clean

# Install necessary R packages
RUN R -e "install.packages(c('dplyr','stringr','tidyr','ggplot2','httr','git2r','data.table','magrittr','devtools','DT'), dependencies = T, repos = 'http://cran.us.r-project.org')"
RUN R -e "source(\"https://bioconductor.org/biocLite.R\"); biocLite(c('deconstructSigs', 'KEGGREST','VariantAnnotation','BSgenome.Hsapiens.UCSC.hg19','GenomeInfoDb','GenomicRanges','S4Vectors'))"
RUN R -e "library(devtools); devtools::install_github('mjkallen/rlogging')"


# Install tools used for compilation
RUN pip install -U setuptools
RUN pip install $PYTHON_MODULES

# Install vcfanno version 0.0.11
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.0.11/vcfanno_0.0.11_linux_amd64.tar.gz && \
    tar xvzf vcfanno_0.0.11_linux_amd64.tar.gz && \
    mv vcfanno_0.0.11_linux_amd64/vcfanno /usr/local/bin && \
    rm vcfanno_0.0.11_linux_amd64.tar.gz && \
    rm -rf vcfanno_0.0.11_linux_amd64


# Install Ensembl's Vcf-validator
RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.4.2/vcf_validator && \
mv vcf_validator /usr/local/bin && \
chmod 755 /usr/local/bin/vcf_validator

# Install VEP
# Necessary Perl modules
RUN cpanm File::ShareDir::Install \
 && cpanm Data::UUID \
 && cpanm autodie \
 && cpanm Test::Fatal \
 && cpanm IPC::System::Simple \
 && cpanm File::Fetch \
 && cpanm File::Copy::Recursive \
 && cpanm LWP::Simple \
 && cpanm Bio::Root::Version \
 && cpanm Archive::Extract \
 && cpanm Archive::Zip \
 && cpanm CGI \
 && cpanm DBI \
 && cpanm Archive::Tar
# && cpanm Bio::DB::HTS \

VOLUME        /usr/local/share/vep/data
ENV VEP_DATA="/usr/local/share/vep/data"
ENV VEP_DATA_DOCKER="/usr/local/share/vep/data"
ENV VEP_PATH="/usr/local/share/vep/vep"
ENV VEP_PLUGINS="/usr/local/share/vep/Plugins"
ENV VEP_INSTDIR="/ensembl-tools/scripts/variant_effect_predictor/"
ENV PERL5LIB=$PERL5LIB:$VEP_PATH:$VEP_PLUGINS
ENV PATH=$PATH:$VEP_PATH/htslib
ENV SPECIES="homo_sapiens"
ENV ASSEMBLY="GRCh37"
ENV VEPPLUGIN="LoF,TSSDistance"
RUN git clone -b release/85 https://github.com/Ensembl/ensembl-tools.git
WORKDIR $VEP_INSTDIR
RUN printf 'y\n' | perl INSTALL.pl --AUTO a --SPECIES $SPECIES --ASSEMBLY $ASSEMBLY --PLUGINS $VEPPLUGIN --DESTDIR $VEP_PATH --CACHEDIR $VEP_DATA
RUN cp *.pl /usr/local/bin/
WORKDIR /

## Add vcfbreakmulti (for future preprocessing of input)
RUN git clone --recursive https://github.com/ekg/vcflib.git
WORKDIR vcflib
RUN make
RUN cp bin/vcfbreakmulti /usr/local/bin
RUN rm -rf /vcflib

# Add local PCGR R package
WORKDIR /
ADD R/ /
RUN R -e "devtools::install('pcgrr2')"

# Add local PCGR Python scripts/libraries
ADD pcgr.tgz /
ENV PATH=$PATH:/pcgr
ENV PYTHONPATH=:/pcgr/lib:${PYTHONPATH}
ENV VCFANNO_DATA_DOCKER="/data"

# Clean Up
RUN apt-get clean autoclean
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN rm -rf /var/lib/{dpkg,cache,log}

VOLUME /workdir
WORKDIR /workdir/
USER root
RUN mkdir /data && chmod 777 /data
WORKDIR /data
VOLUME /data
WORKDIR /
ADD pcgr.R /
