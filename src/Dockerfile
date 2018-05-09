# FROM ubuntu:16.04

# RUN apt-get update \
#     && apt-get install -y software-properties-common curl \
#     && add-apt-repository ppa:jonathonf/python-3.6 \
#     && apt-get remove -y software-properties-common \
#     && apt autoremove -y \
#     && apt-get update \
#     && apt-get install -y python3.6 \
#     && apt-get install -y python3.6-dev \
#     && curl -o /tmp/get-pip.py "https://bootstrap.pypa.io/get-pip.py" \
#     && python3.6 /tmp/get-pip.py \
#     && apt-get remove -y curl \
#     && apt autoremove -y \
#     && rm -rf /var/lib/apt/lists/*

FROM ubuntu:latest

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get update && apt-get -y install apache2 build-essential cpanminus curl git libmysqlclient-dev libpng-dev libssl-dev manpages mysql-client openssl perl perl-base unzip vim wget sudo 
# install ensembl dependencies
RUN cpanm Test::Object PPI::Document Task::Weaken Test::SubCalls Test::Object DBI DBD::mysql Archive::Zip Perl::Critic

RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo apt-get update \
   && sudo apt-get -y install software-properties-common
RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'

USER root
WORKDIR /

RUN apt-get install apt-transport-https

RUN apt-get update && apt-get install -y --no-install-recommends \
 		littler \
 		r-base \
 		r-base-dev \
 		r-recommended


ENV PACKAGE_BIO="samtools libhts1 bedtools"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf zlib1g-dev liblzma-dev libncurses5-dev libblas-dev liblapack-dev libcurl4-gnutls-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev"
#numpy cython scipy transvar bx-python pyvcf cyvcf cyvcf2 biopython crossmap pandas"
ENV PYTHON_MODULES="numpy cython scipy pandas cyvcf2 toml"


RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		apt-utils nano ca-certificates ed less locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*


## Install vcfanno version 0.2.8
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.2.9/vcfanno_linux64 && \
    mv vcfanno_linux64 vcfanno && \
    mv vcfanno /usr/local/bin && \
    chmod 755 /usr/local/bin/vcfanno

## Install Ensembl's Vcf-validator
RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.6/vcf_validator && \
mv vcf_validator /usr/local/bin && \
chmod 755 /usr/local/bin/vcf_validator

## Install pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.1.3/pandoc-2.1.3-1-amd64.deb && \
#RUN wget https://github.com/jgm/pandoc/releases/download/1.19.1/pandoc-1.19.1-1-amd64.deb && \
      dpkg -i pandoc* && \
      rm pandoc* && \
      apt-get clean

USER root
WORKDIR /


## Install tools used for compilation
RUN sudo -H pip install --upgrade pip
RUN sudo -H pip install -U setuptools
RUN sudo -H pip install $PYTHON_MODULES
#
RUN wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz \
		&& gzip -dc nlopt-2.4.2.tar.gz | tar xvf - \
		&& cd nlopt-2.4.2 \
    && ./configure \
		&& make \
		&& make install


RUN apt-get update \
  && apt-get install -y --no-install-recommends libpq-dev libxt-dev libudunits2-dev


###### INSTALL VARIANT EFFECT PREDICTOR (VEP)
# create vep user
RUN useradd -r -m -U -d /home/vep -s /bin/bash -c "VEP User" -p '' vep
RUN usermod -a -G sudo vep
USER vep
ENV HOME /home/vep
WORKDIR $HOME

# clone git repositories
RUN mkdir -p src
WORKDIR $HOME/src
RUN git clone https://github.com/Ensembl/ensembl.git
RUN git clone https://github.com/Ensembl/ensembl-vep.git

# get VEP dependencies
WORKDIR $HOME/src
RUN ensembl-vep/travisci/get_dependencies.sh
ENV PERL5LIB $PERL5LIB:$HOME/src/bioperl-live-release-1-6-924
ENV KENT_SRC $HOME/src/kent-335_base/src
ENV HTSLIB_DIR $HOME/src/htslib
ENV MACHTYPE x86_64
ENV CFLAGS "-fPIC"
ENV DEPS $HOME/src

# and run the complilation/install as root
USER root
RUN ensembl-vep/travisci/build_c.sh

# install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# install bioperl-ext, faster alignments for haplo
WORKDIR $HOME/src
RUN git clone https://github.com/bioperl/bioperl-ext.git
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
RUN perl Makefile.PL
RUN make
RUN make install

# install perl dependencies
WORKDIR $HOME/src
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile .

# switch back to vep user
USER vep

# update bash profile
RUN echo >> $HOME/.profile && \
echo PATH=$HOME/src/ensembl-vep:\$PATH >> $HOME/.profile && \
echo export PATH >> $HOME/.profile

# setup environment
ENV PATH $HOME/src/ensembl-vep:$PATH

# run INSTALL.pl
WORKDIR $HOME/src/ensembl-vep
RUN chmod u+x *.pl
RUN ./INSTALL.pl -a ap -g miRNA -l

##### END INSTALL VEP

USER root
WORKDIR /

RUN R -e "install.packages(c('dplyr','stringr','tidyr','httr','git2r','data.table','magrittr','devtools','htmltools','caret','randomForest','plotly','RcppTOML'), dependencies = T, repos = 'http://cran.us.r-project.org')"
RUN R -e "source(\"https://bioconductor.org/biocLite.R\"); biocLite(c('deconstructSigs', 'VariantAnnotation','BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.UCSC.hg38','GenomeInfoDb','GenomicRanges','S4Vectors'))"
RUN R -e "library(devtools); devtools::install_github('rstudio/DT'); devtools::install_github('hadley/ggplot2'); devtools::install_github('mjkallen/rlogging'); devtools::install_github('kent37/summarywidget')"
#
RUN git clone https://github.com/atks/vt.git
WORKDIR vt
RUN make
RUN make test
RUN cp vt /usr/local/bin
#RUN ln -s /usr/bin/python3.6 /usr/local/bin/python
RUN export PATH=/usr/local/bin:$PATH

## Add local PCGR R package
WORKDIR /
ADD R/ /
RUN R -e "install.packages('configr', dependencies = T, repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_github('rstudio/crosstalk')"
RUN R -e "devtools::install('pcgrr')"

# Add local PCGR Python scripts/libraries
ADD pcgr.tgz /
ENV PATH=$PATH:/pcgr
ENV PYTHONPATH=:/pcgr/lib:${PYTHONPATH}
#ENV VCFANNO_DATA_DOCKER="/data"

WORKDIR /
RUN curl -L -o mskcc-vcf2maf.tar.gz https://api.github.com/repos/mskcc/vcf2maf/tarball/v1.6.16
RUN tar -zxf mskcc-vcf2maf.tar.gz
RUN ln -s /mskcc-vcf2maf-*/vcf2maf.pl /usr/local/bin/vcf2maf.pl

## Clean Up
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
ADD pcgr_predispose.R /

USER root
WORKDIR /
RUN rm -f nlopt-2.4.2.tar.gz
RUN rm -rf $HOME/src/ensembl-vep/t/
RUN rm -f $HOME/src/v335_base.tar.gz
RUN rm -f $HOME/src/release-1-6-924.zip
RUN rm -rf /vt
