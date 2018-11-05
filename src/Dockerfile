FROM ubuntu:xenial

RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

RUN apt-get update && apt-get -y install apache2 apt-utils build-essential cpanminus curl git libmysqlclient-dev libpng-dev libssl-dev manpages mysql-client openssl perl perl-base unzip vim wget sudo
# install ensembl dependencies
RUN cpanm Test::Object PPI::Document Task::Weaken Test::SubCalls Test::Object DBI DBD::mysql Archive::Zip Perl::Critic Set::IntervalTree

RUN apt-get update && apt-get install apt-transport-https


RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN sudo apt-get update \
   && sudo apt-get -y install software-properties-common
RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts1 bedtools"
ENV PACKAGE_DEV="gfortran gcc-multilib autoconf liblzma-dev libncurses5-dev libblas-dev liblapack-dev libssh2-1-dev libxml2-dev vim libssl-dev libcairo2-dev libbz2-dev libcurl4-openssl-dev"
ENV PYTHON_MODULES="numpy cython scipy pandas cyvcf2 toml"
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		nano ed locales vim-tiny fonts-texgyre \
    $PACKAGE_DEV $PACKAGE_BIO \
	&& rm -rf /var/lib/apt/lists/*

RUN apt-get autoremove

RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base-core \
    r-recommended \
 		r-base

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN bunzip2 -dc samtools-1.9.tar.bz2 | tar xvf -
RUN cd samtools-1.9 && ./configure --prefix=/usr/local/bin && make && make install


WORKDIR /

RUN R -e "source(\"https://bioconductor.org/biocLite.R\"); biocLite(c('regioneR','MutationalPatterns','deconstructSigs', 'VariantAnnotation','BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.UCSC.hg38','GenomeInfoDb','GenomicRanges','S4Vectors','karyoploteR'))"
RUN R -e "install.packages(c('configr','rmarkdown','dplyr','stringr','tidyr','httr','git2r','data.table','magrittr','devtools','htmltools','caret','randomForest','plotly','RcppTOML','SeqKat'), dependencies = T, repos = 'http://cran.us.r-project.org')"
RUN R -e "library(devtools); devtools::install_github('rstudio/DT'); devtools::install_github('hadley/ggplot2'); devtools::install_github('mjkallen/rlogging'); devtools::install_github('kent37/summarywidget')"
RUN R -e "library(devtools); devtools::install_github('Francescojm/CELLector'); devtools::install_github('rstudio/crosstalk')"
RUN rm -rf /tmp/R*


## Install vcfanno version 0.3.1
RUN wget https://github.com/brentp/vcfanno/releases/download/v0.3.1/vcfanno_linux64 && \
    mv vcfanno_linux64 vcfanno && \
    mv vcfanno /usr/local/bin && \
    chmod 755 /usr/local/bin/vcfanno

## Install Ensembl's Vcf-validator
RUN wget https://github.com/EBIvariation/vcf-validator/releases/download/v0.6/vcf_validator && \
mv vcf_validator /usr/local/bin/ && \
chmod 755 /usr/local/bin/vcf_validator

## Install pandoc (for HTML report generation)
RUN wget https://github.com/jgm/pandoc/releases/download/2.3.1/pandoc-2.3.1-1-amd64.deb && \
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

### START VEP Install

RUN apt-get update && apt-get -y install \
      apache2 \
      build-essential \
      cpanminus \
      curl \
      git \
      libmysqlclient-dev \
      libpng-dev \
      libssl-dev \
      locales \
      manpages \
      mysql-client \
      openssl \
      perl \
      perl-base \
      unzip \
      vim \
      wget

# install ensembl dependencies
RUN cpanm DBI DBD::mysql

# create vep user
RUN useradd -r -m -U -d /opt/vep -s /bin/bash -c "VEP User" -p '' vep
RUN usermod -a -G sudo vep
USER vep
ENV OPT /opt/vep
WORKDIR $OPT

# clone git repositories
RUN mkdir -p src
WORKDIR $OPT/src
RUN git clone https://github.com/Ensembl/ensembl.git
RUN git clone https://github.com/Ensembl/ensembl-vep.git
# get VEP dependencies

WORKDIR $OPT/src
RUN ensembl-vep/travisci/get_dependencies.sh
ENV PERL5LIB $PERL5LIB:$OPT/src/bioperl-live-release-1-6-924
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT/src/htslib
ENV MACHTYPE x86_64
ENV CFLAGS "-fPIC"
ENV DEPS $OPT/src

# and run the complilation/install as root
USER root
RUN ensembl-vep/travisci/build_c.sh

# install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install

# install bioperl-ext, faster alignments for haplo
WORKDIR $OPT/src
RUN git clone https://github.com/bioperl/bioperl-ext.git
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL
RUN perl Makefile.PL
RUN make
RUN make install

# install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT/src
RUN wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip
RUN unzip -q ensembl-xs.zip
RUN mv ensembl-xs-2.3.2 ensembl-xs
RUN rm ensembl-xs.zip
WORKDIR ensembl-xs
RUN perl Makefile.PL
RUN make
RUN make install

# install perl dependencies
WORKDIR $OPT/src
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl/cpanfile .
RUN cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile .

# configure locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
locale-gen en_US.utf8 && \
/usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# switch back to vep user
USER vep

# update bash profile
RUN echo >> $OPT/.profile && \
echo PATH=$OPT/src/ensembl-vep:\$PATH >> $OPT/.profile && \
echo export PATH >> $OPT/.profile

# setup environment
ENV PATH $OPT/src/ensembl-vep:$PATH

# run INSTALL.pl
WORKDIR $OPT/src/ensembl-vep
RUN chmod u+x *.pl
RUN ./INSTALL.pl -a ap -g miRNA,LoF -l

WORKDIR /
ADD loftee.tgz $OPT/src/ensembl-vep/modules


##### END INSTALL VEP

USER root
WORKDIR /


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
RUN R -e "devtools::install_github('mjkallen/rlogging')"
RUN R -e "devtools::install('pcgrr')"

# Add local PCGR Python scripts/libraries
ADD pcgr.tgz /
ENV PATH=$PATH:/pcgr
ENV PYTHONPATH=:/pcgr/lib:${PYTHONPATH}
#ENV VCFANNO_DATA_DOCKER="/data"

WORKDIR /
RUN curl -L -o mskcc-vcf2maf.tar.gz https://api.github.com/repos/mskcc/vcf2maf/tarball/v1.6.16
RUN tar -zxf mskcc-vcf2maf.tar.gz
RUN chmod 755 /mskcc-vcf2maf-*/*.pl
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
ADD cpsr.R /

USER root
WORKDIR /
WORKDIR /

RUN rm -f nlopt-2.4.2.tar.gz
RUN rm -rf $HOME/src/ensembl-vep/t/
RUN rm -f $HOME/src/v335_base.tar.gz
RUN rm -f $HOME/src/release-1-6-924.zip
RUN rm -rf /vt
RUN rm -rf /samtools-1.9.tar.bz2
RUN rm -rf /mskcc-vcf2maf.tar.gz

ENV PATH=/usr/local/bin/bin:$PATH

# WORKDIR /
# ENV HOME=/usr/local
# RUN curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# RUN bash Miniconda3-latest-Linux-x86_64.sh -b
# ENV PATH=$PATH:/usr/local/miniconda3/bin:
# RUN conda config --add channels bioconda
# RUN conda install mosdepth
# RUN rm -rf /Miniconda3-latest-Linux-x86_64.sh
