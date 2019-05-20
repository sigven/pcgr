

###################################################
# Stage 1 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:18.04 as builder

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    git \
    libpng-dev \
    perl \
    perl-base \
    unzip \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV HTSLIB_DIR $OPT_SRC/htslib
# Branch to clone, e.g. "-b release/96"
ENV BRANCH ""

# Working directory
WORKDIR $OPT_SRC

# Clone/download repositories/libraries
# Clone ensembl git repository and extract useful ensemb core file
RUN git clone $BRANCH --depth 1 https://github.com/Ensembl/ensembl.git && \
    cp ensembl/cpanfile ensembl_cpanfile && \
    rm -rf ensembl && \
    # Clone ensembl-vep git repository
    git clone $BRANCH --depth 1 https://github.com/Ensembl/ensembl-vep.git && chmod u+x ensembl-vep/*.pl && \
    # Clone ensembl-variation git repository and compile C code
    git clone $BRANCH --depth 1 https://github.com/Ensembl/ensembl-variation.git && \
    mkdir var_c_code && \
    cp ensembl-variation/C_code/*.c ensembl-variation/C_code/Makefile var_c_code/ && \
    rm -rf ensembl-variation && \
    chmod u+x var_c_code/* && \
    # Clone bioperl-ext git repository
    git clone --depth 1 https://github.com/bioperl/bioperl-ext.git && \
    # Download ensembl-xs
    wget https://github.com/Ensembl/ensembl-xs/archive/2.3.2.zip -O ensembl-xs.zip && \
    unzip -q ensembl-xs.zip && mv ensembl-xs-2.3.2 ensembl-xs && rm -rf ensembl-xs.zip && \
    # Clone/Download other repositories: bioperl-live is needed so the cpanm dependencies installation from the ensembl-vep/cpanfile file takes less disk space
    ensembl-vep/travisci/get_dependencies.sh && \
    # Only keep the bioperl-live "Bio" library
    mv bioperl-live bioperl-live_bak && mkdir bioperl-live && mv bioperl-live_bak/Bio bioperl-live/ && rm -rf bioperl-live_bak && \
    ## A lot of cleanup on the imported libraries, in order to reduce the docker image ##
    rm -rf Bio-HTS/.??* Bio-HTS/Changes Bio-HTS/DISCLAIMER Bio-HTS/MANIFEST* Bio-HTS/README Bio-HTS/scripts Bio-HTS/t Bio-HTS/travisci \
           bioperl-ext/.??* bioperl-ext/Bio/SeqIO bioperl-ext/Bio/Tools bioperl-ext/Makefile.PL bioperl-ext/README* bioperl-ext/t bioperl-ext/examples \
           ensembl-vep/.??* ensembl-vep/docker \
           ensembl-xs/.??* ensembl-xs/Changes ensembl-xs/INSTALL ensembl-xs/MANIFEST ensembl-xs/README ensembl-xs/t ensembl-xs/travisci \
           htslib/.??* htslib/INSTALL htslib/NEWS htslib/README* htslib/test && \
    # Only keep needed kent-335_base libraries for VEP
    mv kent-335_base kent-335_base_bak && mkdir -p kent-335_base/src && \
    cp -R kent-335_base_bak/confs kent-335_base/ && \
    cp -R kent-335_base_bak/src/lib kent-335_base_bak/src/inc kent-335_base_bak/src/jkOwnLib kent-335_base/src/ && \
    cp kent-335_base_bak/src/*.sh kent-335_base/src/ && \
    rm -rf kent-335_base_bak

# Setup bioperl-ext
WORKDIR bioperl-ext/Bio/Ext/Align/
RUN perl -pi -e"s|(cd libs.+)CFLAGS=\\\'|\$1CFLAGS=\\\'-fPIC |" Makefile.PL

# Install htslib binaries (need bgzip, tabix)
WORKDIR $HTSLIB_DIR
RUN make install && rm -f Makefile *.c cram/*.c

# Compile Variation LD C scripts
WORKDIR $OPT_SRC/var_c_code
RUN make && rm -f Makefile *.c


###################################################
# Stage 2 - docker container to build ensembl-vep #
###################################################
FROM ubuntu:18.04

# Update aptitude and install some required packages
# a lot of them are required for Bio::DB::BigFile
RUN apt-get update && apt-get -y install \
    build-essential \
    cpanminus \
    curl \
    libmysqlclient-dev \
    libpng-dev \
    libssl-dev \
    locales \
    openssl \
    perl \
    perl-base \
    unzip \
    vim && \
    apt-get -y purge manpages-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Setup VEP environment
ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PERL5LIB_TMP $PERL5LIB:$OPT_SRC/ensembl-vep:$OPT_SRC/ensembl-vep/modules
ENV PERL5LIB $PERL5LIB_TMP:$OPT_SRC/bioperl-live
ENV KENT_SRC $OPT/src/kent-335_base/src
ENV HTSLIB_DIR $OPT_SRC/htslib
ENV MACHTYPE x86_64
ENV DEPS $OPT_SRC
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
ENV LANG_VAR en_US.UTF-8

# Create vep user
RUN useradd -r -m -U -d "$OPT" -s /bin/bash -c "VEP User" -p '' vep && usermod -a -G sudo vep && mkdir -p $OPT_SRC
USER vep

# Copy downloaded libraries (stage 1) to this image (stage 2)
COPY --chown=vep:vep --from=builder $OPT_SRC $OPT_SRC
#############################################################

# Change user to root for the following complilations/installations
USER root

# Install bioperl-ext, faster alignments for haplo (XS-based BioPerl extensions to C libraries)
WORKDIR $OPT_SRC/bioperl-ext/Bio/Ext/Align/
RUN perl Makefile.PL && make && make install && rm -f Makefile*

# Install ensembl-xs, faster run using re-implementation in C of some of the Perl subroutines
WORKDIR $OPT_SRC/ensembl-xs
RUN perl Makefile.PL && make && make install && rm -f Makefile*

WORKDIR $OPT_SRC
# Install/compile more libraries
RUN ensembl-vep/travisci/build_c.sh && \
    # Install ensembl perl dependencies (cpanm)
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl_cpanfile . && \
    cpanm --installdeps --with-recommends --notest --cpanfile ensembl-vep/cpanfile . && \
    # Delete bioperl after the cpanm installs as it will be reinstalled by the INSTALL.pl script
    rm -rf bioperl-live && \
    # Configure "locale", see https://github.com/rocker-org/rocker/issues/19
    echo "$LANG_VAR UTF-8" >> /etc/locale.gen && locale-gen en_US.utf8 && \
    /usr/sbin/update-locale LANG=$LANG_VAR && \
    # Copy htslib executables
    cp $HTSLIB_DIR/bgzip $HTSLIB_DIR/tabix $HTSLIB_DIR/htsfile /usr/local/bin/

ENV LC_ALL $LANG_VAR
ENV LANG $LANG_VAR

# Switch back to vep user
USER vep
ENV PERL5LIB $PERL5LIB_TMP

# Final steps
WORKDIR $OPT_SRC/ensembl-vep
# Update bash profile
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile && \
    # Run INSTALL.pl and remove the ensemb-vep tests and travis
    ./INSTALL.pl -a ap -g miRNA,LoF -l && rm -rf t travisci .travis.yml

WORKDIR /
ADD loftee.tgz $OPT/src/ensembl-vep/modules


####################################################
# Stage 3 - adding functionality for PCGR analysis #
####################################################

USER root

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

ENV TZ=Europe/Minsk
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
#RUN sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update && apt-get -y install r-base

USER root
WORKDIR /

ENV PACKAGE_BIO="libhts2 bedtools"
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

RUN apt-get update \
  && apt-get install -y --no-install-recommends libpq-dev libxt-dev libudunits2-dev

RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
RUN bunzip2 -dc samtools-1.9.tar.bz2 | tar xvf -
RUN cd samtools-1.9 && ./configure --prefix=/usr/local/bin && make && make install


WORKDIR /

RUN R -e "install.packages(c('BiocManager','devtools'))"
RUN R -e "BiocManager::install(\"ComplexHeatmap\")"
RUN R -e "BiocManager::install(\"VariantAnnotation\")"
RUN R -e "BiocManager::install(\"Biostrings\")"

## Install maftools from github repository.
RUN R -e "library(\"devtools\")"
RUN R -e "devtools::install_github(repo = \"PoisonAlien/maftools\")"

## Install other bioconductor/cran packages
RUN R -e "BiocManager::install(c('regioneR','MutationalPatterns','deconstructSigs', 'BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.UCSC.hg38','GenomeInfoDb','GenomicRanges','S4Vectors','karyoploteR'))"
RUN R -e "install.packages(c('configr','rmarkdown','httr','git2r','data.table','tidyverse','htmltools','caret','randomForest','plotly','RcppTOML','SeqKat'), dependencies = T, repos = 'http://cran.us.r-project.org')"
RUN R -e "library(devtools); devtools::install_github('rstudio/DT'); devtools::install_github('mjkallen/rlogging'); devtools::install_github('kent37/summarywidget')"
RUN R -e "library(devtools); devtools::install_github('rstudio/crosstalk')"
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
RUN wget https://github.com/jgm/pandoc/releases/download/2.6/pandoc-2.6-1-amd64.deb && \
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
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake

RUN wget https://github.com/stevengj/nlopt/archive/v2.5.0.tar.gz \
		&& gzip -dc v2.5.0.tar.gz | tar xvf - \
		&& cd nlopt-2.5.0 \
    && cmake . \
		&& make \
		&& make install

USER root
WORKDIR /


RUN git clone https://github.com/atks/vt.git
WORKDIR vt
RUN make
RUN make test
RUN cp vt /usr/local/bin
RUN export PATH=/usr/local/bin:$PATH

## Install mosdepth through miniconda
WORKDIR /
ENV HOME=/usr/local
RUN curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
ENV PATH=$PATH:/usr/local/miniconda3/bin:
RUN conda config --add channels bioconda
RUN conda install mosdepth
RUN rm -rf /Miniconda3-latest-Linux-x86_64.sh

## Add local PCGR R package
WORKDIR /
ADD R/ /
RUN R -e "devtools::install_github('mjkallen/rlogging')"
RUN R -e "devtools::install('pcgrr')"
RUN R -e "devtools::install_github('hms-dbmi/UpSetR')"
RUN R -e "devtools::install_github('kassambara/ggpubr')"

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

RUN rm -f nlopt-2.5.0.tar.gz
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
