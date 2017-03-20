Getting started
---------------

Prerequisites
~~~~~~~~~~~~~

Docker
^^^^^^

-  Running PCGR requires that Docker is set up on your host. Docker has
   very complete installation instructions for different platforms:

   -  installing `Docker on
      Linux <https://docs.docker.com/engine/installation/linux/>`__
   -  installing `Docker on Mac
      OS <https://docs.docker.com/engine/installation/mac/>`__

-  TODO: Check that Docker is running

-  Adjust the computing resources dedicated to the Docker, i.e.:

   -  Memory: minimum 5GB
   -  CPUs: minimum 4

Installation of PCGR
^^^^^^^^^^^^^^^^^^^^

Below follows step-by-step instructions for installation, using standard
Unix/command line utilities:

-  Download and unpack the `latest
   release <https://github.com/sigven/pcgr/releases/tag/v0.2>`__

-  Download and unpack the data bundle (approx. 16Gb) in the PCGR
   directory

   -  Download `data
      bundle <https://drive.google.com/open?id=0B8aYD2TJ472mUFVXcmo1ZXY0OWM>`__
      from Google Drive to ``~/pcgr-X.X`` (replace *X.X* with the
      version number)
   -  Decompress and untar the bundle, e.g.
      ``gzip -dc pcgr.bundle.latest.grch37.tgz | tar xvf -``

-  Pull the PCGR Docker image from DockerHub:

   -  ``docker pull sigven/pcgr:latest`` (PCGR annotation engine)

Test PCGR package - generation of clinical report for a cancer genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generate report for the example tumor genome (SNVs/InDels only) present
in the *examples/* folder:

``python run_pcgr.py --input_vcf test.vcf ~/pcgr-X.X ~/pcgr-X.X/examples test_sample``
