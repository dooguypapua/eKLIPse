# Base build, inspired by
# https://github.com/crosbymichael/python-docker/blob/master/Dockerfile
FROM dennishazelett/circos
MAINTAINER David Goudenege, dooguy@tuta.io

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --fix-missing \
    build-essential \
    ca-certificates \
    gcc \
    git \
    libpq-dev \
    make \
    python-pip \
    python2.7 \
    python2.7-dev \
    samtools \
    unzip \
    && apt-get autoremove \
    && apt-get clean
# apt install ncbi-blast+ is the version 2.2.28, we need >2.3.0 to work with eklipse

RUN pip install numpy
RUN pip install biopython
RUN pip install tqdm

# retrieving eKLIPse
RUN git clone https://github.com/dooguypapua/eKLIPse.git
RUN cd /eKLIPse

# Getting good version of blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz
RUN tar -xvzf ncbi-blast-2.7.1+-x64-linux.tar.gz

ENV PATH="/eKLIPse/ncbi-blast-2.7.1+/bin/:${PATH}"

#CMD ["python eKLIPse.py --test"]
