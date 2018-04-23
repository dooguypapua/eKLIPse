# Base build, inspired by
# https://github.com/crosbymichael/python-docker/blob/master/Dockerfile
FROM dennishazelett/circos
MAINTAINER David Goudeneige, dooguy@tuta.io

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
    ncbi-blast+ \
    && apt-get autoremove \
    && apt-get clean

RUN pip install numpy
RUN pip install biopython
RUN pip install tqdm
RUN pip install blast

CMD []
