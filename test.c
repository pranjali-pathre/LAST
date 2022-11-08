FROM ubuntu:18.04
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get upgrade -y && apt-get install git python-pip -y
RUN apt-get install wget curl -y

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN conda --version

RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

RUN pip install psutil
RUN pip install scikit-learn --no-cache-dir
RUN pip install openmm --no-cache-dir
RUN pip install mdtraj --no-cache-dir
RUN pip install MDAnalysis --no-cache-dir
RUN pip install tensorflow --no-cache-dir
RUN pip install statsmodels --no-cache-dir
RUN pip install openmm
WORKDIR /opt/src/

ENV PATH="/opt/venv/bin:$PATH"
COPY LAST /opt/
RUN chmod +x ./run.sh
ENTRYPOINT ["./run.sh"]