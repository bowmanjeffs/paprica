FROM python:3.9

LABEL maintainer="Jeff Bowman"

# Install packages
RUN apt-get update && \
    apt-get install -qy --no-install-recommends \
    make \
    git \
    cmake \
    autotools-dev \
    libtool \
    flex \
    bison \
    cmake \
    automake \
    autoconf \
    build-essential \ 
    git \
    zip

# Install python dependencies, including external python tools
RUN pip3 install numpy biopython joblib pandas seqmagick termcolor
    
RUN cd /

# Install RAxML-ng
RUN wget --no-verbose https://github.com/amkozlov/raxml-ng/releases/download/0.9.0/raxml-ng_v0.9.0_linux_x86_64.zip && \
    unzip raxml-ng_v0.9.0_linux_x86_64.zip && \
    rm raxml-ng_v0.9.0_linux_x86_64.zip

# Install infernal
RUN wget --no-verbose http://eddylab.org/infernal/infernal-1.1.2-linux-intel-gcc.tar.gz && \
    tar -xzvf infernal-1.1.2-linux-intel-gcc.tar.gz && \
    mv infernal-1.1.2-linux-intel-gcc infernal && \
    rm infernal-1.1.2-linux-intel-gcc.tar.gz

# Install gappa
RUN git clone --recursive https://github.com/lczech/gappa.git && \
    make -C /gappa

# Install epa-ng
RUN git clone https://github.com/Pbdas/epa-ng.git && \ 
    make -C /epa-ng

# Modify PATH
ENV PATH="/pplacer:${PATH}"
ENV PATH="/.local/bin:${PATH}"
ENV PATH="/infernal/binaries:${PATH}"
ENV PATH="/infernal/easel:${PATH}"
ENV PATH="/raxml-ng:${PATH}"
ENV PATH="/paprica:${PATH}"
ENV PATH="/epa-ng/bin:${PATH}"
ENV PATH="/gappa/bin:${PATH}"

# Install paprica
RUN git clone https://github.com/bowmanjeffs/paprica.git && cd paprica && chmod a+x *py && chmod a+x *sh

# Run bash on container startup
CMD "/bin/bash"
