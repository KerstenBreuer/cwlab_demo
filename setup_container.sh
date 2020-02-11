apt-get -y install --fix-missing \
    wget \
    tar \
    unzip \
    ca-certificates \
    openjdk-8-jdk \
    pigz \
    unzip \
    bzip2 \
	cmake \
	zlib1g-dev  \
	bzip2 \
	libghc-bzlib-dev \
	liblzma-dev  \
	libncurses5-dev \
    libhts-dev \
    libtbb-dev \
    autoconf \
    nodejs

pip3 install cutadapt==2.6

cd /usr/local \
	&& wget -O fastqc.zip http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip \
	&& unzip fastqc.zip \
	&& rm fastqc.zip \
	&& cd FastQC \
	&& chmod a+x fastqc
echo 'export PATH="/usr/local/FastQC:${PATH}"' >> ~/.bashrc

cd /usr/local \
  && wget -O trim_galore.tar.gz https://github.com/FelixKrueger/TrimGalore/archive/0.6.4.tar.gz \
  && tar -xvzf trim_galore.tar.gz \
  && mv TrimGalore*/trim_galore /usr/local/bin \
  && rm -r TrimGalore* trim_galore.tar.gz

cd /usr/local/ \
    && wget -O bowtie2.zip \
    https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip \
	&& unzip bowtie2.zip \
	&& rm bowtie2.zip \
	&& mv bowtie2* bowtie2
echo 'export PATH="/usr/local/bowtie2:${PATH}"' >> ~/.bashrc

cd /usr/local/ \
    && wget -O samtools.tar.bz2 \
    https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
	&& tar -jxf samtools.tar.bz2 \
	&& rm samtools.tar.bz2 \
	&& mv samtools* samtools \
	&& cd samtools \
    && autoheader \
    && autoconf -Wno-syntax \
    && ./configure \
	&& make \
	&& make install
echo 'export PATH="/usr/local/samtools:${PATH}"' >> ~/.bashrc

source ~/.bashrc