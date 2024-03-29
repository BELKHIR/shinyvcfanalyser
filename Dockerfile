FROM rocker/r-ver:4.2.3
	ENV PATH /opt/biotools/bin:$PATH
	ENV ROOTSYS /opt/biotools/root
	ENV LD_LIBRARY_PATH '$LD_LIBRARY_PATH:$ROOTSYS/lib'
    ENV  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib \
    ENV  export CFLAGS="-I/usr/local/include" \
    ENV  export LDFLAGS="-L/usr/local/lib" \
	RUN apt-get update

	RUN apt-get install -yq tzdata
	RUN apt-get install -y locales
	RUN locale-gen "en_US.UTF-8"
	RUN export LC_ALL=en_US.UTF-8
	RUN export LANG=en_US.UTF-8
	
	RUN apt-get update && apt-get install -y curl wget apt-utils
	RUN apt-get install -y gcc fort77 aptitude
	RUN aptitude install -y g++ xorg-dev libreadline-dev  gfortran
	#RUN apt-get install -y software-properties-common &&  apt-add-repository universe
	#RUN apt-get update
	RUN apt-get install -y libssl-dev libxml2-dev libpcre3-dev liblzma-dev libbz2-dev libcurl4-openssl-dev liblapack3 git nano graphviz libproj-dev
	RUN apt-get update
	RUN apt-get install -y  autotools-dev automake cmake grep sed dpkg fuse zip build-essential pkg-config bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 mercurial subversion zlib1g-dev libncurses5-dev libncursesw5-dev
	RUN apt-get install -y  libgsl-dev  python-tk
	RUN apt-get update
    RUN apt-get install -qy libsodium-dev libcairo2-dev texlive-latex-recommended texlive-fonts-recommended libxml2-dev texlive-latex-extra libnspr4
    RUN apt-get update

	RUN apt-get clean

	RUN if [ ! -d "/opt/biotools" ];then mkdir /opt/biotools; fi
	RUN if [ ! -d "/opt/biotools/bin" ];then mkdir /opt/biotools/bin; fi
	RUN chmod 777 -R /opt/biotools/
	ENV PATH /opt/biotools/bin:$PATH
    
	RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 \
    && tar -jxvf samtools-1.17.tar.bz2 \
    && cd samtools-1.17 \
    &&  ./configure --without-curses --disable-lzma \
    && make -j 8 \
    && make install \
	&& cd .. && rm -rf samtools-1.17*

	RUN git clone https://github.com/vcftools/vcftools.git \
    && cd vcftools \
    &&  ./autogen.sh \
	&& ./configure \
    && make -j 8\
    && make install \
	&& cd .. && rm -rf vcftools*

	#httpuv dont a besoin d3scatter a besoin d'automake 1.16.1 !!
	RUN wget https://ftp.gnu.org/gnu/automake/automake-1.16.1.tar.gz \
    && tar -zxvf automake-1.16.1.tar.gz \
    && cd automake-1.16.1 \
    && ./configure --prefix=/opt/aclocal-1.16.1 \
    && make && make install \
	&& cd .. \
	&& rm -f automake-1.16.1.tar.gz

	ENV PATH "/opt/aclocal-1.16.1/bin:${PATH}" 

	RUN	Rscript -e 'install.packages("calibrate",repos="https://cloud.r-project.org/",Ncpus=8, clean=TRUE)'
	RUN Rscript -e 'install.packages("DT",Ncpus=8, repos="https://cloud.r-project.org/")'	
	RUN Rscript -e 'install.packages("shiny",Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages("shinydashboard",Ncpus=8, repos="https://cloud.r-project.org/")'
    RUN Rscript -e 'install.packages("pairsD3",Ncpus=8, repos="https://cloud.r-project.org/")'
    RUN Rscript -e 'install.packages("ggplot2",Ncpus=8, repos="https://cloud.r-project.org/")'

    RUN Rscript -e 'install.packages("threejs",Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages(c("gridExtra","gtable","tidyr"),dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages(c("ggalt"),dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages(c("Rtsne"),dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages("shinyFiles",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'install.packages("vroom",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'

    RUN Rscript -e 'install.packages("BiocManager",Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'BiocManager::install("tximport",Ncpus=8)'
	RUN Rscript -e 'BiocManager::install("GenomicFeatures",Ncpus=8)'
	RUN Rscript -e 'BiocManager::install("apeglm",Ncpus=8)'
	RUN Rscript -e 'BiocManager::install("BiocParallel",Ncpus=8)'
	RUN	Rscript -e 'BiocManager::install("SNPRelate",Ncpus=8, clean=TRUE)'
	RUN Rscript -e 'install.packages("caTools",Ncpus=8, repos="https://cloud.r-project.org/")'
	RUN Rscript -e 'BiocManager::install("gplots",Ncpus=8, clean=TRUE)'
	RUN Rscript -e 'BiocManager::install("ggtree",Ncpus=8, clean=TRUE)'

    RUN apt-get install -y libgit2-dev

 	RUN Rscript -e 'remotes::install_github("rstudio/crosstalk", quiet=T)'
    RUN Rscript -e 'install.packages("httpuv",Ncpus=8, repos="https://cloud.r-project.org/")'
    RUN Rscript -e 'remotes::install_github("jcheng5/d3scatter", quiet=F)' 
	RUN Rscript -e 'remotes::install_github("andrewsali/shinycssloaders", quiet=T)'
    

	RUN mkdir /samples
	RUN	mkdir /results
	RUN	mkdir /references

	RUN apt-get update && apt-get install -y python3
    RUN apt-get install -y python3-pip  python3-scipy python3-matplotlib
	RUN pip3 install numpy --upgrade
    RUN pip3 install cython

    RUN wget http://gnu.mirror.vexxhost.com/gsl/gsl-2.4.tar.gz && tar -zxvf gsl-2.4.tar.gz \
    && cd gsl-2.4 \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -R gsl-2.4.tar.gz gsl-2.4 \
	&& ln -s /usr/local/lib/libgsl* /usr/lib/

	# RUN cd /opt/biotools/bin  \
	# && pip3 uninstall -y Cython; pip3 install Cython==0.27.3 \
    # && git clone https://github.com/jashapiro/fastStructure.git  \
    # && cd /opt/biotools/bin/fastStructure && git checkout py3  \
    # && cd /opt/biotools/bin/fastStructure/vars  \
    # && python3 setup.py build_ext -f --inplace  \
    # && cd /opt/biotools/bin/fastStructure  \
    # && python3 setup.py build_ext -f  --inplace  \
    # && sed -i '2iimport matplotlib as mpl' /opt/biotools/bin/fastStructure/distruct.py \
    # && sed -i '3impl.use(\"svg\")' /opt/biotools/bin/fastStructure/distruct.py  \
    # && echo export 'PATH=/opt/biotools/bin/fastStructure:$PATH' >> /etc/environment &&\
	
    RUN Rscript -e "remotes::install_github('royfrancis/pophelper', Ncpus=8, upgrade ='never')" \
	# && pip3 uninstall -y Cython; pip3 install Cython

 

# ENV PATH /opt/biotools/bin/fastStructure:$PATH
# RUN cd /opt/biotools/bin \
# && wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20190304.zip \
# && unzip plink_linux_x86_64_20190304.zip \
# && rm -f plink_linux_x86_64_20190304.zip

RUN apt-get install -y libatlas-base-dev
# RUN cd /opt/biotools/bin \
# && git clone https://github.com/chrchang/plink-ng.git \
# && cd plink-ng/1.9 \
# && sed -i s'/MAX_CHROM_TEXTNUM 95/MAX_CHROM_TEXTNUM 950/' plink_common.h \
# && ./plink_first_compile \
# && cp plink /opt/biotools/bin \
# && cd ../.. \
# && rm -rf plink-ng
RUN apt-get update
RUN apt-get install -y imagemagick 
#RUN apt-get install -y python-backports.functools-lru-cache
RUN pip install backports.functools-lru-cache

RUN cd /opt/biotools/bin \
&& wget \https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 \
&& bunzip2 bcftools-1.17.tar.bz2 && tar -xvf bcftools-1.17.tar \
&& cd bcftools-1.17 && make && make install \
&& cd .. && rm -rf bcftools-1.17

RUN sed -i 's|11/2.54|20/2.54|' /usr/local/bin/plot-vcfstats
RUN sed -i 's|10/2.54|14/2.54|' /usr/local/bin/plot-vcfstats
RUN sed -i 's|window_len/2|window_len//2|g' /usr/local/bin/plot-vcfstats



RUN cd /opt/biotools && git clone https://bitbucket.org/gutenkunstlab/dadi.git \
&& cd dadi && python3 setup.py install


RUN apt-get install -y pandoc
RUN apt-get install -y python3-pandas

RUN Rscript -e 'install.packages("reticulate",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
RUN Rscript -e 'install.packages("shinyalert",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'

RUN pip3 install demes
RUN cd /opt/biotools && git clone https://bitbucket.org/simongravel/moments.git 
ADD local_scipy_dual_anneal.py /opt/biotools/moments/moments

ADD moments_inference_dualanneal.py /opt/biotools/moments/moments
RUN cd /opt/biotools/moments pip3 install -r requirements.txt \
&& python3 setup.py build_ext -i && python3 setup.py install 

RUN Rscript -e 'remotes::install_github("yonicd/snapper") '
RUN R --slave -e "install.packages(c('webshot', 'kableExtra','pander'));webshot::install_phantomjs()" 
RUN cp /root/bin/phantomjs /usr/local/bin/ 

#RUN Rscript -e 'BiocManager::install(c("SeqArray", "Rsamtools"))'
#Por
RUN Rscript -e 'remotes::install_github("zhengxwen/gdsfmt"); remotes::install_github("zhengxwen/SeqArray")'
RUN Rscript -e 'install.packages("pheatmap",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'

RUN apt-get install libxt-dev
RUN Rscript -e 'remotes::install_github("jokergoo/circlize", Ncpus=8)';
RUN Rscript -e 'remotes::install_github("jokergoo/ComplexHeatmap", Ncpus=8)'
RUN Rscript -e 'install.packages("shinyWidgets",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
RUN Rscript -e 'install.packages("patchwork",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
RUN Rscript -e 'install.packages("echarts4r",dependencies=T,Ncpus=8, repos="https://cloud.r-project.org/")'
#in case precedent ggtree installation failure
RUN Rscript -e 'install.packages("remotes"); remotes::install_github("YuLab-SMU/yulab.utils"); BiocManager::install("ggtree",Ncpus=8, clean=TRUE) '

#see : https://github.com/rajanil/fastStructure/issues/56 and https://structure-threader.readthedocs.io/en/latest/
RUN python3 -m pip install structure_threader --user

RUN pip3 install psutil
RUN mkdir sagApp
COPY app.R /sagApp/
COPY finalReport.Rmd /sagApp/
COPY local_Misc.py /sagApp/

ADD dadimodels /sagApp/dadimodels
ADD momentsmodels /sagApp/momentsmodels
ADD modelspages /sagApp/modelspages
ADD www /sagApp/www

#Cleaning
RUN rm -rf /opt/biotools/dadi
RUN rm -rf /opt/biotools/moments
#RUN pip3 install demes


EXPOSE 3838
ENTRYPOINT [ "Rscript", "-e", "Sys.setenv(BASH_ENV='/etc/environment'); setwd('/sagApp/'); shiny::runApp('/sagApp/app.R', port=3838 , host='0.0.0.0')"]
# To build
#docker build -t shinyvcfmultisampleanalyser .

#bind local directories to /Data and /Results in the container
#DOCK_VOLUME="--mount type=bind,src=/home/khalid/projets/workspace/testData,dst=/Data --mount type=bind,src=/tmp/Results,dst=/Results"

#Run with binding and mapping host tcp port 90 to port 3838 in the container (named vcfmultisample)
#docker run -d -p 90:3838 $DOCK_VOLUME --name vcfmultisample shinyvcfmultisampleanalyser

#The running container can be accessed in command line this way :
# docker exec -i -t  vcfmultisample /bin/bash
