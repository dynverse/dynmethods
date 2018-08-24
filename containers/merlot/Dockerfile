FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("destiny")'

# install rgl dependencies
RUN apt-get update
RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

# install python
RUN apt-get install -y python3 python3-tk python3-pip
RUN apt-get install -y python3-scipy python3-numpy python3-pandas
RUN pip3 install cython
RUN pip3 install git+https://github.com/soedinglab/csgraph_mod

# install udunits (for elpigraph)
RUN apt-get -y install libudunits2-dev
RUN Rscript -e 'install.packages("udunits2", configure.args =  c(udunits2 = "--with-udunits2-include=/usr/include/udunits2"))'

# install merlot
RUN R -e "devtools::install_github('soedinglab/merlot')"

ADD . /code

ENTRYPOINT /code/run.sh
