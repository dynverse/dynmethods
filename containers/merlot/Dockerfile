FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN R -e 'devtools::install_cran("destiny")'

RUN apt-get update && apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN apt-get install -y python3 python3-tk python3-pip
RUN apt-get install -y python3-scipy python3-numpy python3-pandas

RUN pip3 install cython

RUN pip3 install git+https://github.com/soedinglab/csgraph_mod

RUN apt-get -y install libudunits2-dev

RUN Rscript -e 'devtools::install_cran("udunits2", configure.args =  c(udunits2 = "--with-udunits2-include=/usr/include/udunits2"))'

RUN R -e "devtools::install_github('soedinglab/merlot')"

ADD . /code

ENTRYPOINT /code/run.sh
