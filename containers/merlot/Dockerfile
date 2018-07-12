FROM rocker/tidyverse

# install rgl
RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

# install python
RUN apt-get install -y python3 python3-tk python3-pip
RUN apt-get install -y python3-scipy python3-numpy python3-pandas
RUN pip3 install cython
RUN pip3 install git+https://github.com/soedinglab/csgraph_mod

# install udunits (for elpigraph)
RUN apt-get -y install libudunits2-dev
RUN Rscript -e 'install.packages("udunits2", configure.args =  c(udunits2 = "--with-udunits2-include=/usr/include/udunits2"))'

# install elpigraph
RUN R -e "setRepositories(ind=1:2);devtools::install_github('Albluca/ElPiGraph.R')"

RUN R -e "print('hi')"

# install merlot
RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/merlot')"

# install destiny
RUN R -e "setRepositories(ind=1:4);install.packages('destiny')"

ADD . /code

ENTRYPOINT Rscript /code/run.R
