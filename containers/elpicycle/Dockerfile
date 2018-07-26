FROM rocker/tidyverse

# install udunits
RUN apt-get -y install libudunits2-dev
RUN Rscript -e 'install.packages("udunits2", configure.args =  c(udunits2 = "--with-udunits2-include=/usr/include/udunits2"))'

# install elpigraph from github
RUN R -e "devtools::install_github('Albluca/ElPiGraph.R')"

ADD . /code
ENTRYPOINT Rscript /code/run.R
