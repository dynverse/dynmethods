FROM dynverse/dynwrap:r

LABEL version 0.1.0.1

RUN apt-get -y install libudunits2-dev

RUN Rscript -e 'devtools::install_cran("udunits2", configure.args =  c(udunits2 = "--with-udunits2-include=/usr/include/udunits2"))'

RUN R -e "devtools::install_github('Albluca/ElPiGraph.R')"

ADD . /code
ENTRYPOINT /code/run.sh
