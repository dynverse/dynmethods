FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN R -e "devtools::install_github('zji90/TSCAN')"
RUN R -e "install.packages('igraph')"

ADD . /code

ENTRYPOINT Rscript code/run.R