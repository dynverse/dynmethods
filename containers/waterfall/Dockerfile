FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN R -e "devtools::install_github('dynverse/Waterfall')"

ADD . /code

ENTRYPOINT Rscript code/run.R