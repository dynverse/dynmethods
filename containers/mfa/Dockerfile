FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "setRepositories(ind=1:2);devtools::install_github('kieranrcampbell/mfa')"
RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/dyndimred')"

ADD . /code

ENTRYPOINT Rscript code/run.R