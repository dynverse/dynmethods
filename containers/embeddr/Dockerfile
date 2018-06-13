FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/scaterlegacy')"
RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/embeddr')"

ADD . /code

ENTRYPOINT Rscript code/run.R