FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'



RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/dynutils')"
RUN R -e "setRepositories(ind=1:2);install.packages('reshape2')"
RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/destiny')"

ADD . /code

ENTRYPOINT Rscript code/run.R