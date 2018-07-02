FROM rocker/tidyverse

RUN R -e "setRepositories(ind=1:2);install.packages('BiocGenerics')"
RUN R -e "install.packages('igraph')"

RUN R -e "setRepositories(ind=1:2);devtools::install_github('dcellwanger/CellTrails')"

ADD . /code
ENTRYPOINT Rscript /code/run.R
