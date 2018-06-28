FROM rocker/tidyverse

RUN R -e "install.packages('BiocGenerics')"
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('Biobase')"

RUN R -e "devtools::install_github('dcellwanger/CellTrails')"

ADD . /code
ENTRYPOINT Rscript /code/run.R
