FROM rocker/tidyverse

# install rgl
RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

# install STEMNET
RUN R -e "devtools::install_github('dgrun/FateID')"

ADD . /code

ENTRYPOINT Rscript code/run.R
