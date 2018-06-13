FROM rocker/tidyverse

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e "setRepositories(ind=1:2);devtools::install_github('kstreet13/slingshot')"
RUN R -e "setRepositories(ind=1:2);devtools::install_github('dynverse/dyndimred')"

ADD . /code

ENTRYPOINT Rscript code/run.R