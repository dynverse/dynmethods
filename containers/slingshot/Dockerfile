FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'
RUN R -e 'devtools::install_github("dynverse/dyndimred")'

ADD . /code

ENTRYPOINT Rscript /code/run.R