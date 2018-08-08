FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'

RUN apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("rcannood/gng")'
RUN R -e 'install.packages("igraph")'
RUN R -e 'devtools::install_github("dynverse/dyndimred")'

ADD . /code

ENTRYPOINT Rscript /code/run.R