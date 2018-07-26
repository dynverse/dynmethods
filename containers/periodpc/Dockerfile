FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/dynwrap")'
RUN R -e 'devtools::install_github("dynverse/dyndimred")'

RUN R -e 'install.packages("princurve")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
