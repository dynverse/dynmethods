FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages("BiocGenerics")'

RUN R -e 'install.packages("SingleCellExperiment")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
