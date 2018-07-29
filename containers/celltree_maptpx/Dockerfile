FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN apt-get install -y libgsl-dev

RUN R -e 'install.packages("cellTree")'

ADD . /code

ENTRYPOINT Rscript /code/run.R
