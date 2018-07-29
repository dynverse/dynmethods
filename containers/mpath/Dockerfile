FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github("dynverse/Mpath")'
# RUN R -e 'devtools::install_url("https://github.com/JinmiaoChenLab/Mpath/raw/master/Mpath_1.0.tar.gz")'
# dynverse/Mpath is just so much easier to use than the tar gz above

ADD . /code

ENTRYPOINT Rscript /code/run.R
