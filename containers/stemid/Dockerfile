FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit","vegan"))'

RUN git clone https://github.com/dgrun/StemID.git

ADD . /code

ENTRYPOINT Rscript /code/run.R
