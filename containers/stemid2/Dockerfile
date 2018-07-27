FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'install.packages(c("tsne","pheatmap","MASS","cluster","mclust","flexmix","lattice","fpc","RColorBrewer","permute","amap","locfit","vegan"))'

RUN R -e 'install.packages(c("Rtsne","scran","randomForest","rgl", "scran", "DESeq2", "biomaRt"))'

RUN git clone https://github.com/dgrun/RaceID3_StemID2.git

ADD . /code

ENTRYPOINT Rscript /code/run.R
