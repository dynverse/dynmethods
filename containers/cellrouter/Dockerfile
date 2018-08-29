FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN R -e 'devtools::install_cran(c("reshape", "reshape2", "pheatmap", "tsne", "igraph", "ggplot2", "mclust", "grid", "Rtsne", "cccd", "irlba"))'

RUN git clone https://github.com/edroaldo/cellrouter.git && find cellrouter -type f | grep -v "^cellrouter/CellRouter" | xargs rm

RUN apt-get update && apt-get install -y default-jre

ADD . /code
ENTRYPOINT /code/run.sh
