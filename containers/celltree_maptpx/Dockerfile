FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN apt-get install -y libgsl-dev

RUN R -e 'devtools::install_cran("cellTree")'

ADD . /code

ENTRYPOINT /code/run.sh
