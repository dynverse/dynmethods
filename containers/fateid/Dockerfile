FROM dynverse/dynwrap:bioc

RUN R -e 'devtools::install_cran("destiny")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_cran("FateID")'

ADD . /code

ENTRYPOINT /code/run.sh
