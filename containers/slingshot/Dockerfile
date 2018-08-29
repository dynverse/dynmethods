FROM dynverse/dynwrap:bioc

LABEL version 0.1.0.1

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

ADD . /code

ENTRYPOINT /code/run.sh
