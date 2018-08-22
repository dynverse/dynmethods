FROM dynverse/dynwrap:bioc

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'devtools::install_github("kstreet13/slingshot")'

ADD . /code

ENTRYPOINT cat /code/run.R | R --no-save
