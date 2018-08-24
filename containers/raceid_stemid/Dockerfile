FROM dynverse/dynwrap:bioc

RUN R -e 'install.packages("destiny")'

RUN apt-get install -y libcgal-dev libglu1-mesa-dev libglu1-mesa-dev

RUN R -e 'install.packages("FateID")'

RUN R -e 'install.packages("RaceID")'

ADD . /code

ENTRYPOINT /code/run.sh
