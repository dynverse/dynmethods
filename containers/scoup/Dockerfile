FROM dynverse/dynwrap:r

LABEL version 0.1.0.1

RUN git clone https://github.com/hmatsu1226/SCOUP.git && cd SCOUP && make all

ADD . /code

ENTRYPOINT /code/run.sh
