FROM dynverse/dynwrap:py2.7

LABEL version 0.1.0.1

RUN git clone https://github.com/KenLauLab/pCreode.git && pip install pCreode

ADD . /code

ENTRYPOINT /code/run.sh
