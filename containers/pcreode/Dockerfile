FROM dynverse/dynwrap:py2.7

RUN git clone git://github.com/KenLauLab/pCreode && pip install pCreode

ADD . /code

ENTRYPOINT python /code/run.py
