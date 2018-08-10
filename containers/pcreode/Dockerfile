FROM dynverse/dynwrap:py2.7

RUN git clone https://github.com/KenLauLab/pCreode.git && pip install pCreode

ADD . /code

ENTRYPOINT python /code/run.py
