FROM dynverse/dynwrap:py3.6

RUN pip install python-igraph louvain # igraph and louvain do not get installed by scanpy

RUN pip install feather-format

RUN pip install scanpy

ADD . /code

ENTRYPOINT python /code/run.py
