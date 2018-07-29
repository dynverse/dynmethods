FROM library/python
RUN pip install pandas numpy matplotlib

RUN pip install python-igraph louvain # igraph and louvain do not get installed by scanpy
RUN pip install feather-format
RUN pip install scanpy
ADD . /code
ENTRYPOINT python /code/run.py
