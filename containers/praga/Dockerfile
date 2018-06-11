FROM dynverse/paga

RUN pip install pyaml

ADD overwrite_definition.py /code/

RUN python /code/overwrite_definition.py
