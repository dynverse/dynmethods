FROM dynverse/dynwrap:py3.6

RUN pip install Cython

# install wishbone
RUN pip install git+https://github.com/dynverse/pywishbone --upgrade --upgrade-strategy only-if-needed

ADD . /code
ENTRYPOINT /code/run.sh
