FROM dynverse/dynwrap:py2.7


# install GPy
RUN git clone https://github.com/SheffieldML/GPy.git ; \
    cd GPy; \
    python setup.py develop; \
    cd ..

# install topslam
RUN git clone https://github.com/mzwiessele/topslam.git ; \
    cd topslam; \
    python setup.py develop

ADD . /code
ENTRYPOINT python /code/run.py
