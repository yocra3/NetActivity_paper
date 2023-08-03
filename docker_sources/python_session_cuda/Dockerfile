FROM tensorflow/tensorflow:2.7.0-gpu
LABEL authors="Carlos Ruiz-Arenas" \
      description="Docker image containing the python packages for Episignatures project"

## Install python dependencies
RUN pip install -U scikit-learn matplotlib scipy pandas

RUN apt-get install -y git

RUN git clone https://github.com/arnovanhilten/GenNet /usr/GenNet/
RUN pip install -q tensorflow-model-optimization
