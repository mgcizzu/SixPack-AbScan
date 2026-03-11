FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_NO_CACHE_DIR=1
ENV PORT=7860
ENV HOST=0.0.0.0
ENV USER=appuser
ENV HOME=/home/appuser
ENV GRADIO_TEMP_DIR=/home/appuser/app/temp

RUN useradd -m -u 1000 $USER

WORKDIR $HOME/app

COPY requirements.txt $HOME/app/requirements.txt
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . $HOME/app

RUN mkdir -p $HOME/app/runs $GRADIO_TEMP_DIR \
    && chown -R $USER:$USER $HOME

USER $USER

EXPOSE 7860

CMD ["python", "app_gradio.py"]
