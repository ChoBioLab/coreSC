FROM coresc-test

ENV PARAMS_HOST=/home/chris/Downloads/config/params.csv
ENV SAMPLES_HOST=/home/chris/Downloads/config/params.csv
ENV OUT_HOST=/home/chris/Downloads/output

ENV PARAMS_CONTAINER=/app/config/params.csv
ENV SAMPLES_CONTAINER=/app/config/samples.csv
ENV OUT_CONTAINER=/app/output

VOLUME "${PARAMS_HOST}:${PARAMS_CONTAINER}"
VOLUME "${SAMPLES_HOST}:${SAMPLES_CONTAINER}"
VOLUME "${OUT_HOST}:${OUT_CONTAINER}"
