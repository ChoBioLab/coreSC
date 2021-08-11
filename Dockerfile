# syntax=docker/dockerfile:1
FROM r-base
LABEL author="chris tastad"
LABEL email="christopher.tastad@mssm.edu"
COPY . /app
RUN apt-get update
RUN Rscript /app/scripts/install_packages.R
CMD ["/app/run"]
