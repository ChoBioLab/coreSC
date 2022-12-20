FROM public.ecr.aws/chobiolab/seurat:latest

WORKDIR /app
ADD . /app

MAINTAINER christopher.tastad@mssm.edu

ENTRYPOINT ["/app/main"]
