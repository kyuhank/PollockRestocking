IMAGE := dragonflyscience/dragonverse-kyuhan:v2

#all: ModelResults.RData SCRA_output.html

all: ModelResults.RData

ModelResults.RData: assessment.R
	Rscript assessment.R

#SCRA_output.html: SCRA_output.Rmd ModelResults.RData
#	Rscript -e "rmarkdown::render('SCRA_output.Rmd')"

docker:
	docker build -t $(IMAGE) .
	docker push $(IMAGE)

local:
	docker run -it --rm -w /work -v "$$PWD":/work \
    --net=host \
    -e DISPLAY=$$DISPLAY \
    -e RUN= \
    -v /tmp/.X11-unix/X0:/tmp/.X11-unix/X0 \
    -u $$(id -u):$$(id -g) ${IMAGE} bash
