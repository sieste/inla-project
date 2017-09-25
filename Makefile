all: global-temperature

inla-from-scratch: inla-from-scratch.md
global-temperature: global-temperature.md

%.md: %.Rmd
	Rscript -e 'knitr::knit("$<")'

