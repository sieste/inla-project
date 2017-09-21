inla-from-scratch: inla_from_scratch.md

inla_from_scratch.md: inla_from_scratch.Rmd
	Rscript -e 'knitr::knit("inla_from_scratch.Rmd")'


