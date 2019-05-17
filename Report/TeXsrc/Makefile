JobName=Parareal
TeX=pdflatex
LatexFlags= -jobname=$(JobName) --shell-escape

default:
	make build
	make build
	make clean
	cp $(JobName).pdf ../

build:
	$(TeX) $(LatexFlags) main.tex

clean:
	rm -f $(JobName).log $(JobName).aux $(JobName).out $(JobName).toc
	rm -rf _minted-*
