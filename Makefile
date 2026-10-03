.PHONY: paper clean check

paper:
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex
	cd paper && pdflatex -interaction=nonstopmode -halt-on-error main.tex

check:
	python3 checks/verify_audit.py

clean:
	rm -f paper/*.aux paper/*.log paper/*.out paper/*.toc
