version = 1.9
os = alpine
osversion = 3.11

all:
	docker build -t overcookedfrog/bcftools:$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/bcftools:$(version)-$(os)$(osversion) overcookedfrog/bcftools:$(version)
	docker tag overcookedfrog/bcftools:$(version)-$(os)$(osversion) overcookedfrog/bcftools:latest

.PHONY: push
push:
	docker push overcookedfrog/bcftools:$(version)-$(os)$(osversion)
	docker push overcookedfrog/bcftools:$(version)
	docker push overcookedfrog/bcftools:latest
