version = 2.0.6
os = alpine
osversion = 3.11

all:
	docker build -t overcookedfrog/stringtie:$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/stringtie:$(version)-$(os)$(osversion) overcookedfrog/stringtie:$(version)
	docker tag overcookedfrog/stringtie:$(version)-$(os)$(osversion) overcookedfrog/stringtie:latest

.PHONY: push clean
push:
	docker push overcookedfrog/stringtie:$(version)-$(os)$(osversion)
	docker push overcookedfrog/stringtie:$(version)
	docker push overcookedfrog/stringtie:latest
	
clean:
	docker rmi overcookedfrog/stringtie:$(version)-$(os)$(osversion)
	docker rmi overcookedfrog/stringtie:$(version)
	docker rmi overcookedfrog/stringtie:latest
