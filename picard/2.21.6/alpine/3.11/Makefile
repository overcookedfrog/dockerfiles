version = 2.21.6
os = alpine
osversion = 3.11

all: docker

docker: Dockerfile
	docker build -t overcookedfrog/picard:$(version)-$(os)$(osversion) .
	docker tag overcookedfrog/picard:$(version)-$(os)$(osversion) overcookedfrog/picard:$(version)
	docker tag overcookedfrog/picard:$(version)-$(os)$(osversion) overcookedfrog/picard:latest

.PHONY: push clean

push:
	docker push overcookedfrog/picard:$(version)-$(os)$(osversion)
	docker push overcookedfrog/picard:$(version)
	docker push overcookedfrog/picard:latest

clean:
	docker rmi overcookedfrog/picard:$(version)-$(os)$(osversion)
	docker rmi overcookedfrog/picard:$(version)
	docker rmi overcookedfrog/picard:latest
