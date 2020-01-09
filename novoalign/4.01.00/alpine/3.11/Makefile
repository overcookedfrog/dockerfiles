name = novoalign
version = 4.01.00
os = alpine
osversion = 3.11
dockeruser = overcookedfrog

all: docker

docker: Dockerfile
	docker build -t $(dockeruser)/$(name):$(version)-$(os)$(osversion) .
	docker tag $(dockeruser)/$(name):$(version)-$(os)$(osversion) $(dockeruser)/$(name):$(version)
	docker tag $(dockeruser)/$(name):$(version)-$(os)$(osversion) $(dockeruser)/$(name):latest

.PHONY: push clean

push:
	docker push $(dockeruser)/$(name):$(version)-$(os)$(osversion)
	docker push $(dockeruser)/$(name):$(version)
	docker push $(dockeruser)/$(name):latest

clean:
	docker rmi $(dockeruser)/$(name):$(version)-$(os)$(osversion)
	docker rmi $(dockeruser)/$(name):$(version)
	docker rmi $(dockeruser)/$(name):latest



