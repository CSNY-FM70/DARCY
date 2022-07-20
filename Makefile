tests:
	cd tests && ./buildTests.sh
apps:
	cd libdarcy_apps && $(MAKE)
clean:
	cd libdarcy_apps && $(MAKE) clean
doc:
	cd libdarcy_apps && $(MAKE) doc	
cleandoc:
	cd libdarcy_apps && $(MAKE) cleandoc

.PHONY: apps tests clean doc cleandoc
