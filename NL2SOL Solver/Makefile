
Default:

CURDIRBASE := $(shell basename "$(CURDIR)")
TARGET := $(subst -,_,$(subst .,_,$(subst XLL,,$(CURDIRBASE))))

update_from_msvc_build:
	@mkdir -p x64
	@cp -v build/x64/Release/$(TARGET).xll x64/
	@mkdir -p Win32
	@cp -v build/Win32/Release/$(TARGET).xll Win32/

archive:
	cd .. ; mkdir -p 7z; ARCHIVE=7z/"$(TARGET)".$(shell date +%F).7z; mv -f $$ARCHIVE /tmp; rm -f $$ARCHIVE; 7z a -l -t7z -m0=lzma -xr\!.svn -xr\!.vs -xr\!build -xr\!__pycache__ $$ARCHIVE "$(CURDIRBASE)"
