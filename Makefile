#############################################
#   D R O P S   top-level makefile          #
#############################################
DROPS_ROOT = .

# include settings from the config file drops.conf:
include drops.conf

# variables:

PARPACKAGES = parallel levelset partests poisson DiST
SERPACKAGES = geom num out misc poisson stokes navstokes tests levelset surfactant transport
PACKAGES = $(SERPACKAGES) $(PARPACKAGES)
BUILDPACKAGES = $(if $(PAR_BUILD),$(PARPACKAGES),$(SERPACKAGES))

# rules:

default: dep all

all: $(BUILDPACKAGES:%=all_%)
	@echo "--> All executables generated successfully!"

strip: $(PACKAGES:%=strip_%)

clean: $(PACKAGES:%=clean_%)
	@echo "--> All object files and executables removed!"

distclean: $(PACKAGES:%=distclean_%) distclean_dox
	rm -f $(DEPFILE)
	@echo "--> Everything is cleaned up now!"

dep: deldepend $(PACKAGES:%=depend_%)
	@echo "--> Actual dependencies generated in $(DEPFILE)!"

doc:
	doxygen dox.cfg

stat:
	-ls $(PACKAGES:%=%/*.[ct]pp) $(PACKAGES:%=%/*.h) > $(DROPS_ROOT)/srcFiles.tmp
	cat $(DROPS_ROOT)/srcFiles.tmp | xargs wc --lines
	@rm -f $(DROPS_ROOT)/srcFiles.tmp
	@echo "--> lines of code (LOC) statistics done"

all_%:
	cd $* && $(MAKE) all

strip_%:
	cd $* && $(MAKE) strip

clean_%:
	cd $* && $(MAKE) clean

distclean_%:
	cd $* && $(MAKE) distclean

distclean_dox:
	cd ./doc && rm -rf dox

clean_HYPRE:
	cd $(HYPRE_HOME) && gmake clean && cd $(DROPS_ROOT)

topo:
	cd ./geom && $(MAKE) topo.cpp
	@echo "--> topo.cpp generated!"

deldepend:
	cp -f Dep.in $(DEPFILE)

depend_%:
	cd $* && \
        $(DEPEND) -- $(CXXFLAGS) -- -s"# $* dependencies:" -f- ../$*/*.cpp >> ../$(DEPFILE) 2>/dev/null; \
        echo " " >> ../$(DEPFILE)

prog_%:
	cd $(@D) && $(MAKE) $(*F)

HYPRE:	
	cd $(HYPRE_HOME) && gmake install

.PHONY: all clean distclean distclean_dox default dep deldepend doc stat topo check

	
