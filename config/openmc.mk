$(OPENMC_BUILDDIR)/Makefile: $(OPENMC_DIR)/CMakeLists.txt
	mkdir -p $(OPENMC_BUILDDIR)
	mkdir -p $(LIBMESH_DIR)/include/contrib
	cd $(OPENMC_BUILDDIR) && \
	cmake3 -L \
	-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
	-Dlibmesh=ON \
	-Doptimize=ON \
	-DCMAKE_PREFIX_PATH=$(LIBMESH_DIR) \
	-DCMAKE_INSTALL_PREFIX=$(OPENMC_INSTALL_DIR) \
	-DCMAKE_INSTALL_LIBDIR=$(OPENMC_LIBDIR) \
	$(OPENMC_DIR)

build_openmc: | $(OPENMC_BUILDDIR)/Makefile
	make VERBOSE=1 -C $(OPENMC_BUILDDIR) install

cleanall_openmc: | $(OPENMC_BUILDDIR)/Makefile
	make -C $(OPENMC_BUILDDIR) uninstall clean

clobber_openmc:
	rm -rf $(OPENMC_LIB) $(OPENMC_BUILDDIR) $(OPENMC_INSTALL_DIR)

cleanall: cleanall_openmc

clobberall: clobber_openmc

.PHONY: build_openmc cleanall_openmc clobber_openmc
