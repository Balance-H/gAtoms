# R CMD INSTALL invokes this in src, including during R CMD check.
#
# gAtoms links its bundled, statically built igraph (src/igraph, archived as
# src/igraph/libigraph.a) directly into its own shared object, so exactly one
# file is installed here: gAtoms.dll (gAtoms.so on Unix-like platforms). No
# external igraph runtime or third-party DLL closure is copied any more; the
# former recursive runtime-copy step, which required an environment variable
# pointing at an external igraph DLL directory, has been removed.
dll <- paste0("gAtoms", .Platform$dynlib.ext)
dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
if (!file.copy(dll, dest, overwrite = TRUE)) stop("Cannot install gAtoms DLL")

# Preserve object-level symbols of gAtoms' own sources so that R CMD check can
# distinguish this package's code from compiler startup code and from the
# statically linked igraph code. The object files exist in a development build
# (pkgload), but a release installation removes them; read_symbols_from_object_file
# returns NULL for an absent file, so this degrades gracefully instead of failing.
objects <- c("decom_h.o", "r_interface.o")
present <- objects[file.exists(objects)]
tables <- lapply(present, tools:::read_symbols_from_object_file)
names(tables) <- present
saveRDS(tables, file.path(dest, "symbols.rds"), version = 2)
