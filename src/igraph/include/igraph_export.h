
#ifndef IGRAPH_EXPORT_H
#define IGRAPH_EXPORT_H

#ifdef IGRAPH_STATIC
#  define IGRAPH_EXPORT
#  define IGRAPH_NO_EXPORT
#else
#  ifndef IGRAPH_EXPORT
#    ifdef igraph_EXPORTS
        /* We are building this library */
#      define IGRAPH_EXPORT 
#    else
        /* We are using this library */
#      define IGRAPH_EXPORT 
#    endif
#  endif

#  ifndef IGRAPH_NO_EXPORT
#    define IGRAPH_NO_EXPORT 
#  endif
#endif

/* The deprecation attribute is compiler specific.  CMake's export-header
 * generator always emits __declspec(deprecated), which is the MSVC spelling and
 * is rejected by GCC and Clang on Linux and macOS.  gAtoms ships a single copy
 * of this file for every platform, so the attribute is selected here instead of
 * being taken verbatim from one host's CMake run.  The deprecation is still
 * expressed as a real compiler attribute on every platform; it is never defined
 * away to an empty macro. */
#ifndef IGRAPH_DEPRECATED
#  if defined(_MSC_VER)
#    define IGRAPH_DEPRECATED __declspec(deprecated)
#  elif defined(__GNUC__) || defined(__clang__)
#    define IGRAPH_DEPRECATED __attribute__((deprecated))
#  else
#    define IGRAPH_DEPRECATED
#  endif
#endif

#ifndef IGRAPH_DEPRECATED_EXPORT
#  define IGRAPH_DEPRECATED_EXPORT IGRAPH_EXPORT IGRAPH_DEPRECATED
#endif

#ifndef IGRAPH_DEPRECATED_NO_EXPORT
#  define IGRAPH_DEPRECATED_NO_EXPORT IGRAPH_NO_EXPORT IGRAPH_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef IGRAPH_NO_DEPRECATED
#    define IGRAPH_NO_DEPRECATED
#  endif
#endif

#endif /* IGRAPH_EXPORT_H */
