/* gAtoms compatibility shim.
 *
 * igraph 0.10.x shipped its public headers under include/igraph/, so downstream
 * code includes them as <igraph/igraph.h>. igraph 1.0.x flattened the layout:
 * the headers now live directly in include/.
 *
 * This forwarding header keeps the <igraph/igraph.h> include spelling working
 * without modifying any gAtoms source file. It contains no declarations of its
 * own; it only forwards to the vendored upstream header. The upstream headers
 * use quoted same-directory includes, which resolve because src/Makevars adds
 * -Iigraph/include to the search path.
 */
#include "../igraph.h"
