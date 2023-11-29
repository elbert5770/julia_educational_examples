//------------------------------------------------
// MarchingCubes
//------------------------------------------------
//
// MarchingCubes Algorithm
// Version 0.2 - 12/08/2002
//
// Thomas Lewiner thomas.lewiner@polytechnique.org
// Math Dept, PUC-Rio
//
// Translated to C by Ziad S. Saad November 30/04
//________________________________________________

#ifndef _MARCHINGCUBES_H_
#define _MARCHINGCUBES_H_

#ifdef NII2MESH
 #include "meshtypes.h"
 int marchingCubes(float * img, size_t dim[3], int lo[3], int hi[3], int originalMC, float isolevel, vec3d **vs, vec3i **ts, int *nv, int *nt);
#endif
//_____________________________________________________________________________
// types
#define false 0
#define true 1

typedef unsigned char uchar ;
typedef   signed char schar ;

// Vertex structure
typedef struct
{
  float  x,  y,  z ;  // Vertex' coordinates
  float nx, ny, nz ;  // Vertex' normal
} Vertex ;

// Triangle structure
typedef struct
{
  long v1,v2,v3 ;  // Vertex' indexes
} Triangle ;

//MarchingCubes structure
#define N_MAX  15
typedef struct
{
  int      originalMC ;

  int      size_x,size_y,size_z ;      // size of the size_z planes is size_x*size_y
  float    * data      ;                     // sampled function
  long      * x_verts, * y_verts, * z_verts ; // size of the size_z planes is size_x*size_y

  long      nverts,ntrigs ;               // number of allocated vertices and triangles in the arrays
  long      Nverts,Ntrigs ;               // sizes of the vertices' and triangle's arrays
  Vertex   * vertices  ;                     // array of all vertices
  Triangle * triangles ;                     // array of all triangles as triple of vertices

  long      i, j, k  ;                      // counters
  float    cube[8]  ;                      // vertex values
  int      N[N_MAX]    ;                      // case counters
  uchar    lut_entry;                      // cube configuration
  uchar    _case ,config,subconfig ;     // case and config entry   
}MCB;
//_____________________________________________________________________________

#if 0 /* Not needed */
long nverts(MCB *mcb)  { return (mcb->nverts) ; }
long ntrigs(MCB *mcb)  { return (mcb->ntrigs) ; }
int size_x(MCB *mcb)  { return (mcb->size_x) ; }
int size_y(MCB *mcb)  { return (mcb->size_y) ; }
int size_z(MCB *mcb)  { return (mcb->size_z) ; }

Vertex   *vertices (MCB *mcb) { return (mcb->vertices)  ; }
Triangle *triangles(MCB *mcb) { return (mcb->triangles) ; }

Vertex   * vert( MCB *mcb, const long i )  { if( i < 0  || i >= mcb->nverts ) return (( Vertex *)NULL) ; return (mcb->vertices  + i) ; }
Triangle * trig( MCB *mcb, const long i )  { if( i < 0  || i >= mcb->ntrigs ) return ((Triangle*)NULL) ; return (mcb->triangles + i) ; }

#endif


void set_resolution( MCB *mcb, int size_x,  int size_y,  int size_z );
void set_method    ( MCB *mcb, int originalMC ) ;

  // Data access
float get_data  (  MCB *mcb, long i,  long j,  long k )  ;
void  set_data  (  MCB *mcb, float val,  long i,  long j,  long k ) ;

// Data initialization
void init_temps (MCB *mcb ) ;
void init_all   (MCB *mcb ) ;
void clean_temps(MCB *mcb ) ;
void clean_all  (MCB *mcb ) ;

void run(MCB *mcb ) ;              // Main algorithm
void process_cube (MCB *mcb ) ;              // Process a unit cube
int test_interior( MCB *mcb , schar s )    ;  // Test the interior of a cube
int test_face    ( MCB *mcb , schar face ) ;  // Test a face

// Compute the intersection points
void compute_intersection_points(MCB *mcb ) ;

// Adding triangles
void add_triangle ( MCB *mcb , const char* trig, char n, int v12 ) ; /* default for v12 = -1 */

// Adding vertices
void test_vertex_addition(MCB *mcb ) ;
int add_x_vertex(MCB *mcb ) ;
int add_y_vertex(MCB *mcb ) ;
int add_z_vertex(MCB *mcb ) ;
int add_c_vertex(MCB *mcb ) ;

// Computing gradient
float get_x_grad( MCB *mcb ,  long i,  long j,  long k )  ;
float get_y_grad( MCB *mcb ,  long i,  long j,  long k )  ;
float get_z_grad( MCB *mcb ,  long i,  long j,  long k )  ;

// Pre-computed vertices access
long   get_x_vert(  MCB *mcb , long i,  long j,  long k )  ;
long   get_y_vert(  MCB *mcb , long i,  long j,  long k )  ;
long   get_z_vert(  MCB *mcb , long i,  long j,  long k )  ;
void  set_x_vert(  MCB *mcb , long val,  long i,  long j,  long k ) ;
void  set_y_vert(  MCB *mcb , long val,  long i,  long j,  long k ) ;
void  set_z_vert(  MCB *mcb , long val,  long i,  long j,  long k ) ;

// print cube for debug
void print_cube(MCB *mcb) ;

MCB *MarchingCubes (  int size_x ,  int size_y ,  int size_z ) ; /* defaults are -1 -1 -1 */
void FreeMarchingCubes (MCB *mcb);

void set_suma_debug(int dbg);

#endif // _MARCHINGCUBES_H_