#pragma once

typedef  Mesh::Traits::Kernel Kernel;
typedef  Kernel::Point_3 Point;
typedef  Kernel::Point_3 Point_3; 
typedef  Kernel::Point_2 Point_2; 
typedef  Kernel::Vector_3 Vector;
typedef  Kernel::Vector_3 Vector_3; 

typedef  Mesh::HalfedgeDS HDS;
typedef  Mesh::Vertex_handle Vertex_handle; 
typedef  Mesh::Vertex_const_handle Vertex_const_handle; 
typedef  Mesh::Halfedge_handle Halfedge_handle; 
typedef  Mesh::Halfedge_const_handle Halfedge_const_handle; 
typedef  Mesh::Facet_handle Facet_handle; 
typedef  Mesh::Facet_const_handle Facet_const_handle; 

typedef  Mesh::Vertex_iterator Vertex_iterator; 
typedef  Mesh::Vertex_const_iterator Vertex_const_iterator; 
typedef  Mesh::Facet_iterator Facet_iterator; 
typedef  Mesh::Facet_const_iterator Facet_const_iterator; 
typedef  Mesh::Halfedge_iterator Halfedge_iterator; 
typedef  Mesh::Halfedge_const_iterator Halfedge_const_iterator; 
typedef  Mesh::Edge_iterator Edge_iterator; 
typedef  Mesh::Edge_const_iterator Edge_const_iterator; 

typedef  Mesh::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
typedef  Mesh::Halfedge_around_vertex_circulator  Halfedge_around_vertex_circulator;
typedef  Mesh::Halfedge_around_facet_const_circulator Halfedge_around_facet_const_circulator;
typedef  Mesh::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;


typedef  CGAL::Parameterization_polyhedron_adaptor_3<Mesh> MeshAdaptor;
typedef  CGAL::Parameterization_mesh_feature_extractor<MeshAdaptor> MeshFeatureExtractor;
typedef  MeshFeatureExtractor::Skeleton Skeleton;
typedef  MeshFeatureExtractor::Border Border;




/*
typedef CGAL::Exact_predicates_exact_constructions_kernel		ARRKernel;
typedef ARRKernel::FT											ARRNumberType;
typedef CGAL::Arr_segment_traits_2<ARRKernel>					ARRTraits_2;
typedef ARRTraits_2::Point_2									EPoint_2;
typedef ARRKernel::Vector_2										EVector_2;
typedef ARRTraits_2::X_monotone_curve_2							ESegment_2;
typedef ARRTraits_2::Ray_2										ERay_2;
typedef CGAL::Arr_extended_dcel<ARRTraits_2,
								Mesh::Vertex_handle,
								Mesh::Halfedge_handle,
								Mesh::Facet_handle>				Dcel;
typedef CGAL::Arrangement_2<ARRTraits_2, Dcel>					Arrangement_2;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2>		Landmarks_pl;
typedef CGAL::Arr_default_overlay_traits<Arrangement_2>			Overlay_traits;
typedef CGAL::Polygon_2<ARRKernel>								EPolygon_2;
*/


///////////////////////////////////////////////////////////////////////////////////////////////////////////////


#define for_each_vertex(vi, m) for(Mesh::Vertex_iterator vi = (m).vertices_begin(); vi != (m).vertices_end(); ++vi)
#define for_each_const_vertex(vi, m) for(Mesh::Vertex_const_iterator vi = (m).vertices_begin(); vi != (m).vertices_end(); ++vi)
#define for_each_halfedge(hi, m) for(Mesh::Halfedge_iterator hi = (m).halfedges_begin(); hi != (m).halfedges_end(); ++hi)
#define for_each_const_halfedge(hi, m) for(Mesh::Halfedge_const_iterator hi = (m).halfedges_begin(); hi != (m).halfedges_end(); ++hi)
#define for_each_edge(ei, m) for(Mesh::Edge_iterator ei = (m).edges_begin(); ei != (m).edges_end(); ++ei)
#define for_each_const_edge(ei, m) for(Mesh::Edge_const_iterator ei = (m).edges_begin(); ei != (m).edges_end(); ++ei)
#define for_each_facet(fi, m) for(Mesh::Facet_iterator fi = (m).facets_begin(); fi != (m).facets_end(); ++fi)
#define for_each_const_facet(fi, m) for(Mesh::Facet_const_iterator fi = (m).facets_begin(); fi != (m).facets_end(); ++fi)
