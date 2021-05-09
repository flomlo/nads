#include <fstream>
#include <iostream>
#include <cstdint>

#define USE_CELLS_WITHOUT_DIMENSION

#include "include/argparser.h"
#include "include/parameters.h"
#include "include/usage/flagser-count.h"
#include "include/complex/directed_flag_complex_in_memory_computer.h" 
//#include "include/complex/directed_flag_complex_computer.h" 


typedef std::vector<std::pair<int, vertex_index_t>> coboundaries_t;
typedef std::vector<std::pair<vertex_index_t, vertex_index_t>> edges_t;
typedef struct{} NoneType;

struct ads_computer_t {
  const filtered_directed_graph_t& graph;
  int keep_prob; // in per mille
  //srand(0);
  edges_t edges{};
  coboundaries_t coboundaries{};  //have on coboundary vector per thread, reuse per simplex to save allocs

  void done() {}

  void operator()(vertex_index_t* first_vertex, int size, NoneType) {
    std::vector<size_t> vertex_offsets;
    for (int j = 0; j < size; j++) vertex_offsets.push_back(first_vertex[j] >> 6);

    coboundaries.clear();
    
    // find all coboundaries the the simplex given by 
    for (int i = 0; i <= size; i++) {
      // Check intersections in chunks of 64
      for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
        size_t bits = -1; // All bits set

        for (int j = 0; bits > 0 && j < size; j++) {
          // Remove the vertices already making up the cell
          if (vertex_offsets[j] == offset) bits &= ~(ONE_ << (first_vertex[j] - (vertex_offsets[j] << 6)));

          // Intersect with the outgoing/incoming edges of the current vertex
          bits &= j < i ? graph.get_outgoing_chunk(first_vertex[j], offset)
                        : graph.get_incoming_chunk(first_vertex[j], offset);
        }

        size_t vertex_offset = offset << 6;
        while (bits > 0) {
          // Get the least significant non-zero bit
          auto b = __builtin_ctzl(bits);

          // Unset this bit
          bits &= ~(ONE_ << b);

          coboundaries.push_back({i,vertex_offset+b});
        }
      }
    }
    // find almost-d-simplices
    for (size_t i = 0; i < coboundaries.size(); i++) {
      for (size_t j = i+1; j < coboundaries.size(); j++) {
        if (rand() % 1000 < keep_prob) {  // we drop here instead later for performance reasons
          auto cb1 = coboundaries[i];
          auto cb2 = coboundaries[j];
          if (cb1.second != cb2.second) {
            if (cb1.first <= cb2.first) edges.push_back({cb1.second, cb2.second});
            if (cb1.first >= cb2.first) edges.push_back({cb2.second, cb1.second});
          }
        }
      }
    }
  }
};


std::vector<edges_t> sample_edges(const filtered_directed_graph_t& graph, const flagser_parameters& params, const int drop_rate) {
  int min_dimension = 2;
  if (params.min_dimension >= 2) min_dimension = params.min_dimension;
  int max_dimension = 100; //hacky
  if (params.max_dimension >= 2) max_dimension = params.max_dimension;
  assert(min_dimension <= max_dimension);

  // generate flag complex
  directed_flag_complex_in_memory_t<NoneType> complex(graph, params.nb_threads, max_dimension);
  
  std::cout << "generated flag complex" << std::endl;

	std::vector<edges_t> edges;
  for (int dim = min_dimension-2; dim <= max_dimension-2; dim++) {
    //initialize threads
		std::vector<ads_computer_t> ads_enumerator(params.nb_threads, ads_computer_t{graph, drop_rate});
    //actually compute
    complex.for_each_cell(ads_enumerator, dim); 
    
    //sum results
    edges_t edges_in_dim{};
    for (const auto& thread : ads_enumerator) edges_in_dim.insert(edges_in_dim.end(), thread.edges.begin(), thread.edges.end());

    //debug
    /*
    std::cout << '[';
    for (auto edge : edges_in_dim) std::cout << '(' << edge.first << ',' << edge.second << "), ";
    std::cout << ']' ;
    */

    //std::cout << dim+2 << ": " << edges_in_dim << std::endl;

    // nads_in_dim == 0 => nads in higher dimensions must also be 0, thus we may abort.
    if (edges_in_dim.size() == 0) break;
    edges.push_back(edges_in_dim);
  }
	return edges;
}

int main(int argc, char** argv) {
	try {
		auto arguments = parse_arguments(argc, argv);

		auto positional_arguments = get_positional_arguments(arguments);
		auto named_arguments = get_named_arguments(arguments);
		auto params = flagser_parameters(named_arguments);
		named_arguments_t::const_iterator it;

		const char* input_filename = positional_arguments[0];

		filtered_directed_graph_t graph = read_filtered_directed_graph(input_filename, params);

		auto cell_count = sample_edges(graph, params, 1);
	} catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}
