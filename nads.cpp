#include <fstream>
#include <iostream>
#include <cstdint>

#define USE_CELLS_WITHOUT_DIMENSION

#include "include/argparser.h"
#include "include/parameters.h"
#include "include/usage/flagser-count.h"
#include "include/complex/directed_flag_complex_in_memory_computer.h" 


typedef std::vector<std::pair<int, vertex_index_t>> coboundaries_t;
typedef struct{} NoneType;

struct nads_computer_t {
  const filtered_directed_graph_t& graph;
  uint64_t nads = 0;
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
    // now calculate the number of almost-d-simplices
    for (size_t i = 0; i < coboundaries.size(); i++) {
      for (size_t j = i+1; j < coboundaries.size(); j++) {
        auto cb1 = coboundaries[i];
        auto cb2 = coboundaries[j];
        if (cb1.second != cb2.second){
          if (cb1.first <= cb2.first) nads += 1;
          if (cb1.first >= cb2.first) nads += 1;
        }
      }
    }
  }
};


std::vector<uint64_t> count_nads(const filtered_directed_graph_t& graph, const flagser_parameters& params) {
  int min_dimension = 1;
  if (params.min_dimension >= 1) min_dimension = params.min_dimension;
  int max_dimension = 100; //hacky
  if (params.max_dimension >= 1) max_dimension = params.max_dimension;
  assert(min_dimension <= max_dimension);

  // generate flag complex
  directed_flag_complex_in_memory_t<NoneType> complex(graph, params.nb_threads, max_dimension);
  
  std::cout << "generated flag complex" << std::endl;

	std::vector<uint64_t> nads(max_dimension - min_dimension + 1, 0);
  for (int dim = min_dimension-1; dim <= max_dimension; dim++) {
    //initialize threads
		std::vector<nads_computer_t> nads_counter(params.nb_threads, nads_computer_t{graph});
    //actually compute
    complex.for_each_cell(nads_counter, dim); 
    
    //sum results
    uint64_t nads_in_dim = 0;
    for (const auto& thread : nads_counter) nads_in_dim += thread.nads;
    
    std::cout << dim+1 << ": " << nads_in_dim << std::endl;
    nads.push_back(nads_in_dim);

    if (nads_in_dim == 0) {
      // => nads in higher dimensions must also be 0, thus we may abort.
      nads.resize(dim);
      break;
    }
  }
	return nads;
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

		auto cell_count = count_nads(graph, params);
	} catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}
